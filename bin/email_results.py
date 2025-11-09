#!/usr/bin/env python3

"""
email_results.py
Builds HTML summaries from a samplesheet.csv, optionally embeds an execution
timeline report and execution trace table, then sends them via email.

Expected CSV columns (no header):
0 id, 1 fastq_1, 2 fastq_2, 3 reads_after_trim, 4 kept, 5 fastp_html,
6 pct_mapped, 7 bam_name, 8 bam_size_mb, 9 counts
"""

import argparse
import csv
import html
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from typing import Any, List, Optional


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.2"


@dataclass
class TimelineDocument:
    """Parsed fragments of the Nextflow execution timeline HTML report."""

    head_markup: str
    body_markup: str
    source_path: Optional[str] = None


@dataclass
class TraceTable:
    """Represents the parsed Nextflow execution trace table."""

    headers: List[str]
    rows: List[List[str]]
    source_path: Optional[str] = None


class EmailConfig:
    """Configuration for email sending and report generation."""

    def __init__(
        self,
        email: str = "",
        email_on_fail: str = "",
        plaintext_email: bool = False,
        outdir: str = "",
        samplesheet: str = "",
        interactive: bool = False,
        status: str = "success",
        pipeline_name: str = "MetaAssembly",
        run_name: str = "",
        # New SMTP configuration options
        smtp_server: str = "smtp.gmail.com",
        smtp_port: int = 465,
        smtp_user: str = "",
        smtp_password: str = "",
        from_addr: str = "",
        smtp_security: str = "ssl",
        timeline: str = "",
        trace: str = "",
        use_mailx: bool = False,
    ):
        self.email = email
        self.email_on_fail = email_on_fail
        self.plaintext_email = plaintext_email
        self.outdir = outdir
        self.samplesheet = samplesheet
        self.interactive = interactive
        self.status = status
        self.pipeline_name = pipeline_name
        self.run_name = run_name
        # SMTP configuration
        self.smtp_server = smtp_server
        self.smtp_port = smtp_port
        self.smtp_user = smtp_user
        self.smtp_password = smtp_password
        self.from_addr = from_addr
        self.smtp_security = smtp_security
        self.timeline = timeline
        self.trace = trace
        self.use_mailx = use_mailx


class SampleRow:
    """Represents a single row from the samplesheet with parsed data."""

    def __init__(
        self,
        sample_id: str,
        fastq_1: str,
        fastq_2: str,
        reads_after_trim: int,
        read_after_trim_percent: float,
        read_after_decontamination: float,
        percent_mapped: float,
        bam_size_gb: float,
        assembled_count: int,
    ):
        self.sample_id = sample_id
        self.fastq_1 = fastq_1
        self.fastq_2 = fastq_2
        self.reads_after_trim = reads_after_trim
        self.read_after_trim_percent = read_after_trim_percent
        self.read_after_decontamination = read_after_decontamination
        self.percent_mapped = percent_mapped
        self.bam_size_gb = bam_size_gb
        self.assembled_count = assembled_count


def e(s: str) -> str:
    """
    HTML escape a string, safe for None values.

    Args:
        s: String to escape, can be None

    Returns:
        HTML-escaped string, or empty string if input was None
    """
    return html.escape(s or "")


def safe_int(x: Any) -> int:
    """
    Safely convert value to integer, returning 0 on failure.

    Args:
        x: Value to convert to integer

    Returns:
        Integer value, or 0 if conversion fails
    """
    try:
        return int(x)
    except (ValueError, TypeError):
        return 0


def safe_float(x: Any) -> float:
    """
    Safely convert value to float, returning 0.0 on failure.

    Args:
        x: Value to convert to float

    Returns:
        Float value, or 0.0 if conversion fails
    """
    try:
        return float(x)
    except (ValueError, TypeError):
        return 0.0


def fmt_int(n: int) -> str:
    """
    Format integer with thousands separators.

    Args:
        n: Integer to format

    Returns:
        Formatted string with commas as thousands separators
    """
    return f"{n:,}"


def fmt_float(n: float, nd: int = 3) -> str:
    """
    Format float with specified decimal places.

    Args:
        n: Float number to format
        nd: Number of decimal places (default: 3)

    Returns:
        Formatted string with specified precision
    """
    return f"{n:.{nd}f}"


def compute_max(rows: List[SampleRow], key: str, minimum: float = 1.0) -> float:
    """
    Compute maximum value for a given key across all sample rows.

    Args:
        rows: List of SampleRow objects
        key: Attribute name to get maximum value for
        minimum: Minimum value to return if no valid values found

    Returns:
        Maximum value found, or minimum if no valid values
    """
    vals = [getattr(r, key, 0) or 0 for r in rows]
    return max([minimum] + vals)


def _extract_tag_content(markup: str, tag: str) -> str:
    """
    Extract the inner HTML of the requested tag.

    Args:
        markup: Full HTML document as a string
        tag: Tag name (e.g., "head", "body")

    Returns:
        Inner HTML contained within the tag, or empty string if not found
    """
    pattern = rf"(?is)<{tag}\b[^>]*>(.*?)</{tag}>"
    match = re.search(pattern, markup)
    return match.group(1).strip() if match else ""


def render_timeline_section(timeline: TimelineDocument) -> str:
    """
    Build a toggleable section that renders the execution timeline HTML inline.

    Args:
        timeline: Parsed timeline fragments to embed

    Returns:
        HTML snippet containing the disclosure widget
    """
    if not timeline.body_markup.strip():
        return ""

    safe_suffix = ""
    if timeline.source_path:
        safe_suffix = re.sub(r"[^a-zA-Z0-9_-]", "-", os.path.basename(timeline.source_path))
    toggle_id = f"timeline-toggle-{safe_suffix or 'report'}"
    location = e(os.path.abspath(timeline.source_path)) if timeline.source_path else "provided execution timeline"
    caution = (
        "Note: many email clients block interactive content. "
        "If the chart does not appear, open this report in a web browser."
    )

    return f"""
  <section class="timeline-card">
    <input type="checkbox" id="{toggle_id}" class="timeline-toggle">
    <label for="{toggle_id}" class="timeline-label">Show execution timeline</label>
    <p class="meta timeline-meta">Source: <code>{location}</code></p>
    <div class="timeline-content">
      <div class="timeline-fragment">
{timeline.body_markup}
      </div>
    </div>
    <p class="timeline-note">{caution}</p>
  </section>
""".strip()


def load_timeline_html(timeline_path: str) -> Optional[TimelineDocument]:
    """
    Load the execution timeline HTML from disk, if provided.

    Args:
        timeline_path: Filesystem path to execution_timeline*.html

    Returns:
        TimelineDocument with extracted head/body markup, or None when unavailable
    """
    if not timeline_path:
        return None
    if not os.path.isfile(timeline_path):
        print(
            f"[email_results.py] WARNING: timeline not found: {timeline_path}",
            file=sys.stderr,
        )
        return None
    try:
        with open(timeline_path, "r", encoding="utf-8", errors="replace") as handle:
            raw_html = handle.read()
    except OSError as exc:
        print(
            f"[email_results.py] WARNING: failed to read timeline: {exc}",
            file=sys.stderr,
        )
        return None
    # Remove the doctype to avoid duplicates when embedding
    raw_html = re.sub(r"(?is)<!doctype[^>]*>", "", raw_html).strip()
    head_markup = _extract_tag_content(raw_html, "head")
    body_markup = _extract_tag_content(raw_html, "body") or raw_html
    return TimelineDocument(
        head_markup=head_markup,
        body_markup=body_markup,
        source_path=timeline_path,
    )


def load_trace_table(trace_path: str) -> Optional[TraceTable]:
    """
    Load the Nextflow execution trace table from disk.

    Args:
        trace_path: Filesystem path to execution_trace*.txt

    Returns:
        TraceTable instance, or None when unavailable
    """
    if not trace_path:
        return None
    if not os.path.isfile(trace_path):
        print(
            f"[email_results.py] WARNING: trace not found: {trace_path}",
            file=sys.stderr,
        )
        return None
    try:
        with open(trace_path, newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.reader(handle, delimiter="\t")
            rows = [row for row in reader if row]
    except OSError as exc:
        print(
            f"[email_results.py] WARNING: failed to read trace: {exc}",
            file=sys.stderr,
        )
        return None

    if not rows:
        return None
    headers, *data_rows = rows
    return TraceTable(headers=headers, rows=data_rows, source_path=trace_path)


def render_trace_section(trace_table: TraceTable) -> str:
    """
    Convert the trace table into an HTML section with styled status badges.

    Args:
        trace_table: Parsed trace table data

    Returns:
        HTML snippet containing the trace table wrapped in a details block
    """
    if not trace_table.rows:
        return ""

    status_idx = -1
    for idx, header in enumerate(trace_table.headers):
        if header.lower() == "status":
            status_idx = idx
            break

    table_rows = []
    for row in trace_table.rows:
        cells = []
        for idx, value in enumerate(row):
            text = e(value)
            if idx == status_idx:
                normalized = value.strip().lower()
                if normalized in {"completed", "complete"}:
                    pill_class = "status-pill status-complete"
                else:
                    pill_class = "status-pill status-failed"
                cells.append(f'<td><span class="{pill_class}">{text}</span></td>')
            else:
                cells.append(f"<td>{text}</td>")
        table_rows.append("<tr>" + "".join(cells) + "</tr>")

    headers_html = "".join(f"<th>{e(h)}</th>" for h in trace_table.headers)
    location = (
        e(os.path.abspath(trace_table.source_path))
        if trace_table.source_path
        else "provided execution trace"
    )

    return f"""
  <section class="tasks-card">
    <details>
      <summary>Show Tasks</summary>
      <p class="meta">Source: <code>{location}</code></p>
      <div class="trace-table-wrapper">
        <table class="trace-table">
          <thead><tr>{headers_html}</tr></thead>
          <tbody>
            {"".join(table_rows)}
          </tbody>
        </table>
      </div>
    </details>
  </section>
""".strip()


def to_plaintext(html_str: str) -> str:
    """
    Convert HTML to plain text by stripping tags and normalizing whitespace.

    Args:
        html_str: HTML string to convert

    Returns:
        Plain text version of the HTML content
    """
    txt = re.sub(r"(?is)<style.*?</style>", "", html_str)
    txt = re.sub(r"(?is)<script.*?</script>", "", txt)
    txt = txt.replace("<br>", "\n").replace("<br/>", "\n").replace("<br />", "\n")
    txt = re.sub(r"(?i)</p>", "\n\n", txt)
    txt = re.sub(r"(?i)</h[1-6]>", "\n", txt)
    txt = re.sub(r"(?s)<[^>]+>", "", txt)
    txt = re.sub(r"[ \t]+\n", "\n", txt)
    return txt.strip()


def load_rows(csv_path: str) -> List[SampleRow]:
    """
    Load and parse samplesheet CSV into SampleRow objects.

    Args:
        csv_path: Path to the samplesheet CSV file

    Returns:
        List of SampleRow objects with parsed data

    Raises:
        FileNotFoundError: If CSV file doesn't exist
        PermissionError: If CSV file cannot be read
    """
    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        for cols in reader:
            if not cols or not any(cols):
                continue
            # Pad missing columns to ensure we have at least 10
            cols = (cols + [""] * 10)[:10]

            # Parse with safe fallbacks for conversion errors
            row = SampleRow(
                sample_id=cols[0],
                fastq_1=cols[1],
                fastq_2=cols[2],
                reads_after_trim=safe_int(cols[3]),
                read_after_trim_percent=safe_float(cols[4]),
                read_after_decontamination=safe_float(cols[5]),
                percent_mapped=safe_float(cols[6]),
                bam_size_gb=safe_float(cols[7]),
                assembled_count=safe_int(cols[8]),
            )
            rows.append(row)
    return rows


def build_email_html(
    rows: List[SampleRow],
    subject: str,
    outdir: str,
    timeline: Optional[TimelineDocument] = None,
    trace: Optional[TraceTable] = None,
) -> str:
    """
    Build email-safe HTML table report with progress bars.

    Args:
        rows: List of SampleRow objects to include in report
        subject: Email subject for the report title
        outdir: Output directory path for reference in report
        timeline: Optional TimelineDocument to embed
        trace: Optional TraceTable to include

    Returns:
        Complete HTML document as string
    """
    rows_html = []
    for r in rows:
        reads_pct =  r.read_after_trim_percent
        reads_after_deacon_percent = r.read_after_decontamination
        reads_after_deacon = int(r.reads_after_trim * reads_after_deacon_percent / 100)

        rows_html.append(
            f"""
      <tr>
        <td>{e(r.sample_id)}</td>
        <td>{e(r.fastq_1)}</td>
        <td>{e(r.fastq_2)}</td>
        <td class="num">
          {fmt_int(r.reads_after_trim)} ({fmt_float(r.read_after_trim_percent)}%)
          <div class="bar"><div class="fill" style="width:{reads_pct:.1f}%"></div></div>
        </td>
        <td class="num">
          {fmt_int(reads_after_deacon)} ({fmt_float(reads_after_deacon_percent)}%)
          <div class="bar"><div class="fill" style="width:{reads_after_deacon_percent:.1f}%"></div></div>
        </td>
        <td class="num">{r.percent_mapped:.2f}</td>
        <td class="num">{fmt_float(r.bam_size_gb)}GB</td>
        <td class="num">{fmt_int(r.assembled_count)}</td>
      </tr>
    """.rstrip()
        )

    sentry = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{e(subject)}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
{timeline.head_markup if timeline else ""}
<style>
  body{{font-family:Arial,Helvetica,sans-serif;line-height:1.45;color:#222;margin:16px;}}
  h1{{margin:0 0 6px 0;font-size:20px}}
  p.meta{{color:#666;font-size:12px;margin:0 0 12px 0}}
  table{{border-collapse:collapse;width:100%;max-width:1200px}}
  th,td{{padding:8px;border-bottom:1px solid #eee;vertical-align:top}}
  thead th{{background:#f7f7f7;border-bottom:1px solid #ddd;text-align:left}}
  td.num{{text-align:right;white-space:nowrap}}
  .bar{{width:140px;height:10px;background:#eee;border-radius:6px;margin-top:4px;float:right}}
  .fill{{height:10px;background:#4b9ae8;border-radius:6px}}
  a{{color:#1e6bd6;text-decoration:none}}
  code{{background:#f3f3f3;padding:2px 4px;border-radius:4px}}
  .timeline-card{{margin-top:24px;border:1px solid #e6e6e6;border-radius:10px;padding:14px;background:#fafafa}}
  .timeline-toggle{{position:absolute;opacity:0}}
  .timeline-label{{display:inline-block;font-weight:600;cursor:pointer}}
  .timeline-card .timeline-content{{display:none;margin-top:14px;overflow:auto}}
  .timeline-card .timeline-toggle:checked ~ .timeline-content{{display:block}}
  .timeline-meta{{font-size:12px;margin:8px 0 0 0;color:#555}}
  .timeline-note{{font-size:12px;color:#777;margin-top:12px}}
  .timeline-fragment{{background:#fff;padding:12px;border:1px solid #ddd;border-radius:8px}}
  .tasks-card{{margin-top:24px}}
  .tasks-card details{{border:1px solid #e6e6e6;border-radius:10px;padding:12px;background:#fafafa}}
  .tasks-card summary{{font-weight:600;cursor:pointer;outline:none}}
  .trace-table-wrapper{{overflow:auto;margin-top:12px}}
  .trace-table{{width:100%;border-collapse:collapse;font-size:13px}}
  .trace-table th,.trace-table td{{padding:6px 8px;border-bottom:1px solid #eee;text-align:left;white-space:nowrap}}
  .status-pill{{display:inline-block;padding:2px 10px;border-radius:999px;font-size:12px;font-weight:600;color:#fff}}
  .status-pill.status-complete{{background:#2d9d4d}}
  .status-pill.status-failed{{background:#c0392b}}
</style>
</head>
<body>
  <h1>MetaAssembly – Run summary</h1>
  <p class="meta">Generated: {e(sentry)} &middot; Reports directory: <code>{e(os.path.join(outdir, "reports"))}</code></p>

  <table>
    <thead>
      <tr>
        <th>Sample</th>
        <th>FASTQ 1</th>
        <th>FASTQ 2</th>
        <th>After fastp</th>
        <th>After deacon</th>
        <th>Mapped %</th>
        <th>BAM size</th>
        <th>Assembled</th>
      </tr>
    </thead>
    <tbody>
      {"".join(rows_html)}
    </tbody>
  </table>

  {render_trace_section(trace) if trace else ""}

  {render_timeline_section(timeline) if timeline else ""}

  <p class="meta">This email is self-contained (no JS, no external CSS).</p>
</body>
</html>"""
    return html_doc


def build_interactive_html(
    rows: List[SampleRow],
    subject: str,
    outdir: str,
    timeline: Optional[TimelineDocument] = None,
    trace: Optional[TraceTable] = None,
) -> str:
    """
    Build interactive HTML report with expandable cards for each sample.

    Args:
        rows: List of SampleRow objects to include in report
        subject: Report title
        outdir: Output directory path for reference in report
        timeline: Optional TimelineDocument to embed
        trace: Optional TraceTable to include

    Returns:
        Complete HTML document as string
    """
    panels = []
    for r in rows:
        reads_after_deacon_percent = r.read_after_decontamination
        reads_after_deacon = int(r.reads_after_trim * reads_after_deacon_percent)

        panels.append(
            f"""
<details class="card">
  <summary><b>{e(r.sample_id)}</b> — {e(r.fastq_1)}{(" / " + e(r.fastq_2)) if r.fastq_2 else ""}</summary>
  <div class="grid">
    <div><div class="k">After fastp</div><div class="v">{fmt_int(r.reads_after_trim)}</div></div>
    <div><div class="k">After deacon</div><div class="v">{fmt_int(reads_after_deacon)}</div></div>
    <div><div class="k">Mapped %</div><div class="v">{r.percent_mapped:.2f}</div></div>
    <div><div class="k">BAM size</div><div class="v">{fmt_float(r.bam_size_gb)}</div></div>
    <div><div class="k">Assembled</div><div class="v">{fmt_int(r.assembled_count)}</div></div>
  </div>
</details>
""".rstrip()
        )

    sentry = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{e(subject)}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
{timeline.head_markup if timeline else ""}
<style>
  body{{font-family:Inter,Arial,Helvetica,sans-serif;line-height:1.5;color:#222;margin:24px}}
  h1{{margin:0 0 8px 0;font-size:22px}}
  p.meta{{color:#666;font-size:13px;margin:0 0 16px 0}}
  .cards{{display:grid;grid-template-columns:1fr;gap:12px}}
  @media(min-width:900px){{.cards{{grid-template-columns:1fr 1fr}}}}
  .card{{border:1px solid #e9e9e9;border-radius:12px;padding:10px;background:#fff;box-shadow:0 1px 3px rgba(0,0,0,.05)}}
  summary{{cursor:pointer;outline:none}}
  .grid{{display:grid;grid-template-columns:repeat(2, minmax(0,1fr));gap:10px;margin-top:10px}}
  .k{{color:#555;font-size:12px}}
  .v{{font-weight:600}}
  a{{color:#1e6bd6;text-decoration:none}}
  code{{background:#f3f3f3;padding:2px 4px;border-radius:4px}}
  .timeline-card{{margin-top:32px;border:1px solid #e6e6e6;border-radius:12px;padding:16px;background:#fafafa}}
  .timeline-toggle{{position:absolute;opacity:0}}
  .timeline-label{{font-weight:600;font-size:14px;cursor:pointer}}
  .timeline-card .timeline-content{{display:none;margin-top:16px;overflow:auto}}
  .timeline-card .timeline-toggle:checked ~ .timeline-content{{display:block}}
  .timeline-meta{{font-size:12px;margin:10px 0 0 0;color:#555}}
  .timeline-note{{font-size:12px;color:#777;margin-top:10px}}
  .timeline-fragment{{background:#fff;padding:14px;border:1px solid #ddd;border-radius:10px}}
  .tasks-card{{margin-top:28px}}
  .tasks-card details{{border:1px solid #e6e6e6;border-radius:12px;padding:14px;background:#fafafa}}
  .tasks-card summary{{font-weight:600;cursor:pointer}}
  .trace-table-wrapper{{overflow:auto;margin-top:14px}}
  .trace-table{{width:100%;border-collapse:collapse;font-size:13px}}
  .trace-table th,.trace-table td{{padding:6px 10px;border-bottom:1px solid #eee;text-align:left;white-space:nowrap}}
  .status-pill{{display:inline-block;padding:3px 12px;border-radius:999px;font-size:12px;font-weight:600;color:#fff}}
  .status-pill.status-complete{{background:#2d9d4d}}
  .status-pill.status-failed{{background:#c0392b}}
</style>
</head>
<body>
  <h1>MetaAssembly – Interactive Summary</h1>
  <p class="meta">Generated: {e(sentry)} &middot; Reports at <code>{e(os.path.join(outdir, "reports"))}</code>. Expand a card to see per-sample details.</p>
  <div class="cards">
    {"".join(panels)}
  </div>
  {render_trace_section(trace) if trace else ""}
  {render_timeline_section(timeline) if timeline else ""}
</body>
</html>"""
    return doc


def send_via_smtplib(
    to_addr: str,
    subject: str,
    html_body: str,
    smtp_server: str = "smtp.gmail.com",
    smtp_port: int = 465,
    smtp_user: str = "",
    smtp_password: str = "",
    from_addr: str = "",
    smtp_security: str = "ssl",
) -> None:
    """
    Send email using Python's smtplib with configurable security.

    Args:
        to_addr: Recipient email address
        subject: Email subject
        html_body: HTML email content
        smtp_server: SMTP server hostname
        smtp_port: SMTP server port
        smtp_user: SMTP username/email
        smtp_password: SMTP password or app password
        from_addr: From email address
        smtp_security: Either "ssl" or "tls"

    Raises:
        smtplib.SMTPException: If SMTP communication fails
        Exception: For other connection issues
    """
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    # Create message
    msg = MIMEMultipart()
    msg["From"] = from_addr or smtp_user
    msg["To"] = to_addr
    msg["Subject"] = subject

    # Add HTML body
    msg.attach(MIMEText(html_body, "html", "utf-8"))

    security = (smtp_security or "ssl").lower()
    if security not in {"ssl", "tls"}:
        raise ValueError(f"Unsupported SMTP security mode: {smtp_security}")

    if security == "ssl":
        connector = smtplib.SMTP_SSL
    else:
        connector = smtplib.SMTP

    with connector(smtp_server, smtp_port) as server:
        if security == "tls":
            server.starttls()
        server.login(smtp_user, smtp_password)
        server.send_message(msg)


def send_via_mailx(
    to_addr: str,
    subject: str,
    body: str,
    from_addr: str = "",
    content_type: str = "plain",
) -> None:
    """
    Send email using the system mailx command.

    Args:
        to_addr: Recipient email address
        subject: Email subject line
        body: Email body sent through STDIN
        from_addr: Optional From address passed with -r
        content_type: Either "plain" or "html" to set MIME headers
    """
    cmd = ["mailx"]
    if from_addr:
        cmd.extend(["-r", from_addr])
    if content_type == "html":
        cmd.extend(
            [
                "-S",
                "mime=1",
                "-S",
                "content-type=text/html",
                "-S",
                "charset=UTF-8",
            ]
        )
    cmd.extend(["-s", subject, to_addr])

    subprocess.run(cmd, input=body, text=True, check=True)


def send_email(
    to_addr: str,
    subject: str,
    html_body: str,
    plaintext: bool = False,
    smtp_server: str = "smtp.gmail.com",
    smtp_port: int = 465,
    smtp_user: str = "",
    smtp_password: str = "",
    from_addr: str = "",
    smtp_security: str = "ssl",
    use_mailx: bool = False,
) -> None:
    """
    Send email using Python's smtplib with configurable security.

    Args:
        to_addr: Recipient email address
        subject: Email subject
        html_body: HTML email content
        plaintext: Whether to force plaintext email
        smtp_server: SMTP server hostname
        smtp_port: SMTP server port
        smtp_user: SMTP username/email
        smtp_password: SMTP password or app password
        from_addr: From email address
        smtp_security: Either "ssl" or "tls"
        use_mailx: When True, route through the system mailx command

    Raises:
        smtplib.SMTPException: If SMTP communication fails
        Exception: For other connection issues
    """
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    # Choose content type
    if plaintext:
        email_body = to_plaintext(html_body)
        content_type = "plain"
    else:
        email_body = html_body
        content_type = "html"

    if use_mailx:
        send_via_mailx(
            to_addr=to_addr,
            subject=subject,
            body=email_body,
            from_addr=from_addr or smtp_user,
            content_type=content_type,
        )
        return

    # Create message
    msg = MIMEMultipart()
    msg["From"] = from_addr or smtp_user
    msg["To"] = to_addr
    msg["Subject"] = subject

    # Add body
    msg.attach(MIMEText(email_body, content_type, "utf-8"))

    security = (smtp_security or "ssl").lower()
    if security not in {"ssl", "tls"}:
        raise ValueError(f"Unsupported SMTP security mode: {smtp_security}")

    connector = smtplib.SMTP_SSL if security == "ssl" else smtplib.SMTP

    with connector(smtp_server, smtp_port) as server:
        if security == "tls":
            server.starttls()
        server.login(smtp_user, smtp_password)
        server.send_message(msg)


def parse_args() -> EmailConfig:
    """
    Parse command line arguments and return EmailConfig object.

    Returns:
        EmailConfig object with all parsed arguments

    Raises:
        SystemExit: If argument parsing fails
    """
    ap = argparse.ArgumentParser(description="Send MetaAssembly results by email.")

    # Email configuration arguments
    ap.add_argument(
        "--email",
        dest="email",
        default="",
        help="Recipient email address for successful runs",
    )
    ap.add_argument(
        "--email-on-fail",
        dest="email_on_fail",
        default="",
        help="Recipient email used only on failure when --email is empty",
    )
    ap.add_argument(
        "--plaintext-email",
        action="store_true",
        help="Send plaintext instead of HTML email",
    )

    # SMTP configuration arguments
    ap.add_argument(
        "--smtp-server",
        dest="smtp_server",
        default="smtp.gmail.com",
        help="SMTP server hostname (default: smtp.gmail.com)",
    )
    ap.add_argument(
        "--smtp-port",
        dest="smtp_port",
        type=int,
        default=465,
        help="SMTP server port (default: 465 for SSL)",
    )
    ap.add_argument(
        "--smtp-user", dest="smtp_user", default="", help="SMTP username/email address"
    )
    ap.add_argument(
        "--smtp-password",
        dest="smtp_password",
        default="",
        help="SMTP password or app password",
    )
    ap.add_argument(
        "--from-addr",
        dest="from_addr",
        default="",
        help="From email address (defaults to smtp-user if not provided)",
    )
    ap.add_argument(
        "--smtp-security",
        dest="smtp_security",
        choices=["ssl", "tls"],
        default="ssl",
        help="SMTP security mode (default: ssl)",
    )
    ap.add_argument(
        "--use-mailx",
        dest="use_mailx",
        action="store_true",
        help="Send email via the system mailx CLI instead of smtplib",
    )

    # File and directory arguments
    ap.add_argument(
        "--outdir", required=True, help="Output directory (to store copies of reports)"
    )
    ap.add_argument("--samplesheet", required=True, help="Path to samplesheet.csv")
    ap.add_argument(
        "--timeline",
        dest="timeline",
        default="",
        help="Path to Nextflow execution timeline HTML to embed",
    )
    ap.add_argument(
        "--trace",
        dest="trace",
        default="",
        help="Path to Nextflow execution trace table (.txt) to embed",
    )

    # Report generation arguments
    ap.add_argument(
        "--interactive",
        action="store_true",
        help="Also write an interactive (browser) report",
    )
    ap.add_argument(
        "--status",
        choices=["success", "failed"],
        default="success",
        help="Overall pipeline status (affects email_on_fail routing)",
    )
    ap.add_argument(
        "--pipeline-name",
        default="MetaAssembly",
        help="Pipeline name used in email subject/title",
    )
    ap.add_argument("--run-name", default="", help="Run name for email subject/title")

    args = ap.parse_args()

    return EmailConfig(
        email=args.email,
        email_on_fail=args.email_on_fail,
        plaintext_email=args.plaintext_email,
        outdir=args.outdir,
        samplesheet=args.samplesheet,
        interactive=args.interactive,
        status=args.status,
        pipeline_name=args.pipeline_name,
        run_name=args.run_name,
        smtp_server=args.smtp_server,
        smtp_port=args.smtp_port,
        smtp_user=args.smtp_user,
        smtp_password=args.smtp_password,
        from_addr=args.from_addr,
        smtp_security=args.smtp_security,
        timeline=args.timeline,
        trace=args.trace,
        use_mailx=args.use_mailx,
    )


def main() -> int:
    """
    Main execution function for email results processing.

    Returns:
        Exit code (0 for success, non-zero for errors)
    """
    # Parse command line arguments
    config = parse_args()

    # Determine recipient based on status
    to_addr = config.email
    if (not to_addr) and config.email_on_fail and (config.status == "failed"):
        to_addr = config.email_on_fail

    # Build email subject
    status_label = "Successful" if config.status == "success" else "FAILED"
    subject = f"[{config.pipeline_name}] {status_label}: {config.run_name or 'run'}"

    # Validate and load samplesheet
    if not os.path.isfile(config.samplesheet):
        print(
            f"[email_results.py] ERROR: samplesheet not found: {config.samplesheet}",
            file=sys.stderr,
        )
        return 2

    rows = load_rows(config.samplesheet)

    # Ensure output directories exist
    reports_dir = os.path.join(config.outdir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    # Optionally load execution timeline HTML for embedding
    timeline_doc = load_timeline_html(config.timeline)
    trace_table = load_trace_table(config.trace)

    # Build HTML reports
    email_html = build_email_html(
        rows, subject, config.outdir, timeline=timeline_doc, trace=trace_table
    )
    email_html_path = os.path.join(reports_dir, "metaassembly_email_summary.html")
    with open(email_html_path, "w", encoding="utf-8") as f:
        f.write(email_html)

    # Build interactive report if requested
    if config.interactive:
        interactive_html = build_interactive_html(
            rows,
            subject,
            config.outdir,
            timeline=timeline_doc,
            trace=trace_table,
        )
        interactive_html_path = os.path.join(
            reports_dir, "metaassembly_interactive_report.html"
        )
        with open(interactive_html_path, "w", encoding="utf-8") as f:
            f.write(interactive_html)

    if to_addr:
        try:
            send_email(
                to_addr=to_addr,
                subject=subject,
                html_body=email_html,
                plaintext=config.plaintext_email,
                smtp_server=config.smtp_server,
                smtp_port=config.smtp_port,
                smtp_user=config.smtp_user,
                smtp_password=config.smtp_password,
                from_addr=config.from_addr,
                smtp_security=config.smtp_security,
                use_mailx=config.use_mailx,
            )
        except Exception as ex:
            print(
                f"[email_results.py] ERROR: failed to send email: {ex}", file=sys.stderr
            )
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
