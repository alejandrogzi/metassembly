/**
 * Parse samplesheet.csv (no header) into a list of maps.
 * Expected columns:
 * 0 id, 1 fastq_1, 2 fastq_2, 3 reads_after_trim, 4 kept, 5 fastp_html,
 * 6 pct, 7 bam_name, 8 bam_size_mb, 9 counts
 */
def parseSamplesheet(File f) {
    def rows = []
    f.eachLine { line ->
        def s = line?.trim()
        if (!s) return
        def cols = s.split(/,(?=(?:[^"]*"[^"]*")*[^"]*$)/, -1) // robust-ish CSV split
        rows << [
            id            : cols[0],
            fastq_1       : cols[1],
            fastq_2       : cols[2],
            reads_after   : (cols[3] ?: "0") as long,
            kept          : (cols[4] ?: "0") as long,
            fastp_html    : cols[5],
            pct_mapped    : (cols[6] ?: "0") as double,
            bam_name      : cols[7],
            bam_size_mb   : (cols[8] ?: "0") as double,
            counts        : (cols[9] ?: "0") as long
        ]
    }
    return rows
}

/**
 * Email-safe HTML (no JS, no external resources).
 * Uses a styled table and inline "bar" cells using simple DIV widths.
 */
def buildEmailHtml(List rows) {
    def maxReads = Math.max(1L, rows.collect{ it.reads_after }.max() ?: 1L)
    def maxCounts = Math.max(1L, rows.collect{ it.counts }.max() ?: 1L)
    def maxBam = Math.max(0.0001, rows.collect{ it.bam_size_mb }.max() ?: 0.0001)

    def esc = { s -> (s ?: '').replace('&','&amp;').replace('<','&lt;').replace('>','&gt;') }

    def bar = { pct ->
        def w = Math.max(1, Math.min(100, (int)Math.round(pct)))
        return "<div style='background:#eee;width:100%;height:10px;border-radius:6px;'><div style='height:10px;border-radius:6px;width:${w}%;background:#4b9ae8;'></div></div>"
    }

    def html = new StringBuilder()
    html << """
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MetaAssembly – Summary</title>
</head>
<body style="font-family:Arial,Helvetica,sans-serif; line-height:1.4; color:#222;">
  <h2 style="margin:0 0 8px 0;">MetaAssembly – Run summary</h2>
  <p style="margin:0 0 16px 0;">Below is a compact overview of your samples. Links to fastp reports are included per sample.</p>
  <table cellpadding="6" cellspacing="0" border="0" style="border-collapse:collapse;width:100%;max-width:1000px;">
    <thead>
      <tr style="background:#f5f5f5;border-bottom:1px solid #ddd;">
        <th align="left">Sample</th>
        <th align="left">FASTQ 1</th>
        <th align="left">FASTQ 2</th>
        <th align="right">Reads After Trim</th>
        <th align="right">Kept</th>
        <th align="right">Mapped %</th>
        <th align="left">BAM</th>
        <th align="right">BAM MB</th>
        <th align="right">Counts</th>
        <th align="left">fastp</th>
      </tr>
    </thead>
    <tbody>
"""
    rows.each { r ->
        def readsPct = r.reads_after ? (100.0 * r.reads_after / maxReads) : 0
        def countsPct = r.counts ? (100.0 * r.counts / maxCounts) : 0
        def bamPct = r.bam_size_mb ? (100.0 * r.bam_size_mb / maxBam) : 0

        html << """
      <tr style="border-bottom:1px solid #eee;">
        <td>${esc(r.id)}</td>
        <td>${esc(r.fastq_1)}</td>
        <td>${esc(r.fastq_2)}</td>
        <td align="right">
          ${String.format("%,d", r.reads_after)}
          <div style="margin-top:4px;">${bar(readsPct)}</div>
        </td>
        <td align="right">${String.format("%,d", r.kept)}</td>
        <td align="right">${String.format(Locale.US, "%.2f", r.pct_mapped)}</td>
        <td>${esc(r.bam_name)}</td>
        <td align="right">
          ${String.format(Locale.US, "%.3f", r.bam_size_mb)}
          <div style="margin-top:4px;">${bar(bamPct)}</div>
        </td>
        <td align="right">
          ${String.format("%,d", r.counts)}
          <div style="margin-top:4px;">${bar(countsPct)}</div>
        </td>
        <td>${r.fastp_html ? "<a href='${esc(r.fastp_html)}'>fastp</a>" : ""}</td>
      </tr>
"""
    }

    html << """
    </tbody>
  </table>
  <p style="margin-top:18px;font-size:12px;color:#555;">Tip: this version avoids JavaScript so it renders well in most email clients. A richer interactive report is attached/linked.</p>
</body>
</html>
"""
    return html.toString()
}

/**
 * Interactive browser report with Plotly (bars for mapped %, BAM size, counts).
 * This is meant to be written to disk and opened in a browser (most email clients block JS).
 */
def buildInteractiveHtml(List rows, String title='MetaAssembly – Interactive Report') {
    // Build data arrays
    def ids         = rows.collect { it.id }
    def pctMapped   = rows.collect { it.pct_mapped }
    def bamSizes    = rows.collect { it.bam_size_mb }
    def counts      = rows.collect { it.counts }
    def readsAfter  = rows.collect { it.reads_after }

    // Escape JSON safely
    def json = new groovy.json.JsonBuilder([
        ids: ids, pctMapped: pctMapped, bamSizes: bamSizes, counts: counts, readsAfter: readsAfter
    ]).toString()

    return """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>${title}</title>
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <!-- Plotly via CDN (fine for browsers, not for email) -->
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body{font-family:Inter,system-ui,Segoe UI,Arial,sans-serif;margin:24px;color:#222;}
    h1{margin:0 0 8px 0;font-size:22px;}
    .grid{display:grid;grid-template-columns:1fr;gap:16px;}
    @media (min-width: 900px){ .grid{grid-template-columns:1fr 1fr;}}
    .card{border:1px solid #eee;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.05);padding:16px;}
    .muted{color:#666;font-size:14px;margin:0 0 12px 0;}
    a{color:#1e6bd6;text-decoration:none;}
  </style>
</head>
<body>
  <h1>${title}</h1>
  <p class="muted">Open this file in a browser for interactivity. The email version is simplified for compatibility.</p>
  <div class="grid">
    <div class="card"><div id="pct"></div></div>
    <div class="card"><div id="bam"></div></div>
    <div class="card"><div id="counts"></div></div>
    <div class="card"><div id="reads"></div></div>
  </div>
  <script>
    const data = ${json};
    function bar(id, x, y, xtitle, ytitle){
      Plotly.newPlot(id, [{
        x: x, y: y, type: 'bar', hovertemplate: '%{x}: %{y}<extra></extra>'
      }], {
        margin:{l:50,r:10,t:10,b:80},
        xaxis:{title: xtitle, automargin:true},
        yaxis:{title: ytitle},
        bargap: 0.25,
        responsive: true
      });
    }
    bar('pct',   data.ids, data.pctMapped, 'Sample', 'Mapped %');
    bar('bam',   data.ids, data.bamSizes,  'Sample', 'BAM size (MB)');
    bar('counts',data.ids, data.counts,    'Sample', 'Counts');
    bar('reads', data.ids, data.readsAfter,'Sample', 'Reads after trim');
  </script>
</body>
</html>
"""
}

/**
 * Build both HTMLs and write to disk.
 * Returns a map with { email_html_str, email_html_file, interactive_html_file }.
 */
def metassemblySummary(File samplesheetFile, String outdir) {
    def rows = parseSamplesheet(samplesheetFile)
    def reportsDir = new File("${outdir}/reports"); reportsDir.mkdirs()

    def emailHtml = buildEmailHtml(rows)
    def emailFile = new File(reportsDir, "metaassembly_email_summary.html")
    emailFile.text = emailHtml

    def interactiveHtml = buildInteractiveHtml(rows)
    def interactiveFile = new File(reportsDir, "metaassembly_interactive_report.html")
    interactiveFile.text = interactiveHtml

    return [
      email_html_str       : emailHtml,
      email_html_file      : emailFile.getAbsolutePath(),
      interactive_html_file: interactiveFile.getAbsolutePath()
    ]
}


//
// Construct and send completion email (template-free, HTML or plaintext)
//
def completionEmail(
    String email,                 // primary recipient (send on success/fail)
    String email_on_fail,         // alt recipient used only on failure when 'email' is empty
    boolean plaintext_email,      // force plaintext body
    String outdir,                // output root to store a copy of the sent report
    boolean monochrome_logs=true, // for logging colours only
    String html_body=''           // already-built HTML body (eg. from metassemblySummary)
) {

    // Subject – keep same style as your original
    def subject = "[${workflow.manifest.name}] Successful: ${workflow.runName}"
    if (!workflow.success) subject = "[${workflow.manifest.name}] FAILED: ${workflow.runName}"

    // Decide who to email
    def email_address = email
    if (!email && email_on_fail && !workflow.success) {
        email_address = email_on_fail
    }

    // If caller didn't give HTML, make a tiny default shell
    def htmlShell = { inner ->
        return """<!doctype html>
<html><head><meta charset="utf-8"><title>${subject}</title></head>
<body style="font-family:Arial,Helvetica,sans-serif;line-height:1.45;color:#222;">
${inner ?: "<p>No summary available.</p>"}
<hr style="margin-top:20px;border:none;border-top:1px solid #eee;">
<p style="font-size:12px;color:#666;margin:8px 0 0 0;">
Run: <b>${workflow.runName}</b> &middot; ${workflow.manifest.name}<br>
Started: ${workflow.start} &middot; Completed: ${workflow.complete}<br>
Nextflow: ${workflow.nextflow.version} (${workflow.nextflow.build})
</p>
</body></html>"""
    }

    def email_html = htmlShell(html_body)

    // Crude plaintext from HTML (works fine for “just send me text”)
    def email_txt = email_html
        .replaceAll("(?s)<style.*?</style>", "")
        .replaceAll("(?s)<script.*?</script>", "")
        .replaceAll("<br\\s*/?>", "\n")
        .replaceAll("</p>", "\n\n")
        .replaceAll("</h\\d>", "\n")
        .replaceAll("<[^>]+>", "")
        .replaceAll("[ \\t]+\\n", "\n")
        .trim()

    // Write copies alongside other pipeline info
    def infoDir = new File("${outdir}/pipeline_info"); infoDir.mkdirs()
    new File(infoDir, "pipeline_report.html").text = email_html
    new File(infoDir, "pipeline_report.txt").text  = email_txt

    // Nothing to send? bail quietly
    if (!email_address) return

    def colors = logColours(monochrome_logs) as Map

    // Try HTML via sendmail (headers + body)
    try {
        if (plaintext_email) throw new org.codehaus.groovy.GroovyException('Force plaintext')

        def msg = new StringBuilder()
        msg << "To: ${email_address}\n"
        msg << "Subject: ${subject}\n"
        msg << "MIME-Version: 1.0\n"
        msg << "Content-Type: text/html; charset=UTF-8\n"
        msg << "Content-Transfer-Encoding: 8bit\n"
        msg << "\n"
        msg << email_html

        ['sendmail','-t'].execute() << msg.toString()
        log.info("-${colors.purple}[${workflow.manifest.name}]${colors.green} Sent summary e-mail to ${email_address} (sendmail HTML)-")
        return
    }
    catch (Exception msg) {
        log.debug(msg.toString())
        log.debug("Falling back to plaintext or 'mail'")
    }

    // Fallback: try HTML via 'mail' (many distros support --content-type)
    try {
        if (!plaintext_email) {
            def mail_cmd = ['mail','-s', subject, '--content-type=text/html', email_address]
            mail_cmd.execute() << email_html
            log.info("-${colors.purple}[${workflow.manifest.name}]${colors.green} Sent summary e-mail to ${email_address} (mail HTML)-")
            return
        }
    } catch (Exception msg) {
        log.debug(msg.toString())
        log.debug("Falling back to plaintext")
    }

    // Final fallback: plaintext via 'mail'
    try {
        def mail_cmd = ['mail','-s', subject, email_address]
        mail_cmd.execute() << email_txt
        log.info("-${colors.purple}[${workflow.manifest.name}]${colors.green} Sent summary e-mail to ${email_address} (mail TXT)-")
    } catch (Exception msg) {
        log.error("Failed to send completion email to ${email_address}: ${msg}")
    }
}
