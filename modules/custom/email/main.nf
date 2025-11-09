process EMAIL_RESULTS {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    val email
    val email_on_fail
    val plaintext_email
    val outdir
    val use_mailx
    path samplesheet

    script:
    // INFO: if use_mailx, all smpt options are ignored
    if (use_mailx) {
        """
        email_results.py \\
            --email ${email} \\
            --email-on-fail ${email_on_fail} \\
            --outdir ${outdir} \\
            --samplesheet ${samplesheet} \\
            --use-mailx
        """
    } else {
        """
        email_results.py \\
            --email ${email} \\
            --email-on-fail ${email_on_fail} \\
            --outdir ${outdir} \\
            --samplesheet ${samplesheet} \\
            --smtp-server ${params.smtp_server} \\
            --smtp-port ${params.smtp_port} \\
            --smtp-user ${params.smtp_user} \\
            --smtp-password ${params.smtp_password} \\
            --smtp-security ${params.smtp_security}
        """
    }
}
