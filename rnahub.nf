#!/usr/bin/env nextflow

params.fasta = '/home/rnahub/rnahub/example/example.fa' // Adjust as needed
params.blast_db = 'pdbnt'
params.rfam_db = '/home/rnahub/rnahub/db/rfam/Rfam.cm'
params.rfam_seed = '/home/rnahub/rnahub/db/rfam/Rfam.seed.gz'

process check_blast {
    input:
    path fasta

    output:
    path "blast_result.txt"

    ///home/rnahub/rnahub/
    """
    python3 ../../../search_blast.py /home/rnahub/rnahub/example/example.fa --db ${params.blast_db} -v > blast_result.txt
    """
}

process check_rfam {
    input:
    path fasta

    output:
    path "rfam_result.txt"

    when:
    // Define the condition to check the output of check_blast
    // For example, checking if the blast_result.txt is empty or contains specific content
    // The 'file' variable here should refer to the output of check_blast

    script:
    """
    python3 /home/rnahub/rnahub/search_rfam.py --rfam ${params.rfam_db} --fasta $fasta --seed ${params.rfam_seed} > rfam_result.txt
    """
}

workflow {
    fasta_file = file(params.fasta)
    blast_result = check_blast(fasta_file)
    if ('no hit') {
        check_rfam(fasta_file)
    }
}
