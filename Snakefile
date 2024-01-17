
configfile: "config.yaml"

rule all:
    input:
        "results/rfam_result.txt",
        "results/blast_result.txt"

rule check_blast:
    input:
        fasta=expand("{seq}", seq=config["fasta"])
    output:
        "results/blast_result.txt"
    shell:
        'python3 ./search_blast.py {input.fasta} --db pdbnt -v > {output}'

rule check_rfam:
    input:
        rfam="/home/rnahub/rnahub/db/rfam/Rfam.cm",
        fasta=expand("{seq}", seq=config["fasta"]),
        seed="/home/rnahub/rnahub/db/rfam/Rfam.seed.gz"
    output:
        "results/rfam_result.txt"
    shell:
        'python3 ./search_rfam.py --rfam {input.rfam} --fasta {input.fasta} --seed {input.seed} > {output}'

