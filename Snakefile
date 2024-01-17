# run first blast on pdb and then rfam on the results
configfile: "config.yaml"

rule all:
    input:
        "results/blast_result.txt",
        "results/rfam_result.txt"

checkpoint run_blast:
    input:
        fasta=expand("{seq}", seq=config["fasta"])
    output:
        "results/blast_result.txt"
    shell:
        'python3 ./search_blast.py {input.fasta} --db pdbnt -v > {output}'

def should_run_rfam(wildcards):
    with open(checkpoints.run_blast.get(**wildcards).output[0], 'r') as f:
        # Define the condition for a no-hit scenario. For instance, checking if the file is empty or contains a specific string.
        if 'No hits found' in f.read():
            return expand("{seq}", seq=config["fasta"])
        else:
            return "dummy_file.txt"  # This should be a dummy file that exists.

rule run_rfam:
    input:
        rfam="/home/rnahub/rnahub/db/rfam/Rfam.cm",
        fasta=should_run_rfam,
        seed="/home/rnahub/rnahub/db/rfam/Rfam.seed.gz"
    output:
        "results/rfam_result.txt"
    shell:
        "python3 ./search_rfam.py --rfam {input.rfam} --fasta {input.fasta} --seed {input.seed} > {output};"
        "head {output}"



