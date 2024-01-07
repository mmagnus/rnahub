"""
File to hold classes for checking if a specific search has results
before running next search in snakemake
"""




class SnakemakeSearchRules():
    
    def __init__(self) -> None:
       pass
    
    def do_blast_run(self):
        return True
    
    def do_rfam_run(self):
        return True