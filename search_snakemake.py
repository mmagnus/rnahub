"""
File for using snakemake to manage search flows
"""

import os





class RunSnakemake():
    
    def __init__(self) -> None:
        pass
    
    def run_snakemake(self):
        command:str = f'snakemake --cores 2'
        result = os.popen(cmd=command)
        raw_result:str = result.read()
        print(raw_result)
        
        
if __name__ == "__main__":
    
    snakemake:RunSnakemake = RunSnakemake()
    snakemake.run_snakemake()
    