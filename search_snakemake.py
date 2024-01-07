"""
File for using snakemake to manage search flows
"""

import os
from pathlib import Path
import shutil





class RunSnakemake():
    
    def __init__(self) -> None:
        pass
    
    def run_snakemake(self):
        
        #first clean up old files
        results_dir:Path = Path('./results')
        is_present:bool = os.path.isdir(results_dir)
        if is_present is True:
            try:
                shutil.rmtree(results_dir, ignore_errors=False)
                print("Successfully cleaned up old files")
            except:
                raise Exception(f'Could not delete folder {results_dir}')
        
        
        command:str = f'snakemake --cores 2'
        result = os.popen(cmd=command)
        raw_result:str = result.read()
        print(raw_result)
        
        
if __name__ == "__main__":
    
    snakemake:RunSnakemake = RunSnakemake()
    snakemake.run_snakemake()
    