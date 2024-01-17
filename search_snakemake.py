#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File for using snakemake to manage search flows.

The provided code defines a Python class, RunSnakemake, that is used for running Snakemake workflows. When the run_snakemake() method is invoked, the "/results" directory is deleted if it exists. This is done to clean any old files. If the deletion fails, an exception is raised, otherwise, a success message is printed on the console. Afterward, the Snakemake workflow is executed with a command sent to the system's shell. The command runs Snakemake with two cores. The output of the command is read and then printed to the console. In the main part of the script, an instance of RunSnakemake is created and the run_snakemake() method is called.
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
        
        command:str = f'snakemake -p --cores 2'
        result = os.popen(cmd=command)
        raw_result:str = result.read()
        print(raw_result)
        
        
if __name__ == "__main__":
    
    snakemake:RunSnakemake = RunSnakemake()
    snakemake.run_snakemake()
    
