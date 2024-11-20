#!/bin/bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/rnahub/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/rnahub/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/rnahub/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/rnahub/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
python search_rfam.py example/example.fa example/gly1_query.fa example/xrrna.fa
#python search_blast.py --db pdbnt example/gly1_query.fa example/u6.fa
exit
#./rnahub.py --create-job-folder --cpus 16 --repeatmasker example/xrrna.fa --db /home/rnahub/db/Riboviria_genome_db/20240826_all_Riboviria_db.fa --dev-skip-search
./rnahub.py --create-job-folder --repeatmasker --cpus 16 example/gly1_query.fa --flanked example/gly1_flanked.fa --db /home/rnahub/db/1409_Acomycota_genomes-may19.fa
exit
cd /home/rnahub/rnahub-web
python daemon.py --debug xrRNA-06be51de
exit
./rnahub.py --create-job-folder --repeatmasker --cpus 16 example/gly1_query.fa --flanked example/gly1_flanked.fa --db /home/rnahub/db/1409_Acomycota_genomes-may19.fa
exit
#./rnahub.py  --cpus 16 --fasta example/xrrna.fa --flank --db /home/rnahub/db/Riboviria_genome_db/20240826_all_Riboviria_db.fa
#./rnahub.py  --cpus 16 --fasta example/xrrna_casp.fa  --db /home/rnahub/db/Riboviria_genome_db/20240826_all_Riboviria_db.fa
#./rnahub.py  --cpus 16 --fasta example/xrna_casp_flank.fa --db /home/rnahub/db/Riboviria_genome_db/20240826_all_Riboviria_db.fa
#./rnahub.py  --cpus 16 --flank --fasta example/grc1_intron_exon.fa  --db /home/rnahub/db/1409_Acomycota_genomes-may19.fa
#./rnahub.py  --cpus 16 --fasta example/grc1_intron2.fa  --db /home/rnahub/db/1409_Acomycota_genomes-may19.fa

./rnahub.py --flank --cpus 16 --db /home/rnahub/db/1409_Acomycota_genomes-may19.fa --fasta /home/rnahub/rnahub/example/gly1_flanked.fa --dev-skip-nhmmer123 --dev-skip-nhmmer0 --dry
exit
#
rsync -zarv /Users/magnus/work/src/rnahub/ ody:/n/home06/mmagnus/m/rnahub/  --exclude 'config_local.py' --exclude '.git/' --include="*" #--dry-run #
ssh ody '~/m/rnahub/test_ody.sh' # && python unflank.py db/1409_Acomycota_genomes-may19.fa example/grc1_intron.fa'
#rsync -rv odx:/n/home06/mmagnus/m/rnahub/jobs ~/d/
exit


#!/bin/bash
#./rnahub.py --flank db/Saccharomyces_genomes.fa example/YAR014C_plus_IGR.fa
#./rnahub.py --flank db/Saccharomyces_genomes.fa example/YAR015W_plus_IGR.fa


#rsync -zarv /Users/magnus/work/src/rnahub/ odx:/n/home06/mmagnus/m/rnahub/ --exclude '.git/' --include="*/" --include="*.py" --exclude="*" --exclude=='__pycache__' #--dry-run --delete
#rsync -zarv /Users/magnus/work/src/rnahub/ odx:/n/home06/mmagnus/m/rnahub/ --exclude '.git/' --include="*/" --include="*.sh" --exclude="*" --exclude=='__pycache__' --dry-run # --delete
#
# ./nextflow rnahub.nf

echo 'fasta: example/random.fa' > config.yaml #example.fa
#echo 'fasta: example/example.fa' > config.yaml
#python search_hmmer.py --db /Users/magnus/work/rfam/Rfam.cm example/random.fa
/Users/magnus/work/opt/infernal/infernal-1.1.2-macosx-intel/binaries/cmscan  --tblout cmscan.tblout -o cmscan.out /Users/magnus/work/rfam/Rfam.cm example/random.fa
cat cmscan.tblout
#  --noali
#/Users/magnus/work/rfam/Rfam.cm

function test1 {
   eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
   conda activate RF2NA
   python ./search_blast.py --db pdbnt example/random.fa
   }

function test2 {       
   eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
   conda activate RF2NA
   snakemake -p --cores 2 -F ##
   #./search_snakemake.py
}

#test2

# eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
# conda activate RF2NA
# exit

# #!/bin/bash
# eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
# conda activate RF2NA
# #source /home/rnahub/miniconda3/bin/activate RF2NA
# set -x

# python ./search_blast.py --db pdbnt example/example.fa
# python ./search_blast.py --db refseq example/example.fa
# #python ./search_blast.py --db nt example/example.fa
# python search_rfam.py --rfam /home/rnahub/rnahub/db/rfam/Rfam.cm --fasta example/example.fa --seed /home/rnahub/rnahub/db/rfam/Rfam.seed.gz
