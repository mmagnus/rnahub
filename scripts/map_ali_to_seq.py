#!/usr/bin/env python3


import sys
import os
import re
import subprocess

def parse_fasta(case, name, ftype, sqfile):
    seq = ""
    with open(sqfile, 'r') as fh: # open the .sto file for reading
        for line in fh: # go through every line of the file
            if ">" in line: 
                seqname = line.split()[0]
            else:
                seq += line.split()[0]

    print("DIR", dirname)
    
    if case == "UI_exon_DI" and ftype == "unflanked":
        m = re.search(r'\.(\d+)-(\d+)/(\d+)-(\d+)/(\d+)-(\d+)$', seqname)
        exonU_start   = 0
        exonU_end     = 0
        intronU_start = int(m.group(1)) 
        intronU_end   = int(m.group(2))
        exon_start    = int(m.group(3)) 
        exon_end      = int(m.group(4))
        intronD_start = int(m.group(5)) 
        intronD_end   = int(m.group(6))
        exonD_start   = 0
        exonD_end     = 0
    elif case == "UE_UI_exon_DI_DE" and ftype == "flanked":
        m = re.search(r'\.(\d+)-(\d+)/(\d+)-(\d+)/(\d+)-(\d+)/(\d+)-(\d+)/(\d+)-(\d+)$', seqname)
        exonU_start   = int(m.group(1)) 
        exonU_end     = int(m.group(2))
        intronU_start = int(m.group(3)) 
        intronU_end   = int(m.group(4))
        exon_start    = int(m.group(5)) 
        exon_end      = int(m.group(6))
        intronD_start = int(m.group(7)) 
        intronD_end   = int(m.group(8))
        exonD_start   = int(m.group(9)) 
        exonD_end     = int(m.group(10))
    elif case == "UE_UI_exon" and ftype == "flanked":
        m = re.search(r'\.(\d+)-(\d+)/(\d+)-(\d+)/(\d+)-(\d+)$', seqname)
        exonU_start   = int(m.group(1)) 
        exonU_end     = int(m.group(2))
        intronU_start = int(m.group(3)) 
        intronU_end   = int(m.group(4))
        exon_start    = int(m.group(5)) 
        exon_end      = int(m.group(6))
        intronD_start = 0
        intronD_end   = 0
        exonD_start   = 0
        exonD_end     = 0
    elif case == "exon_DI_DE" and ftype == "flanked":
        m = re.search(r'\.(\d+)-(\d+)/(\d+)-(\d+)/(\d+)-(\d+)$', seqname)
        exonU_start   = 0 
        exonU_end     = 0
        intronU_start = 0
        intronU_end   = 0
        exon_start    = int(m.group(1)) 
        exon_end      = int(m.group(2))
        intronD_start = int(m.group(3)) 
        intronD_end   = int(m.group(4))
        exonD_start   = int(m.group(5)) 
        exonD_end     = int(m.group(6))
    else:
         print("could not find this dir", dirname)
         exit()
         
    print("exon start/end", exon_start, exon_end)
    print("exonU start/end", exonU_start, exonU_end,"exonD start/end", exonD_start, exonD_end)
    print("intronU start/end", intronU_start, intronU_end,"intronD start/end", intronD_start, intronD_end)
            
    return seq, seqname, exonU_start, exonU_end, intronU_start, intronU_end, exon_start, exon_end, intronD_start, intronD_end, exonD_start, exonD_end


def parse_msa(msafile):
    '''
    Parse the stockholm file made after 1 round of nhmmer
    using the flanked sequence

    Parameters
    ----------
    infile : string
        path to stockholm file made from 1 round of nhmmer
        on the flanked sequence

    Returns
    -------
    seq : string
        reference sequence for fugu in the .sto alignment
    start : integer
        coordinate of start of seq
    end : integer
        coordinate of end of seq

    '''
    query = "" # emtpy strings to append to
    identifier = "" # emtpy strings to append to
    aseq = ""
    with open(msafile, 'r') as fh: # open the .sto file for reading
        for line in fh: # go through every line of the file
            if "ID" in line: # take the reference sequence (first instance of fugu)
                this = line.split() # split on white space
                query += this[2] # take the identifier with coordinates
                break # break to avoid taking second, third, etc. instances of fugu
        
        for line in fh: # go through every line of the file
            if "Homo sapiens" in line: # take the reference sequence (first instance of fugu)
                this = line.split() # split on white space
                identifier += this[1] # take the identifier with coordinates
                break # break to avoid taking second, third, etc. instances of fugu

        for line in fh: # go through the file
            if line.startswith(identifier): # and search for the identifier in the alignment portion
                aseq += line.split()[1] # get the sequence from the alignment
        print("\n^^aseq", aseq);
    
    return aseq   


dirname = sys.argv[1]
print("dirname", dirname)

# sorek.UI_exon_DI/D12166_7.UI_exon_DI_unflanked
# sorek.UE_UI_exon_DI_DE/D12166_7.UE_UI_exon_DI_DE_flanked
m = re.search(r'sorek\.(\S+)\/(\S+\.\S+)_(flanked|unflanked)', dirname)
case  = m.group(1)
name  = m.group(2)
ftype = m.group(3)
print("case", case, "name", name, "ftype", ftype)

m = re.search(r'(\S+)\_(unflanked|flanked)', dirname)
sqfile = m.group(1)+".fa"
print("^^sqfile", sqfile)

sq, sqname, UE_from, UE_to, UI_from, UI_to, exon_from, exon_to, DI_from, DI_to, DE_from, DE_to = parse_fasta(case, name, ftype, sqfile)
print("\n^^sqname", sqname)
print("EU start/end", UE_from, UE_to)
print("intronU start/end", UI_from, UI_to)
print("exon start/end", exon_from, exon_to)
print("intronD start/end", DI_from, DI_to)
print("exonD start/end", DE_from, DE_to)
sq_L = len(sq)
print("\n^^sq", sq_L, sq)

msafile = dirname+"/rm_v3.sto"
print("msafile", msafile)
aseq = parse_msa(msafile)
aseq_alen = len(aseq)
print("\n^^aseq", aseq_alen, aseq)
aseq_seq = ""
aseq_seq_L = 0
for i in range(aseq_alen): 
    if aseq[i].isalpha():
        aseq_seq_L += 1
        aseq_seq+=aseq[i]
print("\n^^aseq_seq", aseq_seq_L, aseq_seq)

seqfile1 = "foo1"
with open(seqfile1, 'w') as f: # open the .sto file for reading
    print('>{}\n{}\n'.format("seq",sq),file=f)

seqfile2 = "foo2"
with open(seqfile2, 'w') as f: # open the .sto file for reading
    print('>{}\n{}\n'.format(name,aseq_seq),file=f)

outfile = "foo3"
os.system("/n/eddy_lab/users/erivas/rscape_v2.0.4.b/lib/hmmer/src/nhmmer foo1 foo2 > foo3")
os.remove(seqfile1)
os.remove(seqfile2)

start = 0
with open(outfile, 'r') as f: 
    for line in f: 
        if   "Query:" in line:
            continue
        elif " seq " in line:
            start = int(line.split()[1])
            break
end = 0
with open(outfile, 'r') as f: 
    for line in f: 
        if   "Query:" in line:
            continue
        elif " seq " in line:
            end = int(line.split()[3])
os.remove(outfile)
print("start", start)
print("end", end)

map = []
mapback = []
coord = start-1
for j in range(coord):
    mapback.append(0)    
for i in range(aseq_alen): 
    if aseq[i].isalpha(): 
        map.append(coord)
        mapback.append(i+1)
        coord += 1 
    else:
        map.append(0)
for j in range (sq_L-coord):
    mapback.append(0)
    
#print(map)
print(len(mapback),mapback)

if case == "UI_exon_DI" and ftype == "unflanked":
    print("UI",   UI_from,   UI_to,   "(", UI_from-UI_from,   UI_to-UI_from,   ")", mapback[UI_from-UI_from],   mapback[UI_to-UI_from])
    print("exon", exon_from, exon_to, "(", exon_from-UI_from, exon_to-UI_from, ")", mapback[exon_from-UI_from], mapback[exon_to-UI_from])
    print("DI",   DI_from,   DI_to,   "(", DI_from-UI_from,   DI_to-UI_from,   ")", mapback[DI_from-UI_from],   mapback[DI_to-UI_from])

if case == "UE_UI_exon_DI_DE" and ftype == "flanked":
    print("UI",   UI_from,   UI_to,   "(", UI_from-UI_from,   UI_to-UI_from,   ")", mapback[UI_from-UI_from],   mapback[UI_to-UI_from])
    print("exon", exon_from, exon_to, "(", exon_from-UI_from, exon_to-UI_from, ")", mapback[exon_from-UI_from], mapback[exon_to-UI_from])
    print("DI",   DI_from,   DI_to,   "(", DI_from-UI_from,   DI_to-UI_from,   ")", mapback[DI_from-UI_from],   mapback[DI_to-UI_from])

if case == "UE_UI_exon" and ftype == "flanked":
      print("UI",   UI_from,   UI_to,   "(", UI_from-UI_from,   UI_to-UI_from,   ")", mapback[UI_from-UI_from],   mapback[UI_to-UI_from])
  
if case == "exon_DI_DE" and ftype == "flanked":
    print("DI",   DI_from,   DI_to,   "(", DI_from-UI_from,   DI_to-UI_from,   ")", mapback[DI_from-UI_from],   mapback[DI_to-UI_from])

