#!/usr/bin/env python
import sys
import os

def main():
    if len(sys.argv) != 3:
        print("Illegal number of parameters")
        print("Should be: python script.py <alignment> <genome specific string>")
        sys.exit(1)

    # names of files, provided based on first two arguments
    alignment = sys.argv[1]
    seq_spec_string = sys.argv[2]

    # directory of alignment file that you provide, the output will go into there
    DIR = os.path.dirname(alignment)

    # Read the alignment file to get the top sequence that contains reference genome
    seq_name = None
    with open(alignment, 'r') as f:
        for line in f:
            if seq_spec_string in line:
                seq_name = line.split()[1]
                break

    if not seq_name:
        print(f"No sequence found matching {seq_spec_string}")
        sys.exit(1)

    # Extract only the sequence (gaps included) of that sequence from the MSA
    sequence = []
    with open(alignment, 'r') as f:
        for line in f:
            if seq_name in line and not '/' in line.split()[1]:
                sequence.append(line.split()[1])

    bar_txt_path = os.path.join(DIR, 'bar.txt')
    with open(bar_txt_path, 'w') as f:
        f.write(''.join(sequence))

    # Now create a counter that loops through and increments as it reads each character.
    bp = []
    array = []
    counter = 0

    with open(bar_txt_path, 'r') as f:
        while True:
            c = f.read(1)
            if not c:
                break
            counter += 1
            if c in ['U', 'A', 'T', 'G', 'C']:
                if c == 'U':
                    c = 'T'
                bp.append(c)
                array.append(counter)

    # Write header and results to bp_col.txt
    bp_col_txt_path = os.path.join(DIR, 'bp_col.txt')
    with open(bp_col_txt_path, 'w') as f:
        f.write("bp \t Col\n")
        for i, each in enumerate(array, start=1):
            f.write(f"{i} \t {each}\n")

    # Remove the temporary bar.txt file
    os.remove(bar_txt_path)

if __name__ == "__main__":
    main()
