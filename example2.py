# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 2020

@author: E. Agustoni (GitHub: agus-struct)
"""

from motif import motif

def parse_fasta(fasta_filepath): 
    fasta_dict = {}
    name   = ""
    seq    = ""
    with open(fasta_filepath, "r") as fh:
        for line in fh:
            line = line.strip().replace("-", "")
            if line.startswith(">"):
                # Add previous to db
                if name != "": fasta_dict[name.split(">")[-1]] = seq.replace(" ", "").replace("*", "")
                # Start saving current entry
                name = line
                seq = ""
            else:
                seq = "".join((seq, line))
                
    # For last sequence, after iteration finished
    fasta_dict[name.split(">")[-1]] = seq.replace(" ", "").replace("*", "")
    
    return fasta_dict, list(fasta_dict), list(fasta_dict.values())


test_sequences = []
with open(r".\input\PF00072_seed.txt", "r") as f:
    lines = f.readlines()
    for l in lines:
        m = l.split(" ")[-1] # To get rid of the  names in the seed file
        test_sequences.append(m[-18:-13])

test_protein_sequence = parse_fasta(r".\input\Q9ABT2.fasta")[2][0]


test_motif = motif(test_sequences)

scan_output = test_motif.scan(test_protein_sequence,
                              plot=True)

sorted_scan_output = sorted(scan_output,
                            key=lambda x: x[0],
                            reverse=True)
    
print("n\tscore\tpos\tfragment")
for i, j in enumerate(sorted_scan_output[:5]):
    print("%i\t%.2f\t%i\t%s"%(i+1, j[0], j[1], j[2]))
