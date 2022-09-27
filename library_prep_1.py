"""
First step of a library preparation
"""

import itertools
from pprint import pprint

def oneStringFile(seq, new_f):
    with open(seq) as seq, open(new_f, 'w') as fasta:
        lines = seq.readlines()
        for id in range(len(lines)-1):
            if lines[id][0] == '>':
                newl = lines[id].replace('\n', '')
                fasta.write(newl)
            elif lines[id][0] != '>' and lines[id+1][0] != '>':
                new = lines[id].replace('\n', '')
                fasta.write(new)
            else:
                fasta.write(lines[id])

run1 = oneStringFile('/home/sola/py_projects/igcat/igcat/data/germline/human/vl.fasta', 'vl_won.fasta')
run1

def splitByReg(seq, reg):
    with open(seq) as seq, open(reg) as reg:
        seq, reg = seq.readlines(), reg.readlines()
        fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4 = [], [], [], [], [], [], []
        for id, item in enumerate(reg):
            seq1 = seq[id]
            lind = seq1.find('IGHJ') + 4
            reg1 = reg[id].split()
            reg1 = list(map(int, reg1[1:]))
            #print(reg1)
            fr1.append(seq1[lind + reg1[0]:lind+ reg1[1]])
            cdr1.append(seq1[lind + reg1[2]:lind + reg1[3]])
            fr2.append(seq1[lind + reg1[4]:lind + reg1[5]])
            cdr2.append(seq1[lind + reg1[6]:lind + reg1[7]])
            fr3.append(seq1[lind + reg1[8]:lind + reg1[9]])
            cdr3.append(seq1[lind + reg1[10]:lind + reg1[11]])
            fr4.append(seq1[lind + reg1[12]:lind + reg1[13]])
    return [[list(set(fr1)), list(set(cdr1)), list(set(fr2)), list(set(cdr2))], [list(set(fr3)), list(set(cdr3))]]

run2 = splitByReg('vl_won.fasta', '/home/sola/py_projects/igcat/igcat/data/nomenclature/human/vl.kabat')
#run2

def genSeq(inputdata, file, file2):
    save = []
    result = list(itertools.product(*inputdata))
    for c in result:
        save.append(c)
    with open(file, 'w') as fasta, open(file2, 'w') as regions:
        for i in save:
            fasta.write('\t'.join(i)+'\t'+'\n')
            regions.write('\t'.join(i))
            #regions.write('\n')

run31 = genSeq(run2[0], 'combinations___.fasta', 'regions___.txt')
run32 = genSeq(run2[1], 'combinations2___.fasta', 'regions2___.txt')
print(run32)
