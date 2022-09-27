"""
A file combining all FR1 - CDR3 regions and creates unique IDs
The output is 1000 random sequences
"""
import uuid
import random

def unique_ids(seqs):
    ids = []
    for i in seqs:
        ids.append(uuid.uuid4())
    #print(ids)
    return ids

def combine_all_seq(beg_file, end_file):
    f1000, regs = [], []
    with open(beg_file) as f1, open(end_file) as f2:
        f1, f2 = f1.readlines(), f2.readlines()
        ff1, ff2 = random.sample(f1, 1000), random.sample(f2, 1000)
        for i in range(1000):
            ap = ''.join([ff1[i][:-1], ff2[i][:-1]]).split()
            #print(ap)
            f1000.append(''.join([ff1[i].replace('\t', '')[:-1], ff2[i].replace('\t', '')[:-1]]))
            regs1 = []

            for i1 in range(6):
                if i1 == 0:
                    a = len(ap[i1])+1
                    regs1.append('\t'.join(['1', str(a)]))
                elif i1 > 0:
                    ai = a
                    a1 = len(ap[i1-1])
                    a2 = len(ap[i1])
                    a = ai + a2+1
                    regs1.append('\t'.join([str(ai+1), str(a)]))
            regs.append(regs1)
    return f1000, regs

def wr_res(file1, file2, beg_file, end_file):
    gen = combine_all_seq(beg_file, end_file)
    seq = gen[0]
    reg = gen[1]
    ids = unique_ids(seq)
    with open(file1, 'w') as file1, open(file2, 'w') as file2:
        for i in range(1000):
            seq1 = [seq[i][y - 80:y] for y in range(80, len(seq[i]) + 80, 80)]
            line = ''.join(['>LIGHT_', str(ids[i]), '\n'])
            print(reg[i])
            file1.write(line)
            for s in seq1:
                file1.write(''.join([s, '\n']))
            file2.write(line)
            file2.write('\t'.join(reg[i]))
            file2.write('\n')

wr_res('LIGHT_id_comb1.fasta', 'LIGHT_id_reg.fasta', 'combinations___.fasta', 'combinations2___.fasta')