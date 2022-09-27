def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

with open("/home/sola/affirmed_germlines_Human_gapped.fa") as fa:
    lines = fa.readlines()
    for line in range(len(lines)):
        if lines[line][0] == '>':
            to_j = lines[line+1:line+6]
            to_j2 = []
            for j in to_j:
                j2 = j.replace('\n', '')
                to_j2.append(j2)
            new_line = ''.join(to_j2)

            # extract regions
            fr1_beg = 1
            frags = list(find_all(new_line, '.'))
            needed_frags = []
            for f in range(len(frags)-1):
                if frags[f+1]-frags[f] != 1:
                    needed_frags +=(frags[f], frags[f+1])
            fr1_end = int(frags[0])+1
            cdr1_beg = int(needed_frags[0]+1)
            cdr1_end = int(needed_frags[1]+1)
            print(fr1_beg, fr1_end, cdr1_beg)

            print(new_line[cdr1_beg:cdr1_end])



