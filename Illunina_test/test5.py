
from pathlib import Path
import re
import time
import sys

time1 = time.time()
FASTQ = Path("/Users/matsuoryuushi/Desktop/SE11/Illumina_SE11.fastq")
with FASTQ.open() as f:
    fastq = f.readlines()

head = []
seq = []
qual = []
for i in range(len(fastq)):
    if i % 4 == 0:
        head.append(fastq[i][1:].rstrip("\n"))
    elif i % 4 == 1:
        seq.append(fastq[i].rstrip("\n"))
    elif i % 4 == 3:
        qual.append(fastq[i].rstrip("\n"))

#参照配列
REF = Path("/Users/matsuoryuushi/Desktop/SE11/ref_SE11.fasta")
with REF.open() as f:
    ref = f.readlines()

ref_seq = {}
name = []
count_before = 0
for i in range(len(ref)):
    if ref[i][0] == ">":
        if i != 0:
            ref_seq[name.group()] = "".join(ref[count_before + 1:i])
            count_before = i
        name = re.search(r".+?(?=\s)", ref[i])
    else:
        ref[i] = ref[i].rstrip("\n")
    
    if i == len(ref) - 1:
        ref_seq[name.group()] = "".join(ref[count_before + 1:i])

chr_name = {}
for i, k in enumerate(ref_seq):
    chr_name[i + 1] = k[1:]


#algorithm 1

compliment = {"A":"T", "C":"G", "G":"C", "T":"A"}
def minimizerskech(s, w, k):
    M = []
    n = len(s)
    s_compliment_list = []
    for i in range(n):
        s_compliment_list.append(compliment[s[i]])
    s_compliment = "".join(s_compliment_list)

    for i in range(1, n - w - k + 2):
        u_kmer = s[i:i + k]
        v_kmer = s_compliment[i:i + k][::-1]
        u = hash(u_kmer)
        v = hash(v_kmer)
        if u < v:
            M += [[u, i, 0]]
        elif v < u:
            M += [[v, i, 1]]


    return M


# algorithm 3
def indextarget(Tseq, w, k):
    H = {}
    for t, t_key in enumerate(Tseq):
        M = minimizerskech(Tseq[t_key], w, k)
        for h, i, r in M:
            H.setdefault(h, []).append([t + 1, i, r])
    
    return H


# algorithm 4
def mapping(H, q, w, k):
    A = []
    M = minimizerskech(q, w, k)
    for h, i, r in M:
        if h in H:
            for t, ih, rh in H[h]:
                if r == rh:
                    A += [[t, 0, i - ih, ih, i]]
                else:
                    A += [[t, 1, i + ih, ih, i]]

    
    A.sort(key=lambda x: (x[0], x[1], x[3]))

    C = []
    C_current = [A[0]]
    for e in A:
        e_last = C_current[-1]
        if e[0] == e_last[0] and e[1] == e_last[1] and abs(e[2] - e_last[2]) <= 5:
            C_current.append(e)
        else:
            C.append(C_current)
            C_current = [e]

    C.append(C_current)

    C_max_0 = max(C, key=len)
    C_max = C_max_0[0]

    # -ストランドの変異あり
    if C_max[1] == 1 and C_max[4] < (len(q) - kmer):
        C_max[3] = C_max[3] - ((len(q) - kmer) - C_max[4])
    return C_max




window = 1
kmer = 20


hash_ref = indextarget(ref_seq, window, kmer)

time6 = time.time()



for i in range(len(seq)):
    ans = mapping(hash_ref, seq[i], window, kmer)
    if ans[1] == 0:
        print(head[i] + "\t" + chr_name[ans[0]] + "\t" + str(abs(ans[2]) + 1) + "\t" + "+")
    else:
        print(head[i] + "\t" + chr_name[ans[0]] + "\t" + str(ans[3] + 1) + "\t" + "-")

time7 = time.time()



print(time7 - time1, time7 - time6, file=sys.stderr)

