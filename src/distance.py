import Bio
print(Bio.__version__)



from Bio.Seq import Seq
my_seq = Seq("AGTCACAGTTTT")
print(my_seq)

from Bio import SeqIO

1
    
from Bio import SeqIO
v = SeqIO.parse("../data/set1c.fasta", "fasta")
for seq_record in v:
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
;

def hamming_distance(seq1, seq2):
    return sum(((n1 != n2) and (n1 != '-') and (n2 != '-')) for n1, n2 in zip(seq1, seq2))
help(zip)

list(zip(x,x))

;

import matplotlib.pyplot as plt

s = list(SeqIO.parse("../data/set1c.fasta", "fasta")); s
s[0]
s[4]

for i in range(0,143):
    a = 1000*i
    b = 1000*(i+1)-1
    print("distance for partition ",i, ": ",hamming_distance(s[1][a:b],s[2][a:b]))

import math
def f(x):
    return (3/4) * (1 - math.exp(-(4/3) * x))

seq1 = s[4]
seq2 = s[3]
for i in range(0,143):
    a = 1000*i
    b = 1000*(i+1)-1
hams = [hamming_distance(seq1[1000*i:1000*(i+1)-1],seq2[1000*i:1000*(i+1)-1]) for i in range(0,143)]
exp_distances = [f(x) for x in hams]

plt.title("hamming distances for partitions of "+ seq1.description + " and " + seq2.description)
plt.hist(exp_distances)


plt.show()
u = list(zip(seq1[a:b],seq2[a:b]))
u[1][0]

;





