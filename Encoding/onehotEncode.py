from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import sys

tbs = pd.Series(list('AGCT'))
hot_table= pd.get_dummies(tbs)
hot_table['X']=[0,0,0,0]
hash_code= { x:y for x, y in zip(hot_table.columns.values,  hot_table.transpose().values) }

def ecd_seq2(t):
    t3=t.apply(lambda x: [i for i in x])
    mat_pos=np.empty((0,4), int)
    for i in t3:
        m = [hash_code[c] for c in i]
        m1=np.array(m)
        mat_pos=np.vstack((mat_pos,m1))
    return(mat_pos)

rt=sys.argv[1]
fl1=sys.argv[2]
fl2=sys.argv[3]
nm=sys.argv[4]

records = list(SeqIO.parse(rt+fl1, "fasta"))
maxn=9000  ###9063-63
mat_seq=list()
stp=200
wnsz=1000
for rcd in records:
    sq0=str(rcd.seq)
    sqa=sq0.upper()
    sq1=sqa.replace("N","")
    if len(sq1)<maxn:
        sq=sq1+(-len(sq1)+maxn)*'X'
    else:
        sq=sq1[0:maxn]
    t=[sq[i:i+wnsz] for i in range(0, len(sq)-wnsz,stp)]
    mat_seq.append(t)
    

m_all=list()

for k in mat_seq:
    mt=list()
    for xx in k:
        t=pd.DataFrame([xx])[0]
        m=ecd_seq2(t)  #.transpose()
        mt.append(m.tolist())
    m_all.append(mt)

mat_all_pos=np.array(m_all)

records = list(SeqIO.parse(rt+fl2, "fasta"))

mat_seq=list()

for rcd in records:
    sq0=str(rcd.seq)
    sqa=sq0.upper()
    sq1=sqa.replace("N","")
    if len(sq1)<maxn:
        sq=sq1+(-len(sq1)+maxn)*'X'
    else:
        sq=sq1[0:maxn]
    t=[sq[i:i+wnsz] for i in range(0, len(sq)-wnsz,stp)]
    mat_seq.append(t)

m_all=list()

for k in mat_seq:
    mt=list()
    for xx in k:
        t=pd.DataFrame([xx])[0]
        m=ecd_seq2(t) #.transpose()
        mt.append(m.tolist())
    m_all.append(mt)

mat_all_neg=np.array(m_all)

np.save(rt+nm+'_pos_1hotFW.npy',mat_all_pos)
np.save(rt+nm+'_neg_1hotFW.npy',mat_all_neg)



