#!/usr/bin/env python3
import os, sys, argparse, warnings, csv
import numpy as np
import pandas as pd
import subprocess
from collections import Counter

#batch entranzは遺伝子IDの重複を考慮しないので、このスクリプトで補正をします。
#sys.argv[1]: geneid.txt <- get_geneid.shで取得した一覧表
#sys.argv[2]: _edge_gene.fasta <- batch entranzで取得した遺伝子群

geneid = pd.read_table(sys.argv[1], names = ['gene'])
c = Counter(geneid.gene)
_ = [i for i in c.items() if i[1] >= 2]
for line in _:
    for j in range(line[1]-1):
        subprocess.run(f"seqkit grep -nrp {line[0]} {sys.argv[2]} >> dupseq.fasta", shell=True)

subprocess.run(f'cat {sys.argv[2]} dupseq.fasta > ./tmp.fasta', shell=True)
subprocess.run(f'seqkit rename tmp.fasta > edge_gene.fasta', shell=True)
subprocess.run(f'rm dupseq.fasta tmp.fasta', shell=True)
