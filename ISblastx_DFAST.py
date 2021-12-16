#!/usr/bin/env python3
import sys, os, argparse
import subprocess
from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_args():
    parser = argparse.ArgumentParser(description='dkato. November, 2021') 
    parser.add_argument('-g' , dest ='genome', required=True,
                        help='path to genome.fasta from DFAST. The contig ID is unified in the form of "sequence~"') 
    parser.add_argument('-f' , dest ='features', required=True,
                        help='path to features.tsv from DFAST') 
    
    parser.add_argument('-e', '--evalue', type=float, default=.001, 
                       help='evalue in blastx.(default:0.001)')
    parser.add_argument('-th', '--threshold', type=int, default=0, 
                       help='minimum length of interval sequence as input of blastx.(default:0)')
    
    parser.add_argument('-db', '--database', default='/home_ssd/local/db/blastdb.20200904/nr', 
                       help='path to your nr database.(default:/home_ssd/local/db/blastdb.20200904/nr)')
    parser.add_argument('-t' , '--num_threads', type=int, default=48,
                       help='num threads in blastx.(default:16)',) 
    parser.add_argument('-nd', '--num_descriptions', type=int, default=50,
                       help='num descriptions in blastx.(default:50)')
    
    parser.add_argument('-x',dest='x',  type=bool, default=False, choices=[True, False],
                        help='If "True", blastx will not be executed.(default:False)') 
    return parser.parse_args()

def get_edited_features(path_to_features = None):
    _df = pd.read_table(path_to_features)
    start , end = [], []
    
    for i, df_i in enumerate(_df.location):
        temp = df_i.split
        if df_i[0] != 'c':
            start += [temp('..')[0]]
            end += [temp('..')[1]]
        else:
            start += [temp('..')[1][:-1]]
            end += [temp('..')[0][len('complement('):]]
    cds = pd.DataFrame(columns = ['start', 'end'])
    cds.start, cds.end = start, end
    pd.concat([cds, _df], axis = 1).to_csv("new_feature.csv")
    return  pd.concat([_df.sequence, cds], axis = 1)

class MyGetIS(object):
    def __init__(self, df = None, genome = None):
        self.df        = df 
        self.genome    = genome
        self.counter   = 0
        self.id_out    = []
        self.seq_out   = []
        self.fasta_out = []
        self.dup_len   = 0
        self.double_interval = 0
    
    @staticmethod
    def judge_inside(i, last_iter):
        return i != 0 and i != last_iter
    
    @staticmethod
    def judge_normal_interval(cds_h, cds_i):
        return max(int(cds_h.start), int(cds_h.end)) < min(int(cds_i.start), int(cds_i.end))
    
    @classmethod
    def get_1st_interval(cls, cds_i, 
                         genome_i, 
                         id_out  = None, 
                         seq_out = None, 
                         counter = None):
        start = 1
        end = min(int(cds_i.start), int(cds_i.end)) -1
        tmp = genome_i.seq[:end]
        seq_out += [tmp] ;counter+=1
        id_out += [f'interval_region{counter}_LEN{len(tmp)}_{genome_i.id}_{start}_{end}']
        return counter

    @classmethod
    def get_interval_inside(cls, cds_h, cds_i, 
                            genome_i, 
                            id_out = None, 
                            seq_out = None, 
                            counter = None):
        start = max(int(cds_h.start), int(cds_h.end))+1
        end = min(int(cds_i.start), int(cds_i.end))-1
        tmp = genome_i.seq[start-1:end]
        seq_out += [tmp] ;counter+=1
        id_out += [f'interval_region{counter}_LEN{len(tmp)}_{genome_i.id}_{start}_{end}']
        return counter

    @classmethod
    def get_last_interval(cls, cds_h, 
                          genome_h, 
                          id_out = None,
                          seq_out = None, 
                          counter = None):
        start = max(int(cds_h.start), int(cds_h.end))+1
        end = len(genome_h.seq)
        tmp = genome_h.seq[start-1:]
        seq_out += [tmp] ;counter+=1
        id_out += [f'interval_region{counter}_LEN{len(tmp)}_{genome_h.id}_{start}_{end}']
        return counter

    def get_IS(self):
        for i in tqdm(range(len(self.df))):#
            if i == 0:
                cds_i = self.df.loc[i, ] ;tmpi = cds_i.sequence
                
                genome_i = [genome_i for genome_i in self.genome if genome_i.id==tmpi][0]
                self.counter = MyGetIS.get_1st_interval(cds_i, genome_i, 
                                                        id_out = self.id_out, 
                                                        seq_out = self.seq_out, 
                                                        counter = self.counter)

            elif MyGetIS.judge_inside(i, len(self.df)-1):
                cds_h = self.df.loc[i-1, ] ;tmph = cds_h.sequence
                cds_i = self.df.loc[i, ]   ;tmpi = cds_i.sequence

                genome_h = [genome_h for genome_h in self.genome if genome_h.id==tmph][0]
                genome_i = [genome_i for genome_i in self.genome if genome_i.id==tmpi][0]
                
                if tmph == tmpi:
                    if MyGetIS.judge_normal_interval(cds_h, cds_i):
                        self.counter = MyGetIS.get_interval_inside(cds_h, cds_i, genome_i, 
                                                                   id_out = self.id_out, 
                                                                   seq_out = self.seq_out, 
                                                                   counter = self.counter)
                    else: # duplicate CDS
                        def dif(x, y):
                            return max(int(x.start), int(x.end)) - min(int(y.start), int(y.end)) + 1
                        self.dup_len += dif(cds_h, cds_i)
                        
                elif tmph != tmpi:
                    def be_last_interval(cds_h, genome_h):
                        return max(int(cds_h.start), int(cds_h.end)) != len(genome_h.seq)
                    if be_last_interval(cds_h, genome_h):
                        self.counter = MyGetIS.get_last_interval(cds_h, genome_h, 
                                                                 id_out = self.id_out, 
                                                                 seq_out = self.seq_out, 
                                                                 counter = self.counter)
                    def be_initial_interval(cds_i, genome_i):
                        return min(int(cds_i.start), int(cds_i.end)) != 1
                    if be_initial_interval(cds_i, genome_i):
                        self.counter = MyGetIS.get_1st_interval(cds_i, genome_i, 
                                                                id_out = self.id_out, 
                                                                seq_out = self.seq_out, 
                                                                counter = self.counter)
                else:
                    print('ERROR!') ;sys.exit()
                    
            elif i == len(self.df)-1:
                cds_i = self.df.loc[i, ] ;tmpi = cds_i.sequence
                genome_i = [genome_i for genome_i in self.genome if genome_i.id==tmpi][0]
                self.counter = MyGetIS.get_last_interval(cds_i, genome_i, 
                                                         id_out = self.id_out, 
                                                         seq_out = self.seq_out, 
                                                         counter = self.counter)
            
        #print(self.double_interval)
        return self.id_out, self.seq_out
    
    def main(self, id_out, seq_out):
        for j in range(len(seq_out)):
            self.fasta_out += [SeqRecord(seq_out[j], 
                                         id=id_out[j], 
                                         description=id_out[j])]
            
        ids = [_.id for _ in self.genome]
        ids_isCDS = list(self.df.sequence.drop_duplicates())
        ids_noCDS = list(set(ids) - set(ids_isCDS))
        self.fasta_out = self.fasta_out + [_ for _ in self.genome if _.id in ids_noCDS]
        
        SeqIO.write(self.fasta_out, os.path.join(f'interval_regions.fasta'), "fasta")
        print( "======> Total length of duplicated CDS          :",
              self.dup_len, "bp")
        
def split_fasta(multifasta = None):
    if 'each_IS' not in os.listdir(path='./'):
        os.system('mkdir each_IS')
    
    for line in list(SeqIO.parse(multifasta, "fasta")):
        name_i = line.description
        if " " in name_i or "/" in name_i:
            name_i = name_i.replace(" ", "_").replace("/", "_")
            
        SeqIO.write(line, f"./each_IS/{name_i}.fasta", "fasta")

def blastx(dir_in = None, 
           num_threads = None, 
           num_descriptions = None, 
           db = None,
           threshold = None,
           evalue = None):
    error1 = "specify the path to your nr database with db option. (default:/home_ssd/local/db/blastdb.20200904/nr)"
    assert os.path.exists('/home_ssd/local/db/blastdb.20200904/'), error1
    
    if 'res_ISblastx' not in os.listdir(path='./'):
        os.system('mkdir res_ISblastx')

    for fasta_path in tqdm(os.listdir(path=dir_in)):
        fasta = list(SeqIO.parse(f"./each_IS/{fasta_path}", "fasta"))[0]
        if len(fasta.seq) >= threshold:
            subprocess.run(f"blastx -query ./each_IS/{fasta_path} -out ./res_ISblastx/out_{fasta_path[:-6]}.txt -num_threads {num_threads} -evalue {evalue} -num_descriptions {num_descriptions} -num_alignments 100 -db {db}"
                          ,shell=True)
            #subprocess.run(f"diamond blastx --query ./each_IS/{fasta_path} --out ./res_ISblastx/out_{fasta_path[:-6]}.txt -threads {num_threads} --evalue {evalue} -db {db}"
             #             ,shell=True)
    
def main():
    #load data
    df     = get_edited_features(path_to_features = get_args().features)
    genome = list(SeqIO.parse(get_args().genome, "fasta"))
    
    error2 = 'The contig ID does not match the ID specified in feature.tsv. Please use a genome.fasta, the ID of contig is unified in the form of "sequence~" which is output by DFAST'
    assert df.sequence[0]==genome[0].id, error2 
    
    #get IS
    print('1.extracting IS..')
    instance = MyGetIS(df = df,
                       genome = genome)
    seq_out, id_out = instance.get_IS()
    instance.main(seq_out, id_out)
    split_fasta(multifasta ='./interval_regions.fasta')
    
    #blastx
    if not get_args().x:
        print('2.blastx now..')
        blastx(dir_in = './each_IS/',
               num_threads = get_args().num_threads,
               num_descriptions= get_args().num_descriptions,
               db = get_args().database, 
               threshold = get_args().threshold,
               evalue = get_args().evalue )

if __name__ == "__main__":
    main()
