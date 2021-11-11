import os, argparse
from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_args():
    parser = argparse.ArgumentParser(description='My blastx') 
    parser.add_argument('-g' , '--genome', required=True,
                        help='path to genome.fasta from DFAST') 
    parser.add_argument('-f' , '--features', required=True,
                        help='path to features.tsv from DFAST') 
    
    parser.add_argument('-m', '--mode', type=str, choices=['strict', 'loose'], required=True, 
                       help='specify the mode of this program',)
    parser.add_argument('-e', '--evalue', type=float, default=0.0001, 
                       help='evalue')
    parser.add_argument('-th', '--threshold', type=int, default=300, 
                       help='threshold')
    
    parser.add_argument('-db', '--database', default='/home_ssd/local/db/blastdb.20200904/nr', 
                       help='nr database')
    parser.add_argument('-t' , '--num_threads', type=int, default=3,
                       help='num_threads',) 
    parser.add_argument('-nd', '--num_descriptions', type=int, default=50,
                       help='num_descriptions')
    
    parser.add_argument('--Without_blast', type=bool, default=False, choices=[True, False],
                        help='True or False')
    #parser.add_argument('--KYOTO') 
    return parser.parse_args()

def get_edited_features(path_to_features = None):
    _df = pd.read_table(path_to_features)
    start, end = [], []
    
    for i, df_i in enumerate(_df.location):
        if df_i[0] != 'c':
            start += [df_i.split('..')[0]]
            end += [df_i.split('..')[1]]
        else:
            start += [df_i.split('..')[1][:-1]]
            end += [df_i.split('..')[0][len('complement('):]]
    cds = pd.DataFrame(columns = ['start', 'end'])
    cds.start, cds.end = start, end
    return  pd.concat([_df.sequence, cds], axis = 1)

class MyGetIS(object):
    def __init__(self, df=None, genome = None):
        self.df        = df 
        self.genome    = genome
        self.counter   = 0
        self.id_out    = []
        self.seq_out   = []
        self.fasta_out = []
    
    @staticmethod
    def inside(i, last_num):
        return i!=0 and i!=last_num-1
    @staticmethod
    def normal_interval(df_h, df_i):
        return max(int(df_h.start), int(df_h.end)) < min(int(df_i.start), int(df_i.end))
    @staticmethod
    def initialCDS(df_h, df_i):
        return df_h.sequence != df_i.sequence
    
    @classmethod
    def get_1st_interval(cls, df_i, 
                         genome_i, 
                         id_out = None, 
                         seq_out = None, 
                         counter = None):
        start = 1
        end = min(int(df_i.start), int(df_i.end))
        seq_out += [genome_i.seq[:end]] ;counter+=1
        id_out += [f'interval_region{counter}_LEN{end-start+1}_{genome_i.id}_{start}_{end}']
        return counter

    @classmethod
    def get_interval_inside(cls, df_h, df_i, 
                            genome_i, 
                            id_out = None, 
                            seq_out = None, 
                            counter = None):
        start = max(int(df_h.start), int(df_h.end))
        end = min(int(df_i.start), int(df_i.end))
        seq_out += [genome_i.seq[start-1:end]] ;counter+=1
        id_out += [f'interval_region{counter}_LEN{end-start+1}_{genome_i.id}_{start}_{end}']
        return counter

    @classmethod
    def get_last_interval(cls, df_h, 
                          genome_h, 
                          id_out = None,
                          seq_out = None, 
                          counter = None):
        start = max(int(df_h.start), int(df_h.end))
        end = len(str(genome_h.seq))
        seq_out += [genome_h.seq[end-1:]] ;counter+=1
        id_out += [f'interval_region{counter}_LEN{end-start+1}_{genome_h.id}_{start}_{end}']
        return counter

    def get_IS(self):
        for i in tqdm(range(len(self.df))):
            if i == 0:
                df_i = self.df.loc[i, ] ;tmpi = df_i.sequence
                
                genome_i = [genome_i for genome_i in self.genome if genome_i.id==tmpi][0]
                self.counter = MyGetIS.get_1st_interval(df_i, genome_i, 
                                                        id_out = self.id_out, 
                                                        seq_out = self.seq_out, 
                                                        counter = self.counter)

            elif MyGetIS.inside(i, len(self.df)):
                df_h = self.df.loc[i-1, ] ;tmph = df_h.sequence
                df_i = self.df.loc[i, ] ;tmpi = df_i.sequence

                genome_h = [genome_i for genome_i in self.genome if genome_i.id==tmph][0]
                genome_i = [genome_i for genome_i in self.genome if genome_i.id==tmpi][0]
                if MyGetIS.normal_interval(df_h, df_i):
                    if MyGetIS.initialCDS(df_h, df_i):
                        self.counter = MyGetIS.get_last_interval(df_h, genome_i, 
                                                                 id_out = self.id_out, 
                                                                 seq_out = self.seq_out, 
                                                                 counter = self.counter)
                        self.counter = MyGetIS.get_1st_interval(df_i, genome_i, 
                                                                id_out = self.id_out, 
                                                                seq_out = self.seq_out, 
                                                                counter = self.counter)

                    else:
                        self.counter = MyGetIS.get_interval_inside(df_h, df_i, genome_i, 
                                                                   id_out = self.id_out, 
                                                                   seq_out = self.seq_out, 
                                                                   counter = self.counter)

            elif i == len(self.df)-1:
                df_h = self.df.loc[i, ] ;tmph = df_h.sequence
                genome_h = [genome_i for genome_i in self.genome if genome_i.id==tmph][0]
                self.counter = MyGetIS.get_last_interval(df_h, genome_h, 
                                                         id_out = self.id_out, 
                                                         seq_out = self.seq_out, 
                                                         counter = self.counter)
        return self.id_out, self.seq_out
    
    def main(self, id_out, seq_out):
        for j in range(len(seq_out)):
            self.fasta_out += [SeqRecord(seq_out[j], 
                                         id=id_out[j], 
                                         description=id_out[j])]
        SeqIO.write(self.fasta_out, os.path.join(f'interval_regions.fasta'), "fasta")
        
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
    
    #assert os.path.isfile(''), "specify the path to your nr database"
    if 'res_ISblastx' not in os.listdir(path='./'):
        os.system('mkdir res_ISblastx')

    for fasta_path in tqdm(os.listdir(path=dir_in)):
        fasta = list(SeqIO.parse(f"./each_IS/{fasta_path}", "fasta"))[0]
        if len(fasta.seq) > threshold:
            os.system(f"blastx -query ./each_IS/{fasta_path} -out ./res_ISblastx/out_{fasta_path[:-6]}.txt -num_threads {num_threads} -evalue {evalue} -num_descriptions {num_descriptions} -num_alignments 100 -db {db}")
        else:
            pass

def main():
    df     = get_edited_features(path_to_features = get_args().features)
    genome = list(SeqIO.parse(get_args().genome, "fasta"))
    
    print('1.extracting IS..')
    instance = MyGetIS(df = df,
                       genome = genome)
    seq_out, id_out = instance.get_IS()
    instance.main(seq_out, id_out)
    split_fasta(multifasta ='./interval_regions.fasta')
    
    if not get_args().Without_blast:
        print('2.blastx now..')
        if get_args().mode == "strict":
            blastx(dir_in = './each_IS/',
                   num_threads = get_args().num_threads,
                   num_descriptions= get_args().num_descriptions,
                   db = get_args().DataBase, 
                   threshold = 300,
                   evalue =.0001 )
            
        elif get_args().mode == "loose":
            blastx(dir_in = './each_IS/',
                   num_threads = get_args().num_threads,
                   num_descriptions= get_args().num_descriptions,
                   db = get_args().DataBase, 
                   threshold =15,
                   evalue =.01 )

if __name__ == "__main__":
    main()
