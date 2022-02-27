# blastxis.py

## Installation
```
$ git clone https://github.com/dkato2021/blastxis.git
$ chmod u+x blastxis.py
```
## Usage
```
$ blastxis.py -g genome.fna

#if you omit DFAST
$ blastxis.py -g genome.fna -f featues.tsv
```

## optional arguments
```
optional arguments:
  -h, --help            show this help message and exit
  -g GENOME             path to genome.fasta from DFAST. The contig ID is
                        unified in the form of "sequence~", with features-
                        option
  -f FEATURES           path to features.tsv from DFAST.
  -e EVALUE, --evalue EVALUE
                        evalue in blastx.(default:0.001)
  -th THRESHOLD, --threshold THRESHOLD
                        minimum length of interval sequence as input of
                        blastx.(default:0)
  -db DATABASE, --database DATABASE
                        path to your nr database.(default:/home_ssd/local/db/b
                        lastdb.20200904/nr)
  -t NUM_THREADS, --num_threads NUM_THREADS
                        num threads in blastx.(default:48)
  -cpu NUM_THREADS_DFAST, --num_threads_dfast NUM_THREADS_DFAST
                        num threads in dfast.(default:48)
  -nd NUM_DESCRIPTIONS, --num_descriptions NUM_DESCRIPTIONS
                        num descriptions in blastx.(default:50)
  -x {True,False}       If "True", blastx will not be executed.(default:False)
```
[DFAST](https://dfast.ddbj.nig.ac.jp "DFAST Home")
