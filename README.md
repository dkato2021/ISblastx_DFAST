# ISblastx_DFAST.py
**DFASTが出力するデータを継承してCDS間領域をblastxするスクリプト**
## Usage
```
$ python3 ISblastx_DFAST.py -g genome_sample.fasta -f featues_sample.tsv
```
## 入出力と依存
- 入力データ
  - DFASTが出力したgenome.fasta. (contig IDが'sequence~'の形で統一されているもの)
  - DFASTが出力したfeatures.tsv. (CDSの開始点と終了点およびそのCDSが抽出されたcontig IDを使用)
- 出力データ
  - CDS間領域の塩基配列データ (-> interval_regions.fasta)
  - CDS間領域のblastxの結果
- 依存
  - python 3.7

## 性能の評価
`２つの作業、つまりCDSの開始点と終了点の精査およびCDS間領域における新CDSの探索について、t吉岡氏がマニュアルで行なった結果と本スクリプトを使用して行なった結果の比較を行なった`  

- CDS領域の開始点と終了点を精査する作業について
  - 正解率(accuracy) : 93.77%  
  - 適合率(precision) : 48.89%  
  - 再現率(recall) : 32.35%  
  - F1スコア : 38.94%  
  - 混合行列↓  
![](./images/A_ver2.png)
(0, 1) = (変更なし, CDS領域を変更)  

  ``comment=> t吉岡氏が68のCDS領域に対して開始点変更を行ったのに対して、本スクリプトを使用した結果は出会った。今回の実験ではblastxに使用する配列の最低塩基長として50bpという閾値を使用したが、この閾値の値を下げることでより多く開始点の変更が必要であるCDSを検出できる可能性がある。``  

- CDS間領域における新CDSの探索について
  - 正解率(accuracy) : 96.57%  
  - 適合率(precision) : 33.3%  
  - 再現率(recall) : 52.00%  
  - F1スコア : 40.63%  
  - 混合行列↓  
![](./images/B.png)
(0, 1) = (変更なし, 新CDSを発見)  

  ``comment=> t吉岡氏のアノテーションと比べて新たなCDSを多く検出することができた。``  

``実用可能かどうかについて言及する``

## 影響
- CDS領域の開始点と終了点の精査について  
本スクリプトを使用すると開始点と終了点を変更する必要があるCDSを一度に列挙することができ、DFASTが抽出したCDSの前後の領域をgeneticsやexpacyなどを使用して参照する作業の回数が減ると考えられる。

- CDS間領域における新CDSの探索について  
特になし


## optional arguments
```
  -h, --help            show this help message and exit
  -g GENOME             path to genome.fasta from DFAST. The contig ID is
                        unified in the form of "sequence~"
  -f FEATURES, --features FEATURES
                        path to features.tsv from DFAST
  -e EVALUE, --evalue EVALUE
                        evalue in blastx.(default:0.0001)
  -th THRESHOLD, --threshold THRESHOLD
                        minimum length of IS sequence as input of
                        blastx.(default:50)
  -db DATABASE, --database DATABASE
                        path to your nr database.(default:/home_ssd/local/db/b
                        lastdb.20200904/nr)
  -t NUM_THREADS, --num_threads NUM_THREADS
                        num_threads in blastx.(default:16)
  -nd NUM_DESCRIPTIONS, --num_descriptions NUM_DESCRIPTIONS
                        num_descriptions in blastx.(default:50)
  --Without_blast {True,False}
                        If "True", blastx will not be executed.(default:False)
```

## 補足事項
- CDS間領域が数十kbpほど抽出されていない可能性があるため、原因が分かり次第プログラムを修正する予定である。


