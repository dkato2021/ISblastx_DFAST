#!/bin/sh
# get the top hit genes id from the results of blast
#$1:que_bin <- res_blastx

resdir=`basename $1`
mkdir geneid_"$resdir"
for file in $1/* ;do
    sed -n '22p' $file | cut -d" " -f3 >> geneid_"$resdir"/geneid_"$resdir".txt
    #sed -n '23p' $file | cut -d" " -f3 >> geneid_"$resdir"/geneid_"$resdir".txt
    #sed -n '24p' $file | cut -d" " -f3 >> geneid_"$resdir"/geneid_"$resdir".txt
done