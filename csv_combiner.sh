#!bin/bash

mkdir -p /work/carrot-zhang/david/local_ancestry/results/all

chr_list=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22")

for chr in ${chr_list[@]}
do
   destination=/work/carrot-zhang/david/local_ancestry/results/all/tp53_${chr}_coeff_zscore_table.csv 
   > ${destination}

   files=$(ls /work/carrot-zhang/david/local_ancestry/results/${chr}/TP53*_zscore.csv)
   echo combining files for ${chr}

   counter=0
   for file in $files
   do
      if [ ${counter} == 0 ]
      then
         cat ${file} >> ${destination}
         counter=1
      else
         cat ${file} | tail -n +2 >> ${destination}
      fi
   done
done

#cancer=("bladder" "breast" "colorectal" "crosscancer" "esophagus" 
#               "glioma" "liver" "lung" "ovary" "pancreas" "prostate" 
#                "renal" "skin" "uterus")

#SEF=("wSEF" "withoutSEF")

#for type in ${cancer[@]}
#do
#   for included in ${SEF[@]}
#   do
#      cat /work/carrot-zhang/david/ancestry_manuscript_tables/tert_${type}_table_${included}.csv >> ${destination}
#   done
#done
