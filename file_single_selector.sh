#!bin/bash

file = /work/carrot-zhang/tomin/for_jian/analysis_tables_030823/bladder_data_subtypes.tsv
echo ${file}
mkdir -p /work/carrot-zhang/david/ancestry_manuscript/0308
python /work/carrot-zhang/david/ancestry_script/ancestry_$(basename ${file} .tsv).py ${file} > /work/carrot-zhang/david/ancestry_manuscript/0308/$(basename ${file} .tsv)_result_summarys.txt   

