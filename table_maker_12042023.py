import patsy
import statsmodels.api as sm
import pandas as pd
import numpy as np
from functools import reduce
from matplotlib import pyplot as plt
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("white")
import warnings
import sys
warnings.filterwarnings("ignore")
pd.set_option("display.max_rows", None)

# big table
#table_dir='/work/carrot-zhang/data_matrix/big_table_112823/solidheme_clinicogenomic_oncokb.tsv'
table_dir='/work/carrot-zhang/data_matrix/big_table_010924/solidheme_clinicogenomic_oncokb.tsv'
df_table=pd.read_csv(f'{table_dir}', delimiter="\t", low_memory=False)
print('table loaded')
for i in range(len(df_table.columns)):
   print(f'index {i}: {df_table.columns[i]}')

#print(df_table["AVERAGE_BMI"].head(20))
#print(df_table["YOST_INDEX_IMPUTED"].head(20))

# Chris's dataset
text_dir='/work/carrot-zhang/david/aacr_2023_ancestry_io_data_sdh.tsv'
df_text = pd.read_csv(f'{text_dir}', delimiter="\t", low_memory=False)

# hotspot mutation file
#hotspot_dir='/work/carrot-zhang/david/TP53_mutations.txt'
#df_hotspot = pd.read_csv(f'{hotspot_dir}', delimiter="\t", low_memory=False)

#res = [i for i in df_final.columns.tolist() if 'SMOKING' in i]
#print(res)

print(df_table.shape)
print(df_text.shape)

print('merging now')

# merging
merged_table=pd.merge(df_table, df_text[['SAMPLE_ID','NCI_SCORE']], left_on='SAMPLE_ID', right_on='SAMPLE_ID', how='inner')
print(merged_table.shape)
#merged_table=pd.merge(merged_table, df_hotspot, left_on='SAMPLE_ID', right_on='SAMPLE_ID', how='left')
#print(merged_table.shape)

res = [i for i in merged_table.columns.tolist() if 'SMOKING' in i]
print(res)

#merged_table['TP53_R175'] = np.where(merged_table['TP53_y'].astype(str).str.contains("R175"), 1, 0)
#merged_table['TP53_G245'] = np.where(merged_table['TP53_y'].astype(str).str.contains("G245"), 1, 0)
#merged_table['TP53_R248'] = np.where(merged_table['TP53_y'].astype(str).str.contains("R248"), 1, 0)
#merged_table['TP53_R273'] = np.where(merged_table['TP53_y'].astype(str).str.contains("R273"), 1, 0)
#merged_table['TP53_hotspot'] = np.where(merged_table['TP53_y'].astype(str).str.contains("R175|G245|R248|R273"), 1, 0)

print(merged_table.shape)
merged_table=merged_table.drop_duplicates(subset=['SAMPLE_ID'])
print('dropping sample id: ',merged_table.shape)
merged_table=merged_table.drop_duplicates(subset=['PATIENT_ID'])
print('dropping patient id: ',merged_table.shape)
merged_table.to_csv('/work/carrot-zhang/david/ancestry_datatable_01092024.csv')
