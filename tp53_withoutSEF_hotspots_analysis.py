# importing libraries
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

# table directories 
table_dir='/work/carrot-zhang/david/tp53_table.csv'
merged_table_200=pd.read_csv(f'{table_dir}', delimiter=",", low_memory=False)

#text_dir='/work/carrot-zhang/david/aacr_2023_ancestry_io_data_sdh.tsv'
#df_text = pd.read_csv(f'{text_dir}', delimiter="\t", low_memory=False)

#merged_table_200=pd.merge(merged_table_200, df_text[['SAMPLE_ID','MOST_RECENT_INSURANCE_CATEGORY','NCI_SCORE']], left_on='SAMPLE_ID', right_on='SAMPLE_ID', how='left')
#merged_table_200.to_csv('/work/carrot-zhang/david/tp53_table.csv')

#merged_table_200.to_parquet('/work/carrot-zhang/david/tp53_table.parquet',engine='pyarrow')

merged_table_200['TP53_hotspot'] = np.where(merged_table_200['TP53_y'].astype(str).str.contains("R175|G245|R248|R273"), 1, 0)

print(merged_table_200['TP53_hotspot'].value_counts())
print(merged_table_200['TP53_R175'].value_counts())
print(merged_table_200['TP53_R273'].value_counts())
print(merged_table_200['TP53_G245'].value_counts())
print(merged_table_200['TP53_R248'].value_counts())

# defining variables
home='/work/carrot-zhang/david/ancestry_manuscript/0308'
filename=home+'/'+'tp53_withoutSEF_only_hotspot_analysis_result_summarys.txt'
gene="TP53_hotspot"
df_final=merged_table_200

# turning yost index into a categorical variable
df_final['YOST_INDEX_STATUS']=np.nan
df_final['YOST_INDEX_STATUS'][df_final['YOST_INDEX_IMPUTED']>=37]="high"
df_final['YOST_INDEX_STATUS'][df_final['YOST_INDEX_IMPUTED']<37]="low"

# turning bmi into a categorical variable
df_final['BMI_STATUS']=np.nan
df_final['BMI_STATUS'][df_final['AVERAGE_BMI']<18.5]="underweight"
df_final['BMI_STATUS'][(df_final['AVERAGE_BMI']>=18.5) & (df_final['AVERAGE_BMI']<24.9)]="Healthy"
df_final['BMI_STATUS'][(df_final['AVERAGE_BMI']>=24.9) & (df_final['AVERAGE_BMI']<29.9)]="overweight"
df_final['BMI_STATUS'][df_final['AVERAGE_BMI']>=30.0]="obese"

#xdisease stage, BMI, insurance, comorbidities
#print(df_final['HAS_STAGE_IV_DX'])
#print(df_final['BMI_STATUS'].value_counts())

#print(df_final['MOST_RECENT_INSURANCE_CATEGORY'].value_counts())
#print(df_final['NCI_SCORE'].value_counts())

#print(df_final['SMOKING_STATUS_PREDICTED'].value_counts())
#print(df_final['SMOKING_STATUS'].value_counts())

# defining race list
race_list = ["AFR", "EAS", "EUR", "SAS", "ASJ", "SAS"]

#print(df_final[gene].value_counts())
#print(df_final['SMOKING_STATUS'].value_counts(dropna=False))
#print(df_final['GENE_PANEL'].value_counts())
#print(df_final['CANCER_TYPE_DETAILED'].value_counts())
#print(df_fina;['AFR'].value_counts(dropna=False))
#print(df_final['YOST_INDEX_STATUS'].value_counts())
#print(df_final['PATIENT_CURRENT_AGE'].value_counts())
#print(df_final['SMOKING_STATUS'].value_counts())

#df_final.rename({'INSURANCE_(MOST_RECENT)': 'INSURANCE'}, axis=1, inplace=True)
#print(df_final['SMOKING_STATUS'].value_counts())
#print(df_final['HAS_STAGE_IV_DX'].value_counts())
#print(df_final['INSURANCE'].value_counts())
#print(df_final['YOST_INDEX_STATUS'].value_counts())
#print(df_final['BMI_STATUS'].value_counts())
#print(df_final['YOST_INDEX_IMPUTED'].value_counts())

#df_final.to_csv('/work/carrot-zhang/david/tp53_table.csv')

# dropping empty columns
#df_final=df_final.dropna(subset=['YOST_INDEX_STATUS'])
#df_final=df_final.dropna(subset=['BMI_STATUS'])
#df_final=df_final.dropna(subset=['SMOKING_STATUS'])
df_final=df_final.dropna(subset=['SEX'])
df_final=df_final.dropna(subset=['PATIENT_CURRENT_AGE'])
df_final=df_final.dropna(subset=['GENE_PANEL'])
df_final=df_final.dropna(subset=['HAS_STAGE_IV_DX'])
df_final=df_final.dropna(subset=['CANCER_TYPE_DETAILED'])

#df_final=df_final.dropna(subset=['MOST_RECENT_INSURANCE_CATEGORY'])
#df_final=df_final.dropna(subset=['NCI_SCORE'])

#filename_details=home+'/'+'tp53_CancerDetails_counts.txt'
#details = open(filename_details,'w')

#print(df_final['CANCER_TYPE_DETAILED'].value_counts(), file=details)

# including top cancer types detailed 
df_final=df_final.loc[(df_final['CANCER_TYPE_DETAILED'] == 'Lung Adenocarcinoma') | (df_final['CANCER_TYPE_DETAILED'] == 'Small Cell Lung Cancer')|
(df_final['CANCER_TYPE_DETAILED'] == 'Pleural Mesothelioma, Epithelioid Type') | (df_final['CANCER_TYPE_DETAILED'] == 'Non-Small Cell Lung Cancer') | (df_final['CANCER_TYPE_DETAILED'] == 'Large Cell Neuroendocrine Carcinoma') | (df_final['CANCER_TYPE_DETAILED'] == 'Pleural Mesothelioma') |
(df_final['CANCER_TYPE_DETAILED'] == 'Lung Squamous Cell Carcinoma') | (df_final['CANCER_TYPE_DETAILED'] == 'Combined Small Cell Lung Carcinoma') | (df_final['CANCER_TYPE_DETAILED'] == 'Poorly Differentiated Non-Small Cell Lung Cancer') | (df_final['CANCER_TYPE_DETAILED'] == 'Lung Adenosquamous Carcinoma')]

#print(df_final.loc[df_final['CANCER_TYPE_DETAILED'] == 'Lung Adenocarcinoma|Small Cell Lung Cancer'])
#df=df_all.dropna(subset=['SMOKING_STATUS'])

#print(df_final.columns.tolist()[:300])
#res = [i for i in df_final.columns.tolist() if 'TP53' in i]
#print(res)
#print(df_final.columns.tolist().str.contains('INSURANCE'))
#res = [i for i in df_final.columns.tolist() if 'insurance' in i]
#print(res)
#print(df_final)

#print(res)
file = open(filename,'w')

# performing the logistic regression
print('logistic regression beginning...')
for race in race_list:
   f= gene+ f"~ {race} + PATIENT_CURRENT_AGE + SEX + GENE_PANEL + HAS_STAGE_IV_DX + CANCER_TYPE_DETAILED" 
   print(f"{race} is being analyzed")
 
   y, X = patsy.dmatrices(f, df_final, return_type='dataframe')
  
   result = sm.Logit(y, X)
   res=result.fit()
   print(f"\n\n\n---------------Table for {race}---------------------", file=file)
   print(res.summary(), file=file)
