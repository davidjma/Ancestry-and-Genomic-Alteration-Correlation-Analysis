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
import os
import sklearn
from sklearn.model_selection import KFold # import KFold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeavePOut
from operator import add
import scipy.stats
from statistics import median

# defining variables
chr_list=['chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']

# defining table variables
table_dir=f'/work/carrot-zhang/david/local_ancestry/tables/ancestry_wLocalAncestry_chr1_table.csv'
#df_final=pd.read_csv(f'{table_dir}', delimiter=",", low_memory=False)
#print('Initial chr1 shape', df_final.shape)
chrN=sys.argv[4]
kf = KFold(n_splits=int(sys.argv[3]))
home='/work/carrot-zhang/david/local_ancestry/'
gene=sys.argv[2]
cancertype='All'

# merging all chromesome tables
#for num in chr_list:

   # defining variables
#   table_dir=f'/work/carrot-zhang/david/local_ancestry/tables/ancestry_wLocalAncestry_{num}_table.csv'
#   df_chr=pd.read_csv(f'{table_dir}', delimiter=",", low_memory=False)

   # combining all columns with snps for each chromosome table
#   keyword='snp'
#   keylist = [i for i in df_chr.columns if keyword in i]
#   keylist.append('SAMPLE_ID')
#   print(f'For chromosome {num}, there are {len(keylist)} snps')
#   df_final=pd.merge(df_final, df_chr[keylist], how="inner",left_on='SAMPLE_ID',right_on='SAMPLE_ID')
#   keylist = [i for i in df_chr.columns if keyword in i]
#   print(f'Merged shape is {df_final.shape}')

# This table has snp columns for all choromosomes for each sample_id
#df_final.to_csv(f'/work/carrot-zhang/david/{sys.argv[1]}_combined_table.csv')
#print(df_final.shape)

df_final=pd.read_csv(f'/work/carrot-zhang/david/{sys.argv[1]}_combined_table.csv')

# finding snp columns to compute local ancestry
subs='snp'
res = [i for i in df_final.columns if subs in i]
snp_loc=res
#print(snp_loc)
print(len(snp_loc))

# defining paths
txtpath = '/work/carrot-zhang/david/local_ancestry/txtResults/kfold' 
csvpath='/work/carrot-zhang/david/local_ancestry/results/kfold'

# making directories
if not os.path.exists(csvpath):
    os.makedirs(csvpath)

if not os.path.exists(f'{csvpath}/{chrN}/{sys.argv[1]}'):
   os.makedirs(f'{csvpath}/{chrN}/{sys.argv[1]}')

if not os.path.exists(f'{csvpath}/{chrN}/results/{sys.argv[1]}'):
   os.makedirs(f'{csvpath}/{chrN}/results/{sys.argv[1]}')

if not os.path.exists(txtpath):
    os.makedirs(txtpath)

if not os.path.exists(f'{txtpath}/{chrN}/{sys.argv[1]}'):
   os.makedirs(f'{txtpath}/{chrN}/{sys.argv[1]}')

#cancer_condition=[
#df_final['CANCER_TYPE'] == 'Endometrial Cancer',
#df_final['CANCER_TYPE'] == 'Hepatobiliary Cancer',
#df_final['CANCER_TYPE'] == 'Non-Small Cell Lung Cancer',
#df_final['CANCER_TYPE'] == 'Ovarian Cancer',
#df_final['CANCER_TYPE'] == 'Pancreatic Cancer',
#df_final['CANCER_TYPE'] == 'Breast Cancer',
#df_final['CANCER_TYPE'] == 'Colorectal Cancer',
#df_final['CANCER_TYPE'] == 'Prostate Cancer',
#df_final['CANCER_TYPE'] == 'Esophagogastric Cancer',
#df_final['CANCER_TYPE'] == 'Glioma'
#]

#df_final=df_final.loc[np.bitwise_or.reduce(cancer_condition)]

#condition=[
#df_final['CANCER_TYPE_DETAILED'] =='Endometrial Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor',
#df_final['CANCER_TYPE_DETAILED'] =='Uterine Clear Cell Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Uterine Endometrioid Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Uterine Mixed Endometrial Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Hepatocellular Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Lung Adenocarcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Lung Squamous Cell Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Clear Cell Ovarian Cancer',
#df_final['CANCER_TYPE_DETAILED'] =='Endometrioid Ovarian Cancer',
#df_final['CANCER_TYPE_DETAILED'] =='High-Grade Serous Ovarian Cancer',
#df_final['CANCER_TYPE_DETAILED'] =='Ovarian Epithelial Tumor',
#df_final['CANCER_TYPE_DETAILED'] =='Serous Ovarian Cancer',
#df_final['CANCER_TYPE_DETAILED'] =='Pancreatic Adenocarcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Breast Invasive Cancer, NOS',
#df_final['CANCER_TYPE_DETAILED'] =='Breast Invasive Carcinoma, NOS',
#df_final['CANCER_TYPE_DETAILED'] =='Breast Invasive Ductal Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Breast Invasive Lobular Carcinoma',
#df_final['CANCER_TYPE_DETAILED'] =='Invasive Breast Carcinoma'
#]

#df_final=df_final.loc[np.bitwise_or.reduce(condition)]

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
race_list = ["AFR"] # "EAS", "EUR", "SAS", "ASJ", "NAM"]

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

#print(df_final['CANCER_TYPE'].value_counts())

# dropping empty columns
df_final=df_final.dropna(subset=['YOST_INDEX_STATUS'])
df_final=df_final.dropna(subset=['BMI_STATUS'])
#df_final=df_final.dropna(subset=['SMOKING_STATUS_PREDICTED'])
df_final=df_final.dropna(subset=['SEX'])
df_final=df_final.dropna(subset=['PATIENT_CURRENT_AGE'])
df_final=df_final.dropna(subset=['GENE_PANEL'])
df_final=df_final.dropna(subset=['HAS_STAGE_IV_DX'])
df_final=df_final.dropna(subset=['MOST_RECENT_INSURANCE_CATEGORY'])
df_final=df_final.dropna(subset=['NCI_SCORE'])

# performing the logistic regression cancer_details is missing
print('logistic regression beginning...')

P={}
C={}
Z={}

count=0
together=pd.DataFrame()
results_output=f'/work/carrot-zhang/david/local_ancestry/results/kfold/{chrN}/results/{sys.argv[1]}/{chrN}_{sys.argv[1]}_results_kfold_GC.csv'
combine=pd.DataFrame()

df_final['risk_sum']=0

chisq1=[]
chisq2=[]
chisq3=[]
chisq4=[]
chisq5=[]
chisq6=[]
chisq7=[]
chisq8=[]
chisq9=[]
chisq10=[]

Z1=[]
Z2=[]
Z3=[]
Z4=[]
Z5=[]
Z6=[]
Z7=[]
Z8=[]
Z9=[]
Z10=[]

df1=[]
df2=[]
df3=[]
df4=[]
df5=[]
df6=[]
df7=[]
df8=[]
df9=[]
df10=[]

# find risk sum 
for snp in snp_loc:
   
   count+=1
   print(f'This {snp} is being analyzed.....{count}/{len(snp_loc)} in progress')
   foldN=0
   for train_index, test_index in kf.split(df_final[gene].values):
      train=df_final.iloc[train_index]
      test=df_final.iloc[test_index]

      if count==1:
         test["risk_sum"]=0
         print('------running analysis----------')

      print(f'training set size: {len(train)}, test set size: {len(test)}')

      foldN+=1

      race='AFR'
      f= gene+ f"~ {race} + {snp} + PATIENT_CURRENT_AGE + YOST_INDEX_STATUS + SEX + GENE_PANEL + HAS_STAGE_IV_DX + NCI_SCORE + MOST_RECENT_INSURANCE_CATEGORY + BMI_STATUS" 
      print(f"{race} is being analyzed......")

      y, X = patsy.dmatrices(f, train, return_type='dataframe') 
      results = sm.Logit(y, X).fit()
      
      P=results.pvalues[11:12].values
      C=results.params[11:12].values
      Z=results.tvalues[11:12].values

      print(f'{foldN}: Z-score={Z}')
    
      if foldN==1:
         chisq=Z[0]**2
#         Z1.append(Z[0])
         chisq1.append(chisq)

      elif foldN==2:
         chisq=Z[0]**2
#         Z2.append(Z[0])
         chisq2.append(chisq)

      elif foldN==3:
         chisq=Z[0]**2
#         Z3.append(Z[0])
         chisq3.append(chisq)
         
      elif foldN==4:
         chisq=Z[0]**2
#         Z4.append(Z[0])
         chisq4.append(chisq)

      elif foldN==5:
         chisq=Z[0]**2
#         Z5.append(Z[0])
         chisq5.append(chisq)

      elif foldN==6:
         chisq=Z[0]**2
#         Z6.append(Z[0])
         chisq6.append(chisq)

      elif foldN==7:
         chisq=Z[0]**2
#         Z7.append(Z[0])
         chisq7.append(chisq)

      elif foldN==8:
         chisq=Z[0]**2
#         Z8.append(Z[0])
         chisq8.append(chisq)

      elif foldN==9:
         chisq=Z[0]**2
#         Z9.append(Z[0])
         chisq9.append(chisq)

      elif foldN==10:
         chisq=Z[0]**2
#         Z10.append(Z[0])
         chisq10.append(chisq)
     
lbda1=median(chisq1)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda2=median(chisq2)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda3=median(chisq3)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda4=median(chisq4)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda5=median(chisq5)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda6=median(chisq6)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda7=median(chisq7)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda8=median(chisq8)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda9=median(chisq9)/scipy.stats.chi2.ppf(1-.05, df=1)
lbda10=median(chisq10)/scipy.stats.chi2.ppf(1-.05, df=1)

for snp in snp_loc:
   count+=1
   print(f'This {snp} is being analyzed for risk score calculation.....{count}/{len(snp_loc)} in progress')

   foldN=0
   for train_index, test_index in kf.split(df_final[gene].values):
      train=df_final.iloc[train_index]
      test=df_final.iloc[test_index]

      if count==1:
         test["risk_sum"]=0
         print('------running analysis----------')

      print(f'training set size: {len(train)}, test set size: {len(test)}')

      foldN+=1

      race='AFR'
      f= gene+ f"~ {race} + {snp} + PATIENT_CURRENT_AGE + YOST_INDEX_STATUS + SEX + GENE_PANEL + HAS_STAGE_IV_DX + NCI_SCORE + MOST_RECENT_INSURANCE_CATEGORY + BMI_STATUS"
      print(f"{race} is being analyzed......")

      y1, X1 = patsy.dmatrices(f, train, return_type='dataframe')
      results = sm.Logit(y, X).fit()

      P=results.pvalues[11:12].values
      C=results.params[11:12].values
      Z=results.tvalues[11:12].values

      print(f'{foldN}: Z-score={Z}')

      if foldN==1:
         Z1.append(Z[0])
         lbda=lbda1

      elif foldN==2:
         Z2.append(Z[0])
         lbda=lbda2

      elif foldN==3:
         Z3.append(Z[0])
         lbda=lbda3

      elif foldN==4:
         Z4.append(Z[0])
         lbda=lbda4

      elif foldN==5:
         Z5.append(Z[0])
         lbda=lbda5

      elif foldN==6:
         Z6.append(Z[0])
         lbda=lbda6

      elif foldN==7:
         Z7.append(Z[0])
         lbda=lbda7

      elif foldN==8:
         Z8.append(Z[0])
         lbda=lbda8

      elif foldN==9:
         Z9.append(Z[0])
         lbda=lbda9

      elif foldN==10:
         Z10.append(Z[0])
         lbda=lbda10

      initial=df_final.loc[test['risk_sum'].index.tolist(), "risk_sum"]
      new=initial+test[snp]*(Z/lbda)
      df_final.loc[new.index.tolist(),'risk_sum']=new.tolist()

df_final.to_csv(f'/work/carrot-zhang/david/{sys.argv[1]}_ancestry_table.csv')
df_final=pd.read_csv(f'/work/carrot-zhang/david/{sys.argv[1]}_ancestry_table.csv', delimiter=",", low_memory=False)

f= gene+ f"~ {race} + risk_sum + PATIENT_CURRENT_AGE + YOST_INDEX_STATUS + SEX + GENE_PANEL + HAS_STAGE_IV_DX + NCI_SCORE + MOST_RECENT_INSURANCE_CATEGORY + BMI_STATUS"
y, X = patsy.dmatrices(f, df_final, return_type='dataframe')
results = sm.Logit(y, X).fit()

coef=pd.DataFrame(results.params, columns=['coefficient']) #regression coeff column
pval=pd.DataFrame(results.pvalues, columns=['p-values']) #p-values column
tval=pd.DataFrame(results.tvalues, columns=['z-score']) #z-score column
conf=pd.DataFrame(results.conf_int(0.05)) #conf interval column
conf['chrN']=chrN
conf['covariates']=conf.index

individual=pd.concat([coef, pval, tval, conf], axis=1)

individual=individual.rename(columns={ 0: 'CI 0.025' })
individual=individual.rename(columns={ 1: 'CI 0.975' })

individual = individual[['chrN','covariates', 'coefficient', 'p-values', 'z-score', 'CI 0.025', 'CI 0.975']]
individual = individual.reset_index(drop=True)
individual.to_csv(results_output)
