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

chrN=sys.argv[6]

# table directories 
table_dir=f'/home/mad1/ancestry_project/tables/{sys.argv[1]}_LocalAncestry_{chrN}_table.csv'
df_final=pd.read_csv(f'{table_dir}', delimiter=",", low_memory=False)
kf = KFold(n_splits=int(sys.argv[5]))

subs='snp'
res = [i for i in df_final.columns if subs in i]

# defining variables
txtpath = '/work/carrot-zhang/david/local_ancestry/txtResults/kfold' 
csvpath='/work/carrot-zhang/david/local_ancestry/results/kfold'

if not os.path.exists(csvpath):
    os.makedirs(csvpath)

if not os.path.exists(f'{csvpath}/{chrN}'):
   os.makedirs(f'{csvpath}/{chrN}')

if not os.path.exists(f'{csvpath}/{chrN}/riskscore'):
   os.makedirs(f'{csvpath}/{chrN}/riskscore')

if not os.path.exists(txtpath):
    os.makedirs(txtpath)

if not os.path.exists(f'{txtpath}/{chrN}'):
   os.makedirs(f'{txtpath}/{chrN}')

home='/work/carrot-zhang/david/local_ancestry/'
gene=sys.argv[2]
cancertype='All'

snp_loc=res

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
#df_final=df_final.dropna(subset=['SMOKING_STATUS'])
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

#print(df_final[snp_loc])
#print(df_final[snp_loc].value_counts())

#print(df_final)

count=0
together=pd.DataFrame()
risk_output=f'/work/carrot-zhang/david/local_ancestry/results/kfold/{chrN}/riskscore/{chrN}_riskscore_kfold.csv'

for snp in snp_loc:
   filename=home+'/txtResults/kfold/'+f'{chrN}/'+f'{sys.argv[3]}_All_LocalAncestry_summarys_{sys.argv[4]}_{chrN}_{snp}_kfold.txt'
   file=open(filename, 'w')

   count+=1
   print(f'This {snp} is being analyzed.....{count}/{len(snp_loc)} in progress')
   
   foldN=0
   for train_index, test_index in kf.split(df_final[gene].values):
      train=df_final.iloc[train_index]
      test=df_final.iloc[test_index]

      test["risk_sum"]=0
      print(f'training set size: {len(train)}, test set size: {len(test)}')

      foldN+=1

#      table_output=f'/work/carrot-zhang/david/local_ancestry/results/{sys.argv[12]}/{sys.argv[3]}_All_LocalAncestry_{sys.argv[4]}_{sys.argv[12]}_{snp}_table_zscore_kfold.csv'

      race='AFR'
      f= gene+ f"~ {race} + {snp} + PATIENT_CURRENT_AGE + YOST_INDEX_STATUS + SEX + GENE_PANEL + HAS_STAGE_IV_DX + NCI_SCORE + MOST_RECENT_INSURANCE_CATEGORY + BMI_STATUS" 
      print(f"{race} is being analyzed......")

      y, X = patsy.dmatrices(f, df_final, return_type='dataframe')
  
      results = sm.Logit(y, X).fit()
      
      P=results.pvalues[11:12].values
      C=results.params[11:12].values
      Z=results.tvalues[11:12].values

#      print(results.pvalues)
#      print('showing......', P)

#      print(results.tvalues)
#      print('showing......', Z)

#      if (float(P)<0.00005) and (float(C)>0):
#         test["risk_sum0"]=test["risk_sum0"]+test[snp_loc].astype(int)*Z.astype(int)
#         print('cond 1')
                        
#      if (float(P)<0.0001) and (float(C)>0):
#         test["risk_sum1"]=test["risk_sum1"]+test[snp_loc].astype(int)*Z.astype(int)
#         print('cond 2')

#      if (float(P)<0.001) and (float(C)>0):
#         test["risk_sum2"]=test["risk_sum2"]+test[snp_loc].astype(int)*Z.astype(int)
#         print('cond 3')

#      if (float(P)<0.05) and (float(C)>0):
#         test["risk_sum3"]=test["risk_sum3"]+test[snp_loc].astype(int)*Z.astype(int)
#         print('cond 4')

#      else:
#       print('cond 5')

      test["risk_sum"]=test["risk_sum"]+test[snp]*Z
         
#      print(test["HAS_STAGE_IV_DX"].value_counts())
#      print(test["GENE_PANEL"].value_counts())
#      print(test["MOST_RECENT_INSURANCE_CATEGORY"].value_counts())	
#      print(test["BMI_STATUS"].value_counts())
#      print(test["SEX"].value_counts())
#      print(test["risk_sum"].value_counts())
#      print(test[snp].value_counts())
#      print(test[race].value_counts())

	# MOST_RECENT_INSURANCE_CATEGORY {snp} + BMI_STATUS	
#      f = gene+ f"~ {race} + risk_sum + PATIENT_CURRENT_AGE + YOST_INDEX_STATUS + SEX + GENE_PANEL + HAS_STAGE_IV_DX + NCI_SCORE + MOST_RECENT_INSURANCE_CATEGORY + BMI_STATUS"
#      y, X = patsy.dmatrices(f, test, return_type='dataframe')
#      result = sm.Logit(y, X).fit()

  # f = gene+ f"~ {race} + {snp_loc}  + PATIENT_CURRENT_AGE + YOST_INDEX_STATUS + SEX + GENE_PANEL + HAS_STAGE_IV_DX + NCI_SCORE + MOST_RECENT_INSURANCE_CATEGORY + BMI_STATUS"
  # y, X = patsy.dmatrices(f, test, return_type='dataframe')
  # result = sm.Logit(y, X).fit()

#      print(f"\n\n\n---------------Table for fold {foldN}---------------------", file=file)
#      print (result.summary(), file=file)

#      print(f"\n\n\n---------------Table for {race}---------------------", file=file)
#      print(res.summary(), file=file)

#         ancestry=pd.DataFrame([race]*len(res.pvalues), columns=['Race']) # Race column
#         type_cancer=pd.DataFrame([cancertype]*len(res.pvalues), columns=['Cancer_type']) #cancer type column

         #snp_name=pd.DataFrame(snp, columns=['p-values']) #p-values column

      #riskscore=pd.DataFrame(test["risk_sum"], columns=['Risk_Score']) #z-score column

      snp_name=pd.DataFrame([snp]*len(test["risk_sum"]), columns=['snpName'])
      test.reset_index(inplace = True)
      chrName=pd.DataFrame([chrN]*len(test["risk_sum"]), columns=['chrN'])

#      print([snp]*len(test["risk_sum"]))
#      print(test["risk_sum"])
#      print(test["SAMPLE_ID"])

#      print(snp_name)
#      print(chrName)

      folds_df=pd.concat([chrName, snp_name, test["SAMPLE_ID"], test["risk_sum"]], axis=1)
      together=pd.concat([together, folds_df], axis=0)
    #  together.reset_index(inplace=True)
      together.to_csv(risk_output)

#         coef=pd.DataFrame(res.params, columns=['coefficient']) #regression coeff column
#         pval=pd.DataFrame(res.pvalues, columns=['p-values']) #p-values column
#         tval=pd.DataFrame(res.tvalues, columns=['z-score']) #z-score column
#         conf=pd.DataFrame(res.conf_int(0.05)) #conf interval column

#         ancestry['row_number'] = ancestry.reset_index().index
#         type_cancer['row_number'] = type_cancer.reset_index().index
#         conf['row_number'] = np.arange(len(conf))
#         conf['covariates']=conf.index

#         individual=pd.concat([coef, pval, tval, conf], axis=1)

#         together=pd.merge(individual,ancestry,left_on='row_number',right_on='row_number', how='left')
#         together=pd.merge(together,type_cancer,left_on='row_number',right_on='row_number', how='left')

#        print(combine) 
#        print(together)
#        print(combine.columns)
#        print(together.columns)

#         combine=pd.concat([combine, together], axis=0)

#      combine=combine.rename(columns={ 0: 'CI 0.025' })
#      combine=combine.rename(columns={ 1: 'CI 0.975' })

#      combine.set_index(["Cancer_type"], append=True, inplace = True, drop = True)
#      combine.reset_index(inplace = True)
#      combine.drop(columns=['row_number', 'level_0'], axis=1, inplace=True)
#      combine["TERT"]=sys.argv[5]
#      combine["TP53"]=sys.argv[6]
#      combine["SEF"]=sys.argv[7]
#      combine["TERT_Promoter_Analysis"]=sys.argv[8]
#      combine["TP53_G245_Analysis"]=sys.argv[9]
#      combine["TP53_R273_Analysis"]=sys.argv[10]
#      combine["TP53_R248_Analysis"]=sys.argv[11]
#      combine = combine[['Cancer_type', 'Race', 'TERT', 'TP53', 'SEF', 'TERT_Promoter_Analysis', 'TP53_G245_Analysis', 'TP53_R273_Analysis',
#                   'TP53_R248_Analysis', 'covariates', 'coefficient', 'p-values', 'z-score', 'CI 0.025', 'CI 0.975']]
#      combine.to_csv(table_output)
#      print(combine)

together.reset_index(inplace=True)
together.to_csv(risk_output)
