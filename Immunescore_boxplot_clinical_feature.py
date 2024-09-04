# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 20:27:43 2022

@author: guixiuqi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import scipy.stats as stats
import seaborn as sns
from scipy.stats import spearmanr
from scipy.stats import pearsonr
#from matplotlib_venn import venn3
#import scipy.stats as stats
#from scipy.stats import chi2_contingency
#import scipy
import matplotlib
import argparse
from scipy.stats import f_oneway
os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\cluster_pca_reverse_zmap")
import reverse_zmap_samplecluster
import reverse_zmap_pca
import os
#matplotlib.use('Agg')
matplotlib.rcParams['font.family'] ='Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from scipy.stats import ranksums
sample_info = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",
                            index_col=0)

sample_info = sample_info.loc[sample_info["sample type"]!="Normal"]

import seaborn as sns



os.chdir(r"G:\cervical_cancer_multiomics\results\expression_table")

ImmuneScore_dict = {}

for index in sample_info.index:
    purity = sample_info.loc[index,"ImmuneScore"]
    sample_type = sample_info.loc[index,"Histological_Diagnosis"]
    
    if sample_type in ImmuneScore_dict.keys():
        if pd.isna(purity):
            pass
        else:
            ImmuneScore_dict[sample_type].append(purity)
    else:
        if pd.isna(purity):
            pass
        else:    
            ImmuneScore_dict[sample_type] = [purity]
    

results = ranksums(ImmuneScore_dict['Cervical Squamous Cell Carcinoma'],
         ImmuneScore_dict['Cervical Adenocarcinoma'])
    



plt.figure(figsize=(3,5))
g=sns.boxplot(x='Histological_Diagnosis',y='ImmuneScore',data=sample_info,color="white")
sns.stripplot(x='Histological_Diagnosis',y='ImmuneScore',data=sample_info)
plt.xticks(g.get_xticks(),["Squamous","Small Cell","Adenocarcinoma","Adenosqamous","Other"],rotation=90)

plt.xlabel("Histological Diagnosis",size=15)
plt.ylabel("ImmuneScore",size=15)
plt.title('Squamous vs Adenocarcinoma\nWilcoxon rank-sum test pvalue=%s'%('%.2e'%(results[1])))

plt.savefig("ImmuneScore_across_cancer_type.pdf",dpi=300,bbox_inches="tight")

plt.savefig("ImmuneScore_across_cancer_type.png",dpi=300,bbox_inches="tight")




grade_ImmuneScore_dict = {}

for index in sample_info.index:
    purity = sample_info.loc[index,"ImmuneScore"]
    sample_type = sample_info.loc[index,"Grade"]
    
    if sample_type in grade_ImmuneScore_dict.keys():
        if pd.isna(purity):
            pass
        else:
            grade_ImmuneScore_dict[sample_type].append(purity)
    else:
        if pd.isna(purity):
            pass
        else:    
            grade_ImmuneScore_dict[sample_type] = [purity]
    

results = ranksums(grade_ImmuneScore_dict['G2'],
         grade_ImmuneScore_dict['G3'])






plt.figure(figsize=(3,5))
g=sns.boxplot(x='Grade',y='ImmuneScore',data=sample_info,color="white",order=["G1","G2","G3"])
sns.stripplot(x='Grade',y='ImmuneScore',data=sample_info,order=["G1","G2","G3"])

plt.xlabel("Grade",size=15)
plt.ylabel("ImmuneScore",size=15)
plt.title('G2 vs G3\nWilcoxon rank-sum test pvalue=%s'%('%.2e'%(results[1])))

plt.savefig("ImmuneScore_across_grade.pdf",dpi=300,bbox_inches="tight")



plt.savefig("ImmuneScore_across_grade.png",dpi=300,bbox_inches="tight")











stage_ImmuneScore_dict = {}

for index in sample_info.index:
    purity = sample_info.loc[index,"ImmuneScore"]
    sample_type = sample_info.loc[index,"Final stage"]
    
    if sample_type in stage_ImmuneScore_dict.keys():
        if pd.isna(purity):
            pass
        else:
            stage_ImmuneScore_dict[sample_type].append(purity)
    else:
        if pd.isna(purity):
            pass
        else:    
            stage_ImmuneScore_dict[sample_type] = [purity]
    

results = f_oneway(stage_ImmuneScore_dict['I'],stage_ImmuneScore_dict['II'],
         stage_ImmuneScore_dict['III'])


results = ranksums(stage_ImmuneScore_dict['I'],
         stage_ImmuneScore_dict['III'])



plt.figure(figsize=(3,5))
g=sns.boxplot(x='Final stage',y='ImmuneScore',data=sample_info,color="white",order=["I","II","III","IV"])
sns.stripplot(x='Final stage',y='ImmuneScore',data=sample_info,order=["I","II","III","IV"])

plt.xlabel("Figo Stage",size=15)
plt.ylabel("ImmuneScore",size=15)
plt.title('I vs III\nWilcoxon rank-sum test pvalue=%s'%('%.2e'%(results[1])))

plt.savefig("ImmuneScore_across_stage.pdf",dpi=300,bbox_inches="tight")


plt.savefig("ImmuneScore_across_stage.png",dpi=300,bbox_inches="tight")







stage_ImmuneScore_dict = {}

for index in sample_info.index:
    purity = sample_info.loc[index,"ImmuneScore"]
    sample_type = sample_info.loc[index,"Final clade"]
    
    if sample_type in stage_ImmuneScore_dict.keys():
        if pd.isna(purity):
            pass
        else:
            stage_ImmuneScore_dict[sample_type].append(purity)
    else:
        if pd.isna(purity):
            pass
        else:    
            stage_ImmuneScore_dict[sample_type] = [purity]
    
#
#results = f_oneway(stage_ImmuneScore_dict['A7'],stage_ImmuneScore_dict['A9'],
#         stage_ImmuneScore_dict['III'])


results = ranksums(stage_ImmuneScore_dict['A7'],
         stage_ImmuneScore_dict['A9'])

hpv_df = sample_info[["ImmuneScore","Final clade"]]

hpv_df = hpv_df.loc[(hpv_df['Final clade']=="A7")|(hpv_df["Final clade"]=="A9")]

hpv_df = hpv_df.dropna(how="any")

hpv_df.to_excel("immunescore_hpv_clade.xlsx")


plt.figure(figsize=(3,5))
g=sns.boxplot(x='Final clade',y='ImmuneScore',data=sample_info,color="white",order=["A7","A9"])
sns.stripplot(x='Final clade',y='ImmuneScore',data=sample_info,order=["A7","A9"])

plt.xlabel("HPV clade",size=15)
plt.ylabel("ImmuneScore",size=15)
plt.title('A7 vs A9\nWilcoxon rank-sum test pvalue=%s'%('%.2e'%(results[1])))

plt.savefig("ImmuneScore_across_HPV.pdf",dpi=300,bbox_inches="tight")


plt.savefig("ImmuneScore_across_HPV.png",dpi=300,bbox_inches="tight")





