# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:01:42 2022

@author: guixi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] ='Arial'
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
import os
import math
from scipy import stats
import seaborn as sns
import scipy
import collections
import numpy as np
import matplotlib.pyplot as plt
#sample_sub_clinical_df = pd.read_excel(r"C:\Users\guixi\Desktop\multi_omics\results\sample_clinical_info_2022_5_13.xlsx",index_col=0)
#sample_sub_clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_add_tumor_purity_estimation_2022_6_30_check_hpv_add_clade.xlsx",index_col=0) 
sample_sub_clinical_df  = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_add_tumor_purity_estimation_check_hpv_add_clade_follow_up_2022_8_31.csv",sep=",",index_col=0)

yujing_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\sample_info_pie\2022_3_20_sample_info_from_yujing.xlsx",index_col=0)
sample_sub_clinical_df["Grade"] = yujing_df.loc[sample_sub_clinical_df.index,"Grade"]

sample_sub_clinical_df.to_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_add_tumor_purity_estimation_check_hpv_add_clade_follow_up_2022_10_10.xlsx")

import os

os.chdir(r"G:\cervical_cancer_multiomics\results\sample_info_pie")
## 183
tumor_sample_list = []  #139
normal_sample_list = []  #40
for sample in sample_sub_clinical_df.index:
    if sample.startswith("T"):
        tumor_sample_list.append(sample)
    else:
        normal_sample_list.append(sample)
        
tumor_sample_df = sample_sub_clinical_df.loc[tumor_sample_list]



 
data_figs1c = tumor_sample_df[["Grade","Final stage","Treatment"]]   

data_figs1c.to_excel("figs1c.xlsx")

#tumor_sample_df["follow_up_state(2022.08)"] = follow_up_df.loc[tumor_sample_df.index,"末次随访状态（截止至2022年8月）"]




normal_sample_df = sample_sub_clinical_df.loc[normal_sample_list]
surgery_df = tumor_sample_df.loc[tumor_sample_df['Treatment']=="Surgery"]
chemo_df = tumor_sample_df.loc[tumor_sample_df['Treatment']!="Surgery"]

progress_df = tumor_sample_df.loc[tumor_sample_df['follow_up_state(2022.08)']!="No-progress"]


no_progress_df = tumor_sample_df.loc[tumor_sample_df['follow_up_state(2022.08)']=="No-progress"]

def bar_stage(tumor_sample_df,column,title="",sample_pop="tumor_sample"):


    
    stage_df = tumor_sample_df[column].to_frame()
    stage_df = stage_df.dropna(how="any")
    stage = list(stage_df[column].values)
    
    new_stage_list = []
    
    for s in stage:
     #   new_stage=s
        new_stage = s.replace("Ⅰ","I").replace("Ⅱ","II").replace("Ⅲ","III").replace("Ⅳ","IV")
        new_stage_list.append(new_stage)
    
    
    stage_dict = collections.Counter(new_stage_list)
    
    sorted_keys = sorted(stage_dict)
    
    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    data = [stage_dict[i] for i in sorted_keys]
    ingredients = sorted_keys
    def func(pct, allvals):
        absolute = int(round(pct*np.sum(allvals)/100))
        return "{:.1f}%\n({:d})".format(pct, absolute)
    
    wedges, texts, autotexts = ax.pie(data, autopct=lambda pct: func(pct, data),
                                  textprops=dict(color="w"))
    
    ax.legend(wedges, ingredients,
          title=title,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts, size=8, weight="bold")
    ax.set_title(sample_pop)
    plt.savefig("%s_%s_pie.pdf"%(sample_pop,title),dpi=300,bbox_inches="tight")
    
bar_stage(tumor_sample_df,"Stage",title="Stage",sample_pop="Tumor sample")
bar_stage(tumor_sample_df,"Final stage",title="Final stage",sample_pop="Tumor sample")
bar_stage(tumor_sample_df,"Treatment",title="Treatment",sample_pop="Tumor sample")    

bar_stage(tumor_sample_df,"Grade",title="Grade",sample_pop="Tumor sample")
bar_stage(tumor_sample_df,"Histological_Diagnosis",title="Histological diagnosis",sample_pop="Tumor sample")     

bar_stage(tumor_sample_df,"follow_up_state(2022.08)",title="follow_up_state(2022.08)",sample_pop="Tumor sample")     


bar_stage(normal_sample_df,"Stage",title="Stage",sample_pop="Normal sample")
bar_stage(normal_sample_df,"Final stage",title="Final stage",sample_pop="Normal sample")
bar_stage(normal_sample_df,"Treatment",title="Treatment",sample_pop="Normal sample")    

bar_stage(normal_sample_df,"Grade",title="Grade",sample_pop="Normal sample")    



    
    
    
bar_stage(surgery_df,"Stage",title="Stage",sample_pop="Surgery sample")
bar_stage(surgery_df,"Final stage",title="Final stage",sample_pop="Surgery sample")
bar_stage(surgery_df,"Treatment",title="Treatment",sample_pop="Surgery sample")    
bar_stage(surgery_df,"Grade",title="Grade",sample_pop="Surgery sample")    
bar_stage(surgery_df,"follow_up_state(2022.08)",title="follow_up_state(2022.08)",sample_pop="Surgery sample")
        
    
bar_stage(chemo_df,"Stage",title="Stage",sample_pop="Chemoradiotherapy sample")
bar_stage(chemo_df,"Final stage",title="Final stage",sample_pop="Chemoradiotherapy sample")
bar_stage(chemo_df,"Treatment",title="Treatment",sample_pop="Chemoradiotherapy sample")    
bar_stage(chemo_df,"Grade",title="Grade",sample_pop="Chemoradiotherapy sample")
bar_stage(chemo_df,"chemo response",title="Chemoradiotherapy response",sample_pop="Chemoradiotherapy sample")     

bar_stage(chemo_df,"follow_up_state(2022.08)",title="follow_up_state(2022.08)",sample_pop="Chemoradiotherapy sample") 

 
chemo_progress_df = chemo_df.loc[chemo_df["follow_up_state(2022.08)"]!="No-progress"]
bar_stage(chemo_progress_df,"Stage",title="Stage",sample_pop="Chemoradiotherapy Progress sample")  
   


bar_stage(progress_df,"Stage",title="Stage",sample_pop="Progress sample")    
bar_stage(progress_df,"Treatment",title="Treatment",sample_pop="Progress sample") 
