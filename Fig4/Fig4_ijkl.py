# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 16:56:29 2023

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
import os
#matplotlib.use('Agg')
matplotlib.rcParams['font.family'] ='Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import seaborn as sns
import collections
import scipy
import gseapy
import matplotlib as mpl
import matplotlib.cm as cm
from sklearn import preprocessing
from scipy import stats


from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test
kmf = KaplanMeierFitter()


y2017_sample_sub_clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\TCGA_DATA\2014_nature_sample_type_hpv.xlsx",
                                           index_col=0)




tcga_df = pd.read_csv(r"G:\cervical_cancer_multiomics\TCGA_DATA\TCGA-CESC.htseq_fpkm-uq_rename.tsv",
                      sep="\t",index_col=0)
survival_df = pd.read_csv(r"F:\project\Cervical_cancer\TCGA_DATA\TCGA-CESC.survival.tsv",sep="\t",index_col=0)



overlap_sample= list(set(tcga_df.columns)&set(survival_df.index))



sample_type_list = []
for index in overlap_sample:
    sample_type_list.append(index.split("-")[-1])




final_sample = []

for sample in overlap_sample:
    if "11A" in sample or "06A" in sample:
        pass
    else:
        final_sample.append(sample)

#survival_df = survival_df.loc[final_sample]


tcga_df = tcga_df[final_sample]
survival_df = survival_df.loc[final_sample]



for index in survival_df.index:
    survival_df.loc[index,'OS.time'] = survival_df.loc[index,'OS.time']/30








from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test




def gene_survival(survival_df,gene="PRKCB"):
    

    TCGA_gene_df = tcga_df.loc[gene,survival_df.index].to_frame()
    TCGA_gene_df["OS.time"] = survival_df.loc[TCGA_gene_df.index,"OS.time"]
    TCGA_gene_df['OS'] = survival_df.loc[TCGA_gene_df.index,"OS"]
    TCGA_gene_df = TCGA_gene_df.sort_values(by=gene)
    
    min_ = 0
    
    min_p = 1
    
    for i in range(10,280):

    
        TCGA_gene_df['group']=["Low"]*i +['High']*(len(TCGA_gene_df)-i)

        T1 =  TCGA_gene_df.loc[TCGA_gene_df["group"]=="High"]["OS.time"].values
        E1 = TCGA_gene_df.loc[TCGA_gene_df["group"]=="High"]["OS"].values
        T2 = TCGA_gene_df.loc[TCGA_gene_df["group"]=="Low"]["OS.time"].values
        E2 = TCGA_gene_df.loc[TCGA_gene_df["group"]=="Low"]["OS"].values
        results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
        
        pvalue = results.p_value
        
        if pvalue < min_p:
            min_=i
            min_p = pvalue
            

    i=min_
    TCGA_gene_df['group']=["Low"]*i +['High']*(len(TCGA_gene_df)-i)
    TCGA_gene_df = TCGA_gene_df
    
    TCGA_gene_df.to_excel(r"G:\cervical_cancer_multiomics\results\TMT_protein\feature_PFS_cox_new_survival_data\TCGA_fig4\%s_TCGA_survival_291.xlsx"%gene)
    
    color_dict={"High":"red","Low":"blue"}

    plt.figure(figsize=(7,3))
    ax = plt.subplot(121)
    sns.violinplot(x='group',y=gene,data=TCGA_gene_df)
    
    
    ax = plt.subplot(122)
    kmf1 = KaplanMeierFitter()
    for name, grouped_df in TCGA_gene_df.groupby("group"):
        kmf1.fit(grouped_df["OS.time"], grouped_df["OS"], label=name+"(%s)"%(len(grouped_df)),
                alpha =0.1)

        kmf1.plot(ax=ax,show_censors=True,ci_show=False,color=color_dict[name],
              censor_styles={'ms': 6},linewidth=0.5)
    plt.rcParams['xtick.labelsize']=12
    plt.rcParams['ytick.labelsize']=12
    plt.xlabel("OS(months)",size=12)
    plt.ylabel("Survival probability",size=12)

    plt.ylim([0,1])
    
    from lifelines.statistics import multivariate_logrank_test
    T1 =  TCGA_gene_df.loc[TCGA_gene_df["group"]=="High"]["OS.time"].values
    E1 = TCGA_gene_df.loc[TCGA_gene_df["group"]=="High"]["OS"].values
    T2 = TCGA_gene_df.loc[TCGA_gene_df["group"]=="Low"]["OS.time"].values
    E2 = TCGA_gene_df.loc[TCGA_gene_df["group"]=="Low"]["OS"].values
    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    
    pvalue = results.p_value
           
    ax.text(10,0.6,"pvalue:%s"%(str("{:.2e}".format(pvalue))),size=12)
    plt.tight_layout()
    plt.savefig(r"\%s_tcga_survival_291.pdf"%(gene),dpi=300,bbox_inches="tight")

gene_list = pd.read_excel(r"fig4_gene.xlsx",index_col=0)


for gene in gene_list.index:

    gene_survival(survival_df,gene=gene)










