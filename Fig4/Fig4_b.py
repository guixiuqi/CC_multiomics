# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 16:01:58 2022

@author: guixiuqi
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from matplotlib_venn import venn3
import os
os.chdir(r"D:\python_module")
import venn
import mkdir_
import get_rgb_color_list

import scipy.stats as stats
from scipy.stats import chi2_contingency
import scipy
import matplotlib
import argparse

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
from scipy.stats import chi2_contingency
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import ranksums
from statannotations.Annotator import Annotator
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import matplotlib.cm as cm
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'






   

def ep300_FOSL2_gene_level_correlation():
    
    imputation_f = r"G:\cervical_cancer_multiomics\results\expression_table\tmt_protein_level_ace_half_quantified_dreamai_imputation_res.csv"
    
    tmt_f = pd.read_csv(imputation_f,sep=",",index_col=0)
    
    normal_sample = []
    
    for column in tmt_f.columns:
        if column.startswith("N"):
            normal_sample.append(column)
        else:
            pass
    
    tmt_f_t = tmt_f.T
    
    # dir_ = r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\%s"%(type_)
    # mkdir_.mkdir(dir_)
    # os.chdir(dir_)
    
    tmt_cluster = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\cluster_3.csv",
                              sep=",",index_col=0)
    tmt_cluster = tmt_cluster.sort_values(by="x")
    
    tmt_f_t = tmt_f_t.loc[tmt_cluster.index]
    
    tmt_f_t["tmt_group"] = tmt_cluster.loc[tmt_f_t.index,"x"]
    
    sub_df = tmt_f_t[["EP300","FOSL2"]]
    sub_df.to_excel(r"G:\cervical_cancer_multiomics\results\FOSL2_k222\all_tumor_sample\EP300_FOSL2_acetylation.xlsx")
    
    EP300_list = []
    
    for column in tmt_f_t.columns:
        if column.startswith("EP300"):
            EP300_list.append(column)
            
    for i in EP300_list:
        plt.figure(figsize=(3,3))
        plt.rcParams['xtick.labelsize']=10
        plt.rcParams['ytick.labelsize']=10
        
        plt.scatter(tmt_f_t["FOSL2"],tmt_f_t[i],s=8)
        
        plt.xlabel("FOSL2 acetylation level",size=10)
        
        plt.ylabel(i+" acetylation level",size=10)
        
        slope, intercept, r, p, std_err = stats.linregress(tmt_f_t["FOSL2"],tmt_f_t[i])
        def myfunc(x):
            return slope * x + intercept
        
        mymodel = list(map(myfunc,tmt_f_t["FOSL2"]))
      #  plt.scatter(arm_df.loc["3q"], arm_df.loc["11q25"])
        plt.plot(tmt_f_t["FOSL2"], mymodel,color="black")
        stats_,pvalue = stats.spearmanr(tmt_f_t["FOSL2"],tmt_f_t[i])
        plt.title("Spearmanr correlation=%s\npvalue=%s"%('%.2f'%(stats_),"%.1e"%(pvalue)))         
       # plt.title("stats:%s\npvalue:%s"%("%.2f"%stats_,"%.2e"%pvalue_))
        plt.savefig(r"%s_%s_spearmanr.pdf"%("FOSL2",i),dpi=300,bbox_inches="tight")
        plt.savefig(r"%s_%s_spearmanr.png"%("FOSL2",i),dpi=300,bbox_inches="tight")
  



ep300_FOSL2_gene_level_correlation()



    
        
   
    
