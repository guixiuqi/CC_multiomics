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
plt.rcParams['xtick.labelsize']=6
plt.rcParams['ytick.labelsize']=6
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

import collections




import scipy

def expression_boxplot(tmt_df,gene_name="TSC2",type_=""):
    
    group_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\cluster_3_add_normal.csv",sep="\t",index_col=0)
    os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\feature_PFS_cox_new_survival_data\boxplot_expression_across_groups\\")
    common_sample = list(set(tmt_df.columns)&set(group_df.index))    
    sub_tmt_df = tmt_df[common_sample].T
    sub_group_df = group_df.loc[common_sample]    
    all_concat_df = pd.concat([sub_tmt_df,sub_group_df],axis=1)
    gene_df = all_concat_df[[gene_name,'x']]
    gene_df = gene_df.dropna(how="any")
    
    a = collections.Counter(gene_df['x'])
    print(a)
    
#    gene_df.to_excel("%s_%s.xlsx"%(gene_name,type_))
    
    
    color_dict = {"Normal":"grey","1":"tab:green","2":"tab:blue","3":"tab:red"}
    hue_order=["Normal","1","2","3"]
    
    fig = plt.figure(figsize=(2,5),dpi=300)   
    plt.rcParams['xtick.labelsize']=6
    plt.rcParams['ytick.labelsize']=6
    ax=fig.gca()
    
    plotting_parameters = {
    'data':    gene_df,
    'x':       'x',
    'y':       gene_name,
    'palette': color_dict,
#    "hue_order": hue_order
    'order': hue_order,
    "showfliers":False
}
    
    
    plotting_parameters1 = {
    'data':    gene_df,
    'x':       'x',
    'y':       gene_name,
    'color': "black",
#    "hue_order": hue_order
    'order': hue_order,
    'size':2
}    
    
    pairs = []
    
    for i in range(0,3):
        for j in range(i,4):
            pair = (hue_order[i],hue_order[j])
            pairs.append(pair)

  
        # Plot with seaborn
        
    ax = sns.boxplot(fliersize=1,**plotting_parameters)
    sns.stripplot(**plotting_parameters1)
    # Add annotations
    annotator = Annotator(ax, pairs, **plotting_parameters)
    annotator.configure(test="t-test_welch",text_format='full',verbose=False).apply_and_annotate()

    ax.set_xticklabels(ax.get_xticklabels(),rotation=30)    
    
    plt.title(type_)
    #ax.yaxis.set_minor_locator(AutoMinorLocator())


    tsc2_dir=r"G:\cervical_cancer_multiomics\results\TMT_protein\feature_PFS_cox_new_survival_data\boxplot_expression_across_groups\\"
#    plt.savefig(tsc2_dir+gene_name+"_%s_z_statistic_boxplot_"%type_+".pdf",dpi=300,bbox_inches="tight")
#    plt.close("all")
    plt.show()




     #   prepare_gsva_r(gene,fianl_z_df)

tmt_df = pd.read_csv(r"tmt_z_statistic_table.csv",index_col=0,sep="\t")

#gene_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\TMT_protein\feature_PFS_cox\tmt_dia_chemo_radio_response_associated_protein.xlsx",index_col=0)

DIA_df = pd.read_csv(r"dia_z_statistic_table.xls",index_col=0,sep="\t")
tmt_phos_df = pd.read_csv(r"tmt_phosprotein_normalized_log2_statistic_table_remove_error_sample.csv",index_col=0,sep=",")


for gene in ['PRKCB']:
    expression_boxplot(tmt_df,gene_name=gene,type_="TMT")
    expression_boxplot(DIA_df,gene_name=gene,type_="DIA")
    expression_boxplot(tmt_phos_df,gene_name=gene,type_="TMT_phos")
    
    













    
        
   
    
