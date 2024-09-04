# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 10:55:53 2022

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
import scipy
import os
import matplotlib
import argparse
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

#matplotlib.use('Agg')
matplotlib.rcParams['font.family'] ='Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import glob
import collections
import seaborn as sns


os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance")
feature_df = pd.read_excel("cluster_heatmap_feature_lymph_node.xlsx",index_col=0)


os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\tmt_cluster_and_feature")
for column in ['Tumor_size']:
    
    
    
    plt.figure(figsize=(2,3))
    column_df = feature_df[[column,"subgroup"]]
    column_df = column_df.dropna(how="any")
    
    palette= {1:"tab:green",2:"tab:blue",3:"tab:red"}
    
 #   sns.boxplot(x="subgroup",y=column,data=feature_df)
    sns.stripplot(x="subgroup",y=column,data=feature_df,palette=palette,alpha=0.3)    
    g1_list = column_df.loc[column_df["subgroup"]==1,column].values
    g2_list = column_df.loc[column_df["subgroup"]==2,column].values
    g3_list = column_df.loc[column_df["subgroup"]==3,column].values
        
    list_ = [g1_list,g2_list,g3_list]
    
    mean_list = [np.mean(i) for i in list_]
    
    std_list = [np.std(i) for i in list_]
    
    plt.bar([0,1,2], mean_list,
       yerr=std_list,
       align='center',
       alpha=0.5,
       ecolor='black',
       capsize=10,
       color = ["tab:green","tab:blue","tab:red"])    
    
    
    
    pvalue1 = '{:.2e}'.format(scipy.stats.ttest_ind(g1_list,g2_list)[1])
    pvalue2 = '{:.2e}'.format(scipy.stats.ttest_ind(g2_list,g3_list)[1])
    pvalue3 = '{:.2e}'.format(scipy.stats.ttest_ind(g1_list,g3_list)[1])
    
    max_value = max(column_df[column])
 #   min_value = min(column_df[column])
    
    interval = int(max_value)/6
  #  interval = 1
    
    size=9  
    plt.plot([0,1],[max_value+interval,max_value+interval],c="black")
    plt.text(0.1,max_value+interval,pvalue1,size=size,verticalalignment='bottom')
    
    plt.plot([1,2],[max_value+2*interval,max_value+2*interval],c="black")
    plt.text(1.1,max_value+2*interval,pvalue2,size=size,verticalalignment='bottom')
    
    plt.plot([0,2],[max_value+3*interval,max_value+3*interval],c="black")
    plt.text(0.5,max_value+3*interval,pvalue3,size=size,verticalalignment='bottom')
 #   ax.set_ylim([-3.5,4])
   # ax.xaxis.set_minor_locator(AutoMinorLocator())
  #  ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xticks([0,1,2],["G1(%s)"%(len(g1_list)),"G2(%s)"%(len(g2_list)),"G3(%s)"%(len(g3_list))],size=size)
    #plt.title("TMT")    
    plt.savefig("subgroup_%s_barplot.pdf"%(column),dpi=300,bbox_inches="tight")
    plt.savefig("subgroup_%s_barplot.png"%(column),dpi=300,bbox_inches="tight")
    
    
    
    
    
for column in [ 'Tumor_size','ImmuneScore',"StromalScore_estimate"]:
    
    
 #   column = 'purity_absolute'
    plt.figure(figsize=(2,3))
    column_df = feature_df[[column,"subgroup"]]
    column_df = column_df.dropna(how="any")
    
    palette= {1:"tab:green",2:"tab:blue",3:"tab:red"}
    
    sns.boxplot(x="subgroup",y=column,data=feature_df,color="white")
    sns.stripplot(x="subgroup",y=column,data=feature_df,palette=palette,alpha=0.3)    
    g1_list = column_df.loc[column_df["subgroup"]==1,column].values
    g2_list = column_df.loc[column_df["subgroup"]==2,column].values
    g3_list = column_df.loc[column_df["subgroup"]==3,column].values
        
#    list_ = [g1_list,g2_list,g3_list]
#    
#    mean_list = [np.mean(i) for i in list_]
#    
#    std_list = [np.std(i) for i in list_]
    
#    plt.bar([0,1,2], mean_list,
#       yerr=std_list,
#       align='center',
#       alpha=0.5,
#       ecolor='black',
#       capsize=10,
#       color = ["tab:green","tab:blue","tab:red"])    
#    
    
    
    pvalue1 = '{:.2e}'.format(scipy.stats.ttest_ind(g1_list,g2_list)[1])
    pvalue2 = '{:.2e}'.format(scipy.stats.ttest_ind(g2_list,g3_list)[1])
    pvalue3 = '{:.2e}'.format(scipy.stats.ttest_ind(g1_list,g3_list)[1])
    
    max_value = max(column_df[column])
    min_value = min(column_df[column])
#    
    interval = int(max_value-min_value)/6
    
    interval = (max_value-min_value)/6
    

    size=9  
    plt.plot([0,1],[max_value+interval,max_value+interval],c="black")
    plt.text(0.1,max_value+interval,pvalue1,size=size,verticalalignment='bottom')
    
    plt.plot([1,2],[max_value+2*interval,max_value+2*interval],c="black")
    plt.text(1.1,max_value+2*interval,pvalue2,size=size,verticalalignment='bottom')
    
    plt.plot([0,2],[max_value+3*interval,max_value+3*interval],c="black")
    plt.text(0.5,max_value+3*interval,pvalue3,size=size,verticalalignment='bottom')
 #   ax.set_ylim([-3.5,4])
   # ax.xaxis.set_minor_locator(AutoMinorLocator())
  #  ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xticks([0,1,2],["G1(%s)"%(len(g1_list)),"G2(%s)"%(len(g2_list)),"G3(%s)"%(len(g3_list))],size=size)
    #plt.title("TMT")    
#    plt.savefig("subgroup_%s_boxplot.pdf"%(column),dpi=300,bbox_inches="tight")
#    plt.savefig("subgroup_%s_boxplot.png"%(column),dpi=300,bbox_inches="tight")    
    



for column in ['SCC(ng/ml)', 'HGB(g/L)',
       'Lymph Node Status(0.Negative; 1.Positive)',
       'Vessel carcinoma embolus(0.Negative; 1.Positive)', 'PIK3CA', 'TP53', 'KRAS', '11q25    Del',
       '3q29     Amp', '2q37.1   Del', '7q11.21  Amp']:
    
    column_df = feature_df[[column,"subgroup"]]
    column_df = column_df.loc[column_df[column]!="white"]
    
    yes_ratio = []
    
    for i in range(1,4):
        g_df = column_df.loc[column_df['subgroup']==i]
        
        n = len(g_df)
        
        yes_df = g_df.loc[g_df[column]=="black"]
        
        m = len(yes_df)
        
        ratio = m/n
        
        yes_ratio.append(ratio)
    
    
    plt.figure(figsize=(2,3))
    

    plt.bar(["G1","G2","G3"],yes_ratio,color=["tab:green","tab:blue","tab:red"])
        
    plt.ylabel("Positive ratio")
    plt.title(column)
    column = column.replace("/","_")
    
    
#    plt.ylim([0,1])
    plt.savefig("subgroup_%s_bar.pdf"%(column),dpi=300,bbox_inches="tight")   
    plt.savefig("subgroup_%s_bar.png"%(column),dpi=300,bbox_inches="tight")         
    
















    