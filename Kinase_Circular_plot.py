# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:04:37 2023

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
os.chdir(r"E:\python_module")
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
matplotlib.rcParams['axes.linewidth'] = 0.5
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



group_specific_dict = {'group1':["PRKAA1","PRKAA2","PRKACA","PRKACB"],
                       'group2':["CDK1","CDK2","EIF2AK2","MAPKAPK5","TTK"],
                       'group3':["PRKCB","PRKCD"]}


kinase_score = r"G:\cervical_cancer_multiomics\results\KSEA\pairwise_comparison\group1\KSEA Kinase Scores.csv"

kinase_substrate = r"G:\cervical_cancer_multiomics\results\KSEA\pairwise_comparison\group1\Kinase-Substrate Links.csv"
    

def bar_color():
    group_dict = {"group1":"tab:green","group2":"tab:blue","group3":"tab:red"}
    color_list = []
    
    for key in group_specific_dict.keys():
        l_ = len(group_specific_dict[key])
        color_list += [group_dict[key]]*l_
    return color_list

color_list = bar_color()
        

def group_score():
    
    ksea_score = []
    
    gene_n_list = []
    for key in group_specific_dict.keys():
        
        kinase_score_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\KSEA\pairwise_comparison\%s\KSEA Kinase Scores.csv"%(key),sep=",",index_col=0)
        
        gene_list = group_specific_dict[key]
        gene_n_list += gene_list
        
        for gene in gene_list:
            score_ = kinase_score_df.loc[gene,'z.score']
            ksea_score.append(score_)
    return ksea_score,gene_n_list

ksea_score_list,gene_n_list = group_score()

def phos_position():
    
    sub_strate_list = []
    
    t_statistic_dict = {}
    
    for key in group_specific_dict.keys():
    
        kinase_substrate_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\KSEA\pairwise_comparison\%s\Kinase-Substrate Links.csv"%(key),sep=",",index_col=0)
        
        gene__list = group_specific_dict[key]
        
        for gene in gene__list:
            
            gene_pos_df = kinase_substrate_df.loc[gene]
            
            gene_pos_df_sort = gene_pos_df.sort_values(by="log2FC",ascending=False)
            
            sub_strate=""
            for i in range(0,8):
                if i+1 < len(gene_pos_df_sort):
                    sub_ = gene_pos_df_sort.iloc[i,0]
                    pos_ = gene_pos_df_sort.iloc[i,1]
                    sub_strate = sub_strate+sub_+":"+pos_+"\n"
                    t_ = gene_pos_df_sort.iloc[i,3]
                    
                    t_statistic_dict[gene+"_"+sub_+":"+pos_] = t_
                else:
                    continue
            sub_strate = sub_strate.strip("\n")            
            sub_strate_list.append(sub_strate)
            
            
    return sub_strate_list,t_statistic_dict


sub_strate_list,t_statistic_dict = phos_position()
                    

polar_df = pd.DataFrame.from_dict(t_statistic_dict, orient='index')

polar_df.to_excel(r"G:\cervical_cancer_multiomics\results\KSEA\pairwise_comparison\group_specific_kinase_polar_data.xlsx")

gene_n = len(color_list)

fig = plt.figure(figsize=(3,3),dpi=300)
ax1 = plt.subplot(111, projection='polar')

theta = np.arange(0,2*np.pi,2*np.pi/gene_n)

width = 2*np.pi/(gene_n)

plt.bar(theta,[3.5]*gene_n,width=width,color = color_list,edgecolor="white")


for i in range(0,gene_n):
    plt.text(theta[i],3.3,gene_n_list[i],rotation = theta[i]/(2*np.pi)*360-90, 
             verticalalignment='center',horizontalalignment='center',size=5)
    
    plt.text(theta[i],4.5,sub_strate_list[i],rotation = theta[i]/(2*np.pi)*360, 
             verticalalignment='center',horizontalalignment='center',size=5)
    
    pos_list = sub_strate_list[i].split("\n")
    
    num = len(pos_list)
    
#    plot_w = width/num
    
    left_ = theta[i]-0.23
    
    right  = theta[i]+0.23
    
    w_ = (right-left_)/(num-1)
    
    x_list_ = [left_] + [left_+m*w_ for m in range(1,num-1)] +[right]
#    plt.scatter(x_list_,[np.log2(t_statistic_dict[gene_n_list[i]+pos]) for pos in pos_list],c="black",s=5)
    
    plt.plot(x_list_[::-1],[np.log2(t_statistic_dict[gene_n_list[i]+"_"+pos]) for pos in pos_list],c="black",linestyle="--",linewidth=0.5,
                      marker='.',markersize=1)
    
    
ax1.yaxis.grid(True)

ax1.xaxis.set_visible(False) 

plt.yticks([1,2,3],[2,4,8])

plt.savefig(r"G:\cervical_cancer_multiomics\results\KSEA\pairwise_comparison\group_specific_kinase_polar.pdf",dpi=300,bbox_inches="tight")




