# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 19:55:26 2023

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


    
tmt_protein_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_pro_half_quantified_dreamai_imputation_res.csv",
                             sep=",",index_col=0)

dia_protein_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\dia_pro_half_quantified_dreamai_imputation_res.csv",
                             sep=",",index_col=0)


tmt_phosite_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_protein_level_phos_half_quantified_dreamai_imputation_res.csv",
                             sep=",",index_col=0)

tmt_ace_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_protein_level_ace_half_quantified_dreamai_imputation_res.csv",
                         sep=",",index_col=0)

rna_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\fpkm_remove_error_sample.xls",
                     sep="\t",index_col=0)


groups_df = pd.read_csv("G:\\cervical_cancer_multiomics\\results\\TMT_protein\\top_variance\\euclidean_km\\cluster_3.csv",sep=',',index_col=0)
groups_df_sort = groups_df.sort_values(by="x")



group_specific_dict = {'group1':["PRKAA1","PRKAA2","PRKACA","PRKACB"],
                       'group2':["CDK1","CDK2","EIF2AK2","MAPKAPK5","TTK"],
                       'group3':["PRKCB","PRKCD"]}



group_specific_gene_list = []
for key in group_specific_dict.keys():    
    list__ = group_specific_dict[key]
    group_specific_gene_list += list__
    

    
def heatmap_(tmt_protein_df):
    
    tumor_sample = []
    
    for sample in tmt_protein_df.columns:
        
        if sample.startswith("T"):
            tumor_sample.append(sample)
        else:
            pass
        
    final_sample = list(set(tumor_sample)&set(groups_df_sort.index))
    
    tumor_group_df = groups_df.loc[final_sample]
    
    tumor_tmt_df = tmt_protein_df[final_sample]
    
    tumor_tmt_df = tumor_tmt_df.sub(tumor_tmt_df.median(axis=1), axis=0)
    
    mean_list2 = []
    for gene in group_specific_gene_list:
        
        if gene in tumor_tmt_df.index:
            mean_list = []
            
            for i in [1,2,3]:
                
                group_sample = list(tumor_group_df.loc[tumor_group_df['x']==i].index)
                
                group_mean = np.mean(tumor_tmt_df.loc[gene,group_sample])
                
                mean_list.append(group_mean)
                
            mean_list2.append(mean_list)
                
        else:
            
            mean_list2.append([np.nan,np.nan,np.nan])
            
    mean_df = pd.DataFrame(mean_list2,index=group_specific_gene_list,columns=["Group1","Group2","Group3"])
    
    return mean_df




tmt_mean_df = heatmap_(tmt_protein_df)

dia_mean_df = heatmap_(dia_protein_df)

tmt_phosite_df = heatmap_(tmt_phosite_df)

os.chdir(r"G:\cervical_cancer_multiomics\results\KSEA")

tmt_mean_df.to_excel("tmt_mean_df.xlsx")

dia_mean_df.to_excel("dia_mean_df.xlsx")

tmt_phosite_df.to_excel("tmt_phos_mean_df.xlsx")



def add_legend(ax):
    remove_spine(ax)
    ax.bar([0],[1],width=1,color="white",edgecolor="grey",linewidth=0.3)
    ax.set_xticks([],[])
    ax.set_yticks([],[])
    ax.set_title('Not detected',x=1,y=-1.5,size=10)



import matplotlib.gridspec as gridspec
fig=plt.figure(figsize=(8,4))
gs=gridspec.GridSpec(nrows=20,ncols=50)

ax0 = fig.add_subplot(gs[0:20,0:6])
#ax00 = fig.add_subplot(gs[5:9,9:10])

ax1 = fig.add_subplot(gs[0:20,10:16])
#ax11 = fig.add_subplot(gs[5:9,19:20])

ax2 = fig.add_subplot(gs[0:20,20:26])
ax22 = fig.add_subplot(gs[5:9,29:30])


color="grey"
width= 0.2
map_color ="coolwarm"


vmin_value = -1
vmax_value = 1

vmin_value = -0.5
vmax_value = 0.5

nan_value_color = "white"


mask = tmt_mean_df.isnull()

g= sns.heatmap(tmt_mean_df,cmap=map_color,vmin=vmin_value,ax=ax0,cbar=False,vmax=vmax_value,
               xticklabels=True,linecolor=color,linewidth=width,mask = mask,square=True)

ax0.set_title("Protein(TMT)",size=10)
g.set_facecolor(nan_value_color)



mask = dia_mean_df.isnull()

g1= sns.heatmap(dia_mean_df,cmap=map_color,vmin=vmin_value,ax=ax1,cbar=False,vmax=vmax_value,
               xticklabels=True,linecolor=color,linewidth=width,mask = mask,square=True)

ax1.set_title("Protein(DIA)",size=10)
g.set_facecolor(nan_value_color)



mask = tmt_phosite_df.isnull()

g2 = sns.heatmap(tmt_phosite_df,cmap=map_color,vmin=vmin_value,ax=ax2,cbar_ax=ax22,vmax=vmax_value,
               xticklabels=True,linecolor=color,linewidth=width,mask = mask,square=True,cbar_kws= {"shrink": 0.6,"label":"Centered z-statistic"})

ax2.set_title("Phospho protein",size=10)
g2.set_facecolor(nan_value_color)

plt.savefig(r"G:\cervical_cancer_multiomics\results\KSEA\group_specific_kinase_heatmap_new.pdf",dpi=300)








