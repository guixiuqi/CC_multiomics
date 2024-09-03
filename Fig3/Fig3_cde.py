# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 17:46:32 2023

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

from random import random
from random import randint
from random import seed
from numpy import arange
from numpy import mean
from numpy import std
from numpy import absolute
from sklearn.datasets import make_regression
from sklearn.linear_model import RANSACRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from matplotlib import pyplot

    
    
tmt_protein = pd.read_excel(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_tmt\E6_spearman_pvalue_tmt.xlsx",index_col=0)

tmt_phos = pd.read_excel(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_phos\E6_spearman_pvalue_phos.xlsx",index_col=0)


# tmt_protein.columns = ["spearman's r(protein)","P-value(protein)"]

# tmt_phos.columns = ["spearman's r(phosphoprotein)","P-value(phosphoprotein)"]

# all_df = pd.concat([tmt_protein,tmt_phos],axis=1)

# all_df.to_csv(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_tmt_spearman")






overlap_gene = list(set(tmt_phos.index)&set(tmt_protein.index))


plt.figure(figsize=(5,5))
phos_specific_gene = []

for gene in overlap_gene:
    plt.scatter(tmt_protein.loc[gene,"r_spearman"],tmt_phos.loc[gene,"r_spearman"],c="blue",s=2)
    if tmt_phos.loc[gene,"r_spearman"]> 0.3:
        print (gene)
        if tmt_protein.loc[gene,"r_spearman"]<0.05 :
            print(gene,tmt_protein.loc[gene,"r_spearman"],tmt_phos.loc[gene,"r_spearman"])
        
            phos_specific_gene.append(gene)  #H2AX,WEE1
        
for gene in phos_specific_gene:
    plt.annotate(gene,(tmt_protein.loc[gene,"r_spearman"],tmt_phos.loc[gene,"r_spearman"]),
                 (tmt_protein.loc[gene,"r_spearman"]+0.05,tmt_phos.loc[gene,"r_spearman"]+0.05),arrowprops=dict(facecolor='black',arrowstyle='-'))

plt.xlim([-0.5,0.5])
plt.ylim([-0.5,0.5])   
plt.xlabel("spearman_r_with_E6(protein)",size=6)     
plt.ylabel("spearman_r_with_E6(phosphosite)",size=6)
plt.savefig(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_phos\scatterplot_spearman_with_e6.pdf",dpi=300)






os.chdir(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression")
hpv_gene_expression_tpm_raw = pd.read_csv("hpv_gene_expression_tpm.csv",sep="\t",index_col=0)

hpv_gene_expression_tpm = np.log2(hpv_gene_expression_tpm_raw+1)
hpv_gene_expression_tpm = hpv_gene_expression_tpm.drop(index=["E10"])


tmt_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_pro_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)

tmt_df.index = [index.replace("-","_").replace("/","_") for index in tmt_df.index]


phos_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_protein_level_phos_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)

phos_df.index= [index.replace("-","_").replace("/","_") for index in phos_df.index]

ace_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_protein_level_ace_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)

ace_df.index= [index.replace("-","_").replace("/","_") for index in ace_df.index]


overlap_sample = set(tmt_df.columns)&set(hpv_gene_expression_tpm.columns)



hpv_e6_df = hpv_gene_expression_tpm.loc["E6",overlap_sample].to_frame()

sub_tmt_df = tmt_df.loc[["H2AX","WEE1"],overlap_sample]




for gene in phos_specific_gene:
    
    plt.figure(figsize=(8,3))
    plt.subplot(121)
    plt.scatter(hpv_gene_expression_tpm.loc["E6",overlap_sample],tmt_df.loc[gene,overlap_sample],s=5)
    x = hpv_gene_expression_tpm.loc["E6",overlap_sample]
    y = tmt_df.loc[gene,overlap_sample]
    slope, intercept = np.polyfit(x, y, 1)
    plt.plot(x, slope * x + intercept, color='red')
    
    plt.xlabel("mRNA abundance of E6")
    plt.ylabel("protein abundance of %s"%gene)
    plt.title("spearman r:%s\np:%s"%("%.2f"%tmt_protein.loc[gene,"r_spearman"],"%.2e"%tmt_protein.loc[gene,"pvalue"]))
    
    plt.subplot(122)
    plt.scatter(hpv_gene_expression_tpm.loc["E6",overlap_sample],phos_df.loc[gene,overlap_sample],s=5)
    x=hpv_gene_expression_tpm.loc["E6",overlap_sample]
    y=phos_df.loc[gene,overlap_sample]
    slope, intercept = np.polyfit(x, y, 1)
    plt.plot(x, slope * x + intercept, color='red')
    
    plt.xlabel("mRNA abundance of E6")
    plt.ylabel("phosphoprotein abundance of %s"%gene)
    plt.title("spearman r:%s\np:%s"%("%.2f"%tmt_phos.loc[gene,"r_spearman"],"%.2e"%tmt_phos.loc[gene,"pvalue"]))    
    
    plt.savefig(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_phos\%s_protein_phos_cor_with_e6.pdf"%gene,dpi=300,bbox_inches="tight")



# theilsen regression on a dataset with outliers
from random import random
from random import randint
from random import seed
from numpy import arange
from numpy import mean
from numpy import std
from numpy import absolute
from sklearn.datasets import make_regression
from sklearn.linear_model import TheilSenRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from matplotlib.pyplot import plt

# prepare the dataset

# evaluate a model
def evaluate_model(X, y, model):
	# define model evaluation method
	cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
	# evaluate model
	scores = cross_val_score(model, X, y, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
	# force scores to be positive
	return absolute(scores)

# plot the dataset and the model's line of best fit
def plot_best_fit(X, y, model):
	# fut the model on all data
	model.fit(X, y)
	# plot the dataset
#	pyplot.scatter(X, y)
	# plot the line of best fit
	xaxis = arange(X.min(), X.max(), 0.01)
	yaxis = model.predict(xaxis.reshape((len(xaxis), 1)))
	plt.plot(xaxis, yaxis, color='r')
	# show the plot
	plt.title(type(model).__name__)
	#pyplot.show()

# load dataset

# define the model


for gene in phos_specific_gene:
    
    plt.figure(figsize=(8,3))
    plt.subplot(121)
    plt.scatter(hpv_gene_expression_tpm.loc["E6",overlap_sample],tmt_df.loc[gene,overlap_sample],s=5)
    x = np.array([[i] for i in list(hpv_gene_expression_tpm.loc["E6",overlap_sample].values)])
    y = tmt_df.loc[gene,overlap_sample].values
    
    
   # slope, intercept = np.polyfit(x, y, 1)
   # plt.plot(x, slope * x + intercept, color='red')
    model = TheilSenRegressor()
# evaluate model
    results = evaluate_model(x, y, model)
    print('Mean MAE: %.3f (%.3f)' % (mean(results), std(results)))
# plot the line of best fit
    plot_best_fit(x, y, model)
    
    plt.xlabel("mRNA abundance of E6")
    plt.ylabel("protein abundance of %s"%gene)
    plt.title("spearman r:%s\np:%s"%("%.2f"%tmt_protein.loc[gene,"r_spearman"],"%.2e"%tmt_protein.loc[gene,"pvalue"]))
    
    plt.subplot(122)
    plt.scatter(hpv_gene_expression_tpm.loc["E6",overlap_sample],phos_df.loc[gene,overlap_sample],s=5)
    x= np.array([[i] for i in hpv_gene_expression_tpm.loc["E6",overlap_sample].values])
    y= phos_df.loc[gene,overlap_sample].values
    
    
    model = TheilSenRegressor()
# evaluate model
    results = evaluate_model(x, y, model)
    print('Mean MAE: %.3f (%.3f)' % (mean(results), std(results)))
# plot the line of best fit
    plot_best_fit(x, y, model)
 
    
    plt.xlabel("mRNA abundance of E6")
    plt.ylabel("phosphoprotein abundance of %s"%gene)
    plt.title("spearman r:%s\np:%s"%("%.2f"%tmt_phos.loc[gene,"r_spearman"],"%.2e"%tmt_phos.loc[gene,"pvalue"]))    
    
    plt.savefig(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_phos\TheilSen_regression_%s_protein_phos_cor_with_e6.pdf"%gene,dpi=300,bbox_inches="tight")












tmt_ace = pd.read_excel(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_ace\E6_spearman_pvalue_ace.xlsx",index_col=0)

overlap_ace = list(set(tmt_ace.index)&set(tmt_protein.index))

ace_specific_gene = []
plt.figure(figsize=(5,5))
for gene in overlap_ace:
    plt.scatter(tmt_protein.loc[gene,"r_spearman"],tmt_ace.loc[gene,"r_spearman"],c="blue",s=2)
    if tmt_protein.loc[gene,"r_spearman"]<0 and tmt_ace.loc[gene,"r_spearman"]> 0.24:
        print(gene,tmt_protein.loc[gene,"r_spearman"],tmt_ace.loc[gene,"r_spearman"])
        
        ace_specific_gene.append(gene)  #H2AX,WEE1
        
for gene in ['RAC1','AP2A1','PAK2','COPA']:
    plt.annotate(gene,(tmt_protein.loc[gene,"r_spearman"],tmt_ace.loc[gene,"r_spearman"]),
                 (tmt_protein.loc[gene,"r_spearman"]+0.05,tmt_ace.loc[gene,"r_spearman"]+0.05),arrowprops=dict(facecolor='black',arrowstyle='-'))

plt.xlim([-0.5,0.5])
plt.ylim([-0.5,0.5])   
plt.xlabel("spearman_r_with_E6(protein)",size=6)     
plt.ylabel("spearman_r_with_E6(acetylsites)",size=6)
plt.savefig(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_ace\scatterplot_spearman_with_e6_ace.pdf",dpi=300)




for gene in ['RAC1','AP2A1','PAK2','COPA']:
    
    plt.figure(figsize=(8,3))
    plt.subplot(121)
 #   plt.scatter(hpv_gene_expression_tpm.loc["E6",overlap_sample],tmt_df.loc[gene,overlap_sample],s=5)
    mean_ = np.mean(hpv_gene_expression_tpm.loc["E6",overlap_sample].values)
    std_= np.std(hpv_gene_expression_tpm.loc["E6",overlap_sample].values)
    
    
    x = np.array([[(i-mean_)/std_] for i in list(hpv_gene_expression_tpm.loc["E6",overlap_sample].values)])
    y = tmt_df.loc[gene,overlap_sample].values
    plt.scatter(x,y,s=5)    
    
    
   # slope, intercept = np.polyfit(x, y, 1)
   # plt.plot(x, slope * x + intercept, color='red')
    model = TheilSenRegressor()
# evaluate model
    results = evaluate_model(x, y, model)
    print('Mean MAE: %.3f (%.3f)' % (mean(results), std(results)))
# plot the line of best fit
    plot_best_fit(x, y, model)
    
    plt.xlabel("mRNA abundance of E6")
    plt.ylabel("protein abundance of %s"%gene)
    plt.title("spearman r:%s\np:%s"%("%.2f"%tmt_protein.loc[gene,"r_spearman"],"%.2e"%tmt_protein.loc[gene,"pvalue"]))
    
    plt.subplot(122)
    #plt.scatter(hpv_gene_expression_tpm.loc["E6",overlap_sample],ace_df.loc[gene,overlap_sample],s=5)
    mean_ = np.mean(hpv_gene_expression_tpm.loc["E6",overlap_sample].values)
    std_= np.std(ace_df.loc[gene,overlap_sample].values)

    x = np.array([[(i-mean_)/std_] for i in list(hpv_gene_expression_tpm.loc["E6",overlap_sample].values)])
    y = ace_df.loc[gene,overlap_sample].values
    plt.scatter(x,y,s=5)    
    
#    x= np.array([[i] for i in hpv_gene_expression_tpm.loc["E6",overlap_sample].values])
#    y= ace_df.loc[gene,overlap_sample].values
    
    
    model = TheilSenRegressor()
# evaluate model
    results = evaluate_model(x, y, model)
    print('Mean MAE: %.3f (%.3f)' % (mean(results), std(results)))
# plot the line of best fit
    plot_best_fit(x, y, model)
 
    
    plt.xlabel("mRNA abundance of E6")
    plt.ylabel("phosphoprotein abundance of %s"%gene)
    plt.title("spearman r:%s\np:%s"%("%.2f"%tmt_ace.loc[gene,"r_spearman"],"%.2e"%tmt_ace.loc[gene,"pvalue"]))    
  
    
    plt.savefig(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_ace\%s_protein_ace_cor_with_e6.pdf"%gene,dpi=300,bbox_inches="tight")










        
protein_specific_gene = []

for gene in overlap_gene:
  #  plt.scatter(tmt_protein.loc[gene,"r_spearman"],tmt_phos.loc[gene,"r_spearman"],c="blue",s=2)
    if tmt_protein.loc[gene,"r_spearman"]>0.25 and tmt_phos.loc[gene,"r_spearman"]<0:
        print(gene,tmt_protein.loc[gene,"r_spearman"],tmt_phos.loc[gene,"r_spearman"])
        protein_specific_gene.append(gene)
        

library_list = ['KEGG_2021_Human','GO_Biological_Process_2023','The_Kinase_Library_2023']
for library in library_list:        
    gp.enrichr(phos_specific_gene,gene_sets=library,description='', outdir=r'G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\E6_correlation_phos_site\protein_specific', 
           background='hsapiens_gene_ensembl', cutoff=0.05, format='png', figsize=(8, 6), top_term=10, no_plot=False, verbose=False)

    
    

