# -*- coding: utf-8 -*-
"""
Created on Tue May 24 00:21:05 2022

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
import glob
import seaborn as sns
import os
import scipy
import collections
os.chdir(r"G:\cervical_cancer_multiomics\results\cnvkit\\")
chro_l_df = pd.read_csv("GRCh38_chromosome_length.txt",sep="\t",header=None,index_col=0)

chromosome_x_dict = {}
x=0

for i in range(1,23):
    
    chro="chr"+str(i)
    chro_l = chro_l_df.loc[chro,1]

    chromosome_x_dict[chro]=x
    x+=chro_l
chromosome_x_dict["chrX"] = x

x_end = chromosome_x_dict["chrX"]+chro_l_df.loc["chrX",1]

gene_location_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\WES\pQTL\gene_location.txt",sep="\t",header=None,index_col=0)



def drew_spearmans_heatmap(spearman_df,figname):
    

     
    gene_list = []

    x_axis_list = [] 
    x_width_list = []
    
    positive_n_list =[]
    negative_n_list = []
    chro_list = ["chr"+str(i) for i in range(1,23)]+["chrX"]
    for chro in chro_list:
        chro_df = gene_location_df.loc[gene_location_df[1]==chro]
        chro_df = chro_df.sort_values(by=2)
        for index in chro_df.index:
            if index in spearman_df.index:
                gene_list.append(index)
                start = chro_df.loc[index,2]
                end=chro_df.loc[index,3]
                length=end-start
                x_axis_list.append(chromosome_x_dict[chro]+start)
                x_width_list.append(length)
                gene_spearman_df = list(spearman_df.loc[index].values)
                gene_dict =collections.Counter(gene_spearman_df)
                
                if 1 in gene_dict.keys():
                    positive_n_list.append(gene_dict[1])
                else:
                    positive_n_list.append(0)
                if -1 in gene_dict.keys():
                    negative_n_list.append(-gene_dict[-1])
                else:
                    negative_n_list.append(0)
            else:
                pass
                    
    fig = plt.figure(figsize=(10,3))
    ax=fig.gca()
    plt.bar(x_axis_list,positive_n_list,width=x_width_list,color="red",linewidth=0,align="edge")
    plt.bar(x_axis_list,negative_n_list,width=x_width_list,color="green",linewidth=0,align="edge")
    plt.xlim([0,x_end])
    ax.spines['top'].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks([])
    
    all_df = pd.DataFrame([gene_list,[gene_location_df.loc[i,1] for i in gene_list],[gene_location_df.loc[i,2] for i in gene_list],positive_n_list,[-i for i in negative_n_list]]).T
    
    all_df.to_excel(figname.split(".")[0]+".xlsx")
    

    for key in list(chromosome_x_dict.keys())[1:]:
        plt.plot([chromosome_x_dict[key],chromosome_x_dict[key]],[-200,200],linestyle='--',color="grey")
    plt.plot([x_end,x_end],[-200,200],linestyle='--',color="grey")
  #  plt.ylim([-200,200])
    plt.ylabel("Number of \nsignificant correlations",size=10)
        
    for key in chromosome_x_dict.keys():
        plt.text(chromosome_x_dict[key]+chro_l_df.loc[key,1]/2,-200,key.strip("chr"),horizontalalignment='center')
    
    
#    plt.savefig(figname,dpi=300,bbox_inches="tight")
    
os.chdir(r"G:\cervical_cancer_multiomics\Code\Fig2")           
     
spearman_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\cnvkit\tumor_gistic2\multiOmicsVic\cna_dia_correlation_fdr_0.01.csv",
                          index_col=0,sep="\t")            
drew_spearmans_heatmap(spearman_df,"cna_dia_correlation_fdr_0.01_cis_trans_number.pdf")


spearman_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\cnvkit\tumor_gistic2\multiOmicsVic\cna_tmt_correlation_fdr_0.01.csv",
                          index_col=0,sep="\t")            
drew_spearmans_heatmap(spearman_df,"cna_tmt_correlation_fdr_0.01_cis_trans_number.pdf")


spearman_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\cnvkit\tumor_gistic2\multiOmicsVic\cna_rna_correlation_fdr_0.01.csv",
                          index_col=0,sep="\t")            
drew_spearmans_heatmap(spearman_df,"cna_rna_correlation_fdr_0.01_cis_trans_number.pdf")


spearman_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\cnvkit\tumor_gistic2\multiOmicsVic\cna_rna_correlation_fdr_0.001.csv",
                          index_col=0,sep="\t")            
drew_spearmans_heatmap(spearman_df,"cna_rna_correlation_fdr_0.001_cis_trans_number.pdf")





spearman_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\cnvkit\tumor_gistic2\multiOmicsVic\cna_Phosphoprotein_correlation_fdr_0.01.csv",
                          index_col=0,sep="\t")            
drew_spearmans_heatmap(spearman_df,"cna_Phosphoprotein_correlation_fdr_0.01_cis_trans_number.pdf")



spearman_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\cnvkit\tumor_gistic2\multiOmicsVic\cna_rna_correlation_fdr_0.01_FPKM.csv",
                          index_col=0,sep="\t")            
drew_spearmans_heatmap(spearman_df,"cna_rna_correlation_fdr_0.01_cis_trans_number_FPKM.pdf")




