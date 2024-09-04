# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 21:10:55 2022

@author: guixiuqi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import scipy.stats as stats
import seaborn as sns
from scipy import stats
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
from scipy.stats import chi2_contingency
import statsmodels.api as sm
from statsmodels.formula.api import ols
import gseapy as gp
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import matplotlib.cm as cm

def mkdir(path):
    
    folder = os.path.exists(path)
    
    if not folder:
        os.makedirs(path)
    else:
        pass
    

def drew_bar_HPV_gene_expression(ax,vmin=-1,vmax=1,gene=""):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.OrRd,norm=pc1_norm,ticks=[vmin,vmax],orientation='horizontal')
    ax.set_title("log$_2$(TPM+1) of %s"%(gene),size=3)
    cb1.ax.set_xticklabels([str(vmin),str(vmax)],rotation=0,fontdict={'fontsize':3}) 


def get_color_list(cmap=cm.bwr,value_list=[],vmin=-3,vmax=8):
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    pc1_cmap = cmap
    pc1_m=cm.ScalarMappable(norm=pc1_norm,cmap=pc1_cmap)
    pc1_color_list = []
    for pc in value_list:
        pc1_color_list.append(pc1_m.to_rgba(pc))
    return pc1_m,pc1_color_list

def hpv_gene_expression_color(gene_epxression_df,ax,vmin=-1,vmax=1):  
    
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    hight=1
    
   # vmin=min(arm_sort_df['3q'])
   # vmax=max(arm_sort_df['3q'])
    pc1_m,pc1_color_list = get_color_list(cmap=cm.bwr,value_list=gene_epxression_df.iloc[:,0].values,vmin=vmin,vmax=vmax)  
    plt.bar(range(0,len(gene_epxression_df)),[hight]*len(gene_epxression_df),width=1,align="edge",color = pc1_color_list)
    
    plt.xlim([0,len(gene_epxression_df)])


def pcc_heatmap(cluster_2d,vmin_=-0.5,vmax_=1,y_ticks="",yticklabels=True):
    g = sns.heatmap(cluster_2d,vmin=vmin_,vmax=vmax_,cmap='bwr',cbar=False,yticklabels=yticklabels,xticklabels=False)
    
    #g = sns.clustermap(cluster_2d,vmin=vmin_,vmax=vmax_,cmap='bwr',cbar=False,yticklabels=yticklabels,xticklabels=False,
    #                   col_cluster=False,row_cluster=False,col_colors=col_colors)
    g.tick_params(right=False,left=False, top=False, labelright=True, labeltop=False,labelleft=False, labelrotation=0)
    plt.plot([0,0],)
    plt.ylabel(y_ticks,size=5)
    
    
def drew_bar_z(ax,vmin = -3,vmax = 3):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.bwr,norm=pc1_norm,ticks=[vmin,vmax],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("z-statistic",size=3)
    cb1.ax.set_xticklabels([str(vmin),str(vmax)],rotation=0,fontdict={'fontsize':3})    


def e6_pathway_gene_heatmap(expression_df,gene,data_,r_cutoff):
    
    
    outdir = r"G:\cervical_cancer_multiomics\results\KSEA\PRKCB\%s_correlation_%s\\"%(gene,data_)
    os.chdir(outdir)
    pathway_f = r"G:\cervical_cancer_multiomics\results\KSEA\PRKCB\%s_correlation_%s\gseapy.prerank.gene_sets.report.csv"%(gene,data_)
    
    spearmanr_f = r"G:\cervical_cancer_multiomics\results\KSEA\PRKCB\%s_correlation_%s\%s_spearman_pvalue_%s.xlsx"%(gene,data_,gene,data_)

    spearmanr_df = pd.read_excel(spearmanr_f,index_col=0)

    heatmap_gene_list =   list(spearmanr_df.loc[(spearmanr_df['r_spearman']>r_cutoff)|(spearmanr_df['r_spearman']<-r_cutoff)].index)

    heatmap_gene_list = list(set(heatmap_gene_list))      
    
    sub_spearmanr_df = spearmanr_df.loc[heatmap_gene_list]
    
    sub_spearmanr_df = sub_spearmanr_df.sort_values(by="r_spearman",ascending=False)
    
    heatmap_gene_list = list(sub_spearmanr_df.index)

    gene_epxression_df = tmt_df.loc[gene,tumor_sample_type_list].to_frame()
    
    gene_epxression_df = gene_epxression_df.sort_values(by=gene)
    
    expression_df = expression_df.loc[heatmap_gene_list,list(gene_epxression_df.index)]
    
    if data_=="rna":
        g = sns.clustermap(expression_df,z_score=0,col_cluster=False,row_cluster=False)
        
        expression_df = g.data2d
    

    fig = plt.figure(figsize=(3,4),dpi=300)
    gs = GridSpec(30,30,figure=fig)
    
    ax3 = fig.add_subplot(gs[1:2,10:20])

    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)            
    vmin_=-2
    vmax_ = 2
    hpv_gene_expression_color(gene_epxression_df,ax3,vmin=vmin_,vmax=vmax_)
                
    ax3 = fig.add_subplot(gs[0:1,10:20])
    
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False) 
    ax3.set_xticks([])
    ax3.set_yticks([])

   
    ax3.text(0,0,"%s_correlated_"%(gene)+"(%s)"%(data_),size=6)
    


    ax = fig.add_subplot(gs[2:20,10:20])

    vmin=-2
    vmax=2
    pcc_heatmap(expression_df,vmin_=vmin,vmax_=vmax,y_ticks="",yticklabels=False)

    plt.rcParams['xtick.labelsize']=6
    plt.rcParams['ytick.labelsize']=6    
    plt.rcParams['ytick.direction'] = 'in'     

    
    ax = fig.add_subplot(gs[22:23,12:14])
    drew_bar_z(ax,vmin = vmin,vmax = vmax)
    
    plt.savefig(outdir+"%s_correlated_spearmanr_%s_pathway_heatmap_scatter.pdf"%(gene,r_cutoff),dpi=300,bbox_inches="tight")
    plt.savefig(outdir+"%s_correlated_spearmanr_%s_pathway_heatmap_scatter.png"%(gene,r_cutoff),dpi=300,bbox_inches="tight")
    
    
    pos_gene_list =   list(spearmanr_df.loc[(spearmanr_df['r_spearman']>r_cutoff)].index)
    
    neg_gene_list =   list(spearmanr_df.loc[(spearmanr_df['r_spearman']<-r_cutoff)].index)
    
    dir_ = outdir
    
    key_list = ["pos","neg"]
    i=0
    for type_ in [pos_gene_list,neg_gene_list]:
        
        key = key_list[i]

        cluster_gene_pos = type_
        cluster_gene_df = pd.DataFrame(cluster_gene_pos,columns=['cluster_%s'%(key)])
        cluster_gene_df.to_excel(r"cluster_%s_gene.xlsx"%(key))
        
        cluster_gene = list(set([i.split(":")[0] for i in cluster_gene_pos]))
        go_result = gseapy.enrichr(cluster_gene,gene_sets="GO_Biological_Process_2021",organism='human', description='', 
                                   outdir=r'%s\cluster_%s'%(dir_,key), 
                                   background='hsapiens_gene_ensembl', cutoff=0.05, 
                                    format='png', figsize=(8, 6), top_term=10)
        kegg_result = gseapy.enrichr(cluster_gene,gene_sets="KEGG_2021_Human", organism='human', description='', 
                                   outdir=r'%s\cluster_%s'%(dir_,key), 
                                   background='hsapiens_gene_ensembl', cutoff=0.05, 
                                   format='png', figsize=(8, 6), top_term=10)
        
        i+=1
    


if __name__=="__main__":
    
    
    os.chdir(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression")

    
    dia_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\dia_pro_half_quantified_dreamai_imputation_res.csv",
                             sep=",",index_col=0)
    
    tmt_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_pro_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)
    
    tmt_df.index = [index.replace("-","_").replace("/","_") for index in tmt_df.index]
    
    rna_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\fpkm_remove_error_sample.xls",sep="\t",index_col=0)
    
    rna_df.index = [index.replace('.',"_").replace("-","_") for index in rna_df.index]
    rna_df = np.log2(rna_df+1)
    
    g = sns.clustermap(rna_df,z_score=1)


    phos_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_protein_level_phos_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)
    
    phos_df.index= [index.replace("-","_").replace("/","_") for index in phos_df.index]
    
    ace_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_protein_level_ace_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)
    
    ace_df.index= [index.replace("-","_").replace("/","_") for index in ace_df.index]
    
    phos_p_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_phos_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)
    
    phos_p_df.index= [index.replace("-","_").replace("/","_") for index in phos_p_df.index]

    

    
    common_sample = set(dia_df.columns)&set(tmt_df.columns)&set(rna_df.columns)
    
    tumor_sample_type_list=[]
    for sample in common_sample:
        if sample.startswith("T"):
            tumor_sample_type_list.append(sample)
    


    os.chdir(r"G:\cervical_cancer_multiomics\results\KSEA\PRKCB")

    e6_pathway_gene_heatmap(tmt_df,gene="PRKCB",data_="tmt",r_cutoff=0.3)
    
    
    e6_pathway_gene_heatmap(dia_df,gene="PRKCB",data_="dia",r_cutoff=0.3)
    
    
    e6_pathway_gene_heatmap(phos_df,gene="PRKCB",data_="phos",r_cutoff=0.3)

    e6_pathway_gene_heatmap(rna_df,gene="PRKCB",data_="rna",r_cutoff=0.3)
    
    e6_pathway_gene_heatmap(phos_p_df,gene="PRKCB",data_="phosp",r_cutoff=0.3)
    
    
    

        
        
        

    
    





 