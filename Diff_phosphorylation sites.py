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




def pcc_heatmap(cluster_2d,vmin_=-0.5,vmax_=1):
    
    sns.heatmap(cluster_2d,vmin=vmin_,vmax=vmax_,cmap='bwr',cbar=False,yticklabels=False,xticklabels=False)
    
def clusterbar(ax,sub_protein_clusters,rank_gene_list):
    bottom=0
    for key in range(1,len(sub_protein_clusters)+1):
        cluster_gene = sub_protein_clusters[key]
        ax.bar(0,-len(cluster_gene),bottom=bottom,width=1,align="edge")
        bottom -= len(cluster_gene)
    ax.set_xlim([0,3])
    ax.set_ylim([-len(rank_gene_list),0])
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

imputation_f = r"G:\cervical_cancer_multiomics\results\expression_table\tmt_phos_half_quantified_dreamai_imputation_res.csv"


def tmt_cluster(imputation_f,type_="tmt_p_pos",cluster_n=4):
    
#    imputation_f = imputation_ace_f
    clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)
    tmt_f = pd.read_csv(imputation_f,sep=",",index_col=0)
    
    normal_sample = []
    
    for column in tmt_f.columns:
        if column.startswith("N"):
            normal_sample.append(column)
        else:
            pass
    
    tmt_f_t = tmt_f.T
    
    dir_ = r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\%s"%(type_)
    mkdir_.mkdir(dir_)
    os.chdir(dir_)
    
    tmt_cluster_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\cluster_3.csv",
                              sep=",",index_col=0)
    tmt_cluster_df['type'] = clinical_df.loc[tmt_cluster_df.index,'Histological_Diagnosis']
    tmt_cluster_df = tmt_cluster_df.sort_values(by=['x','type'])
    
    tmt_f_t = tmt_f_t.loc[tmt_cluster_df.index]
    
    tmt_f_t["tmt_group"] = tmt_cluster_df.loc[tmt_f_t.index,"x"]
  
    pvalue_list = []
    
    for column in tmt_f_t.columns[:-1]:
        
        sub_df = tmt_f_t[[column,"tmt_group"]]
        sub_df.columns = ["gene","tmt_group"]
        a='gene ~ C(tmt_group)'
        
        model = ols(a, data=sub_df).fit()
        
        anova_table = sm.stats.anova_lm(model,typ=2)
        pc1_anova_pvalue = anova_table.loc['C(tmt_group)','PR(>F)']
        
        pvalue_list.append(pc1_anova_pvalue)
        
    final_df = pd.DataFrame([tmt_f_t.columns[:-1],pvalue_list]).T

    final_df.columns = ["p_pos","anova_pvalue"]
    
    final_df["BH-corrected P-value"] =[i if i<1 else 1 for i in final_df["anova_pvalue"]/final_df["anova_pvalue"].rank(axis=0)*len(final_df)]
    
    final_df["bonferroni-corrected P-value"] =[i if i<1 else 1 for i in final_df["anova_pvalue"]*len(final_df)]
    
    final_df_sort = final_df.sort_values(by="anova_pvalue")
    
    final_df_sort.index = final_df_sort["p_pos"]
    
    final_df_sort.to_csv("%s_annova_results.csv",sep="\t")
    final_df_sort.to_excel("%s_annova_results.xlsx")
    
    cuttype = "bonferroni-corrected P-value"
    
    cutoff = 0.05
    
    sig_df = final_df_sort.loc[final_df_sort[cuttype]<cutoff]
    
    
    sig_z_df_raw = tmt_f.loc[sig_df["p_pos"],tmt_cluster_df.index]
    
    sig_z_df = sig_z_df_raw.sub(sig_z_df_raw.median(axis=1), axis=0)
    
    normal_z_df_raw = tmt_f.loc[sig_z_df.index,normal_sample]
    
    normal_z_df = normal_z_df_raw.sub(sig_z_df_raw.median(axis=1),axis=0)
    
    
    
    g = sns.clustermap(sig_z_df,cmap="bwr",col_cluster=False,row_cluster=True,vmin=-2,vmax=2,metric="correlation")
    
    proteins_indexs=scipy.cluster.hierarchy.fcluster(g.dendrogram_row.calculated_linkage,t=cluster_n,criterion='maxclust')
    
    
    sub_protein_clusters={}
    for protein ,cluster in zip(sig_z_df.index,proteins_indexs):
        if cluster in sub_protein_clusters.keys():
            sub_protein_clusters[cluster].append(protein)
        else:
            sub_protein_clusters[cluster] = [protein]
    
    go_result_dict = {}
    kegg_result_dict ={}
    rank_gene_list = []
    for key in [i for i in range(1,cluster_n+1)]:
        cluster_gene_pos = sub_protein_clusters[key]
        cluster_gene_df = pd.DataFrame(cluster_gene_pos,columns=['cluster%s'%(key)])
        cluster_gene_df.to_excel(r"cluster%s_gene.xlsx"%(key))
        
        cluster_gene = list(set([i.split(":")[0] for i in cluster_gene_pos]))
        go_result = gseapy.enrichr(cluster_gene,gene_sets="GO_Biological_Process_2021",organism='human', description='', 
                                   outdir=r'%s\cluster%s'%(dir_,key), 
                                   background='hsapiens_gene_ensembl', cutoff=0.05, 
                                    format='pdf', figsize=(8, 6), top_term=10)
        kegg_result = gseapy.enrichr(cluster_gene,gene_sets="KEGG_2021_Human", organism='human', description='', 
                                   outdir=r'%s\cluster%s'%(dir_,key), 
                                   background='hsapiens_gene_ensembl', cutoff=0.05, 
                                   format='pdf', figsize=(8, 6), top_term=10)
        go_result_dict['cluster%s'%(key)] = go_result.res2d
        kegg_result_dict['cluster%s'%(key)] = kegg_result.res2d
        
        rank_gene_list+= cluster_gene_pos
    
    Z = g.dendrogram_row.linkage
    
    cluster_2d = g.data2d
    cluster_2d=cluster_2d.loc[rank_gene_list]
    
    cluster_2d.to_excel("hyper_%s_tumor_centered_z.xlsx"%(type_))
    
    import argparse
    import matplotlib.gridspec as gridspec
    from matplotlib.gridspec import GridSpec
    
    fig = plt.figure(figsize=(10,10),dpi=300)
    gs = GridSpec(60,30,figure=fig)
    
    
  #  ax2 = fig.add_subplot(gs[0:20,8:25])
    
  #  cluster_map_sample_type(cluster_2d,ax=ax2,add_feature_cluster_df=add_feature_cluster_df)
    
    ax3 = fig.add_subplot(gs[20:60,8:25])
    vmin=-2
    vmax=2
    
    pcc_heatmap(cluster_2d,vmin_=vmin,vmax_=vmax)
    
    ax4 = fig.add_subplot(gs[20:60,4:8])
 #   pcc_heatmap(tmt_z_hyper_normal_sample_df_normalized.loc[rank_gene_list],vmin_=vmin,vmax_=vmax)
    normal_df = normal_z_df.loc[cluster_2d.index]
    normal_df.to_excel("hyper_%s_normal_centered_z.xlsx"%(type_))
    
    
    pcc_heatmap(normal_df,vmin_=vmin,vmax_=vmax)
    
    
    
    ax9 = fig.add_subplot(gs[20:60,25:26])
    clusterbar(ax9,sub_protein_clusters,rank_gene_list)
    
    plt.savefig("hyper_%s_across_tmt_cluster.pdf"%(type_),dpi=300,bbox_inches="tight")
    plt.savefig("hyper_%s_across_tmt_cluster.png"%(type_),dpi=300,bbox_inches="tight")
    
    
imputation_f = r"G:\cervical_cancer_multiomics\results\expression_table\tmt_phos_half_quantified_dreamai_imputation_res.csv"
tmt_cluster(imputation_f,type_="tmt_p_pos",cluster_n=4)


imputation_ace_f = r"G:\cervical_cancer_multiomics\results\expression_table\tmt_ace_half_quantified_dreamai_imputation_res.csv"

tmt_cluster(imputation_ace_f,type_="tmt_ace_pos",cluster_n=3)    
    
    
    
imputation_tmt_f = r"G:\cervical_cancer_multiomics\results\expression_table\tmt_pro_half_quantified_dreamai_imputation_res.csv"

tmt_cluster(imputation_tmt_f,type_="tmt_protein",cluster_n=4)    



imputation_dia_f = r"G:\cervical_cancer_multiomics\results\expression_table\dia_pro_half_quantified_dreamai_imputation_res.csv"

tmt_cluster(imputation_dia_f,type_="dia_protein",cluster_n=4) 










def rna_cluster(imputation_f,type_="tmt_p_pos",cluster_n=4):
    
#    imputation_f = imputation_ace_f
    clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)
    tmt_f = pd.read_csv(imputation_f,sep="\t",index_col=0)
    
    normal_sample = []
    
    for column in tmt_f.columns:
        if column.startswith("N"):
            normal_sample.append(column)
        else:
            pass
    
    tmt_f_t = tmt_f.T
    
    dir_ = r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\%s"%(type_)
    mkdir_.mkdir(dir_)
    os.chdir(dir_)
    
    tmt_cluster = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\cluster_3.csv",
                              sep=",",index_col=0)
    
    
    overlap_sample = list(set(tmt_cluster.index)&set(tmt_f.columns))
    
    tmt_cluster = tmt_cluster.loc[overlap_sample]
    
    tmt_cluster['type'] = clinical_df.loc[tmt_cluster.index,'Histological_Diagnosis']
    tmt_cluster = tmt_cluster.sort_values(by=['x','type'])
    
    
    
    
    tmt_f_t = tmt_f_t.loc[tmt_cluster.index]
    
    tmt_f_t["tmt_group"] = tmt_cluster.loc[tmt_f_t.index,"x"]
  
    pvalue_list = []
    
    for column in tmt_f_t.columns[:-1]:
        
        sub_df = tmt_f_t[[column,"tmt_group"]]
        sub_df.columns = ["gene","tmt_group"]
        a='gene ~ C(tmt_group)'
        
        model = ols(a, data=sub_df).fit()
        
        anova_table = sm.stats.anova_lm(model,typ=2)
        pc1_anova_pvalue = anova_table.loc['C(tmt_group)','PR(>F)']
        
        pvalue_list.append(pc1_anova_pvalue)
        
    final_df = pd.DataFrame([tmt_f_t.columns[:-1],pvalue_list]).T

    final_df.columns = ["p_pos","anova_pvalue"]
    
    final_df["BH-corrected P-value"] =[i if i<1 else 1 for i in final_df["anova_pvalue"]/final_df["anova_pvalue"].rank(axis=0)*len(final_df)]
    
    final_df["bonferroni-corrected P-value"] =[i if i<1 else 1 for i in final_df["anova_pvalue"]*len(final_df)]
    
    final_df_sort = final_df.sort_values(by="anova_pvalue")
    
    final_df_sort.index = final_df_sort["p_pos"]
    
    final_df_sort.to_csv("%s_annova_results.csv"%type_,sep="\t")
    final_df_sort.to_excel("%s_annova_results.xlsx"%type_,)
    
    cuttype = "bonferroni-corrected P-value"
    
    cutoff = 0.05
    
    sig_df = final_df_sort.loc[final_df_sort[cuttype]<cutoff]
    
    
    sig_z_df_raw = tmt_f.loc[sig_df["p_pos"],tmt_cluster.index]
    
    sig_z_df = sig_z_df_raw.sub(sig_z_df_raw.median(axis=1), axis=0)
    
    normal_z_df_raw = tmt_f.loc[sig_z_df.index,normal_sample]
    
    normal_z_df = normal_z_df_raw.sub(sig_z_df_raw.median(axis=1),axis=0)
    
    
    
    g = sns.clustermap(sig_z_df,cmap="bwr",col_cluster=False,row_cluster=True,vmin=-2,vmax=2,metric="correlation",z_score=0)
    
    proteins_indexs=scipy.cluster.hierarchy.fcluster(g.dendrogram_row.calculated_linkage,t=cluster_n,criterion='maxclust')
    
    
    sub_protein_clusters={}
    for protein ,cluster in zip(sig_z_df.index,proteins_indexs):
        if cluster in sub_protein_clusters.keys():
            sub_protein_clusters[cluster].append(protein)
        else:
            sub_protein_clusters[cluster] = [protein]
    
    go_result_dict = {}
    kegg_result_dict ={}
    rank_gene_list = []
    for key in [i for i in range(1,cluster_n+1)]:
        cluster_gene_pos = sub_protein_clusters[key]
        cluster_gene_df = pd.DataFrame(cluster_gene_pos,columns=['cluster%s'%(key)])
        cluster_gene_df.to_excel(r"cluster%s_gene.xlsx"%(key))
        
        cluster_gene = list(set([i.split(":")[0] for i in cluster_gene_pos]))
        go_result = gseapy.enrichr(cluster_gene,gene_sets="GO_Biological_Process_2021",organism='human', description='', 
                                   outdir=r'%s\cluster%s'%(dir_,key), 
                                   background='hsapiens_gene_ensembl', cutoff=0.05, 
                                    format='pdf', figsize=(8, 6), top_term=10)
        kegg_result = gseapy.enrichr(cluster_gene,gene_sets="KEGG_2021_Human", organism='human', description='', 
                                   outdir=r'%s\cluster%s'%(dir_,key), 
                                   background='hsapiens_gene_ensembl', cutoff=0.05, 
                                   format='pdf', figsize=(8, 6), top_term=10)
        go_result_dict['cluster%s'%(key)] = go_result.res2d
        kegg_result_dict['cluster%s'%(key)] = kegg_result.res2d
        
        rank_gene_list+= cluster_gene_pos
    
    Z = g.dendrogram_row.linkage
    
    cluster_2d = g.data2d
    cluster_2d=cluster_2d.loc[rank_gene_list]
    
    cluster_2d.to_excel("hyper_%s_tumor_centered_z.xlsx"%(type_))
    
    import argparse
    import matplotlib.gridspec as gridspec
    from matplotlib.gridspec import GridSpec
    
    fig = plt.figure(figsize=(10,10),dpi=300)
    gs = GridSpec(60,30,figure=fig)
    
    
  #  ax2 = fig.add_subplot(gs[0:20,8:25])
    
  #  cluster_map_sample_type(cluster_2d,ax=ax2,add_feature_cluster_df=add_feature_cluster_df)
    
    ax3 = fig.add_subplot(gs[20:60,8:25])
    vmin=-2
    vmax=2
    
    pcc_heatmap(cluster_2d,vmin_=vmin,vmax_=vmax)
    
    ax4 = fig.add_subplot(gs[20:60,4:8])
 #   pcc_heatmap(tmt_z_hyper_normal_sample_df_normalized.loc[rank_gene_list],vmin_=vmin,vmax_=vmax)
    normal_df = normal_z_df.loc[cluster_2d.index]
    normal_df.to_excel("hyper_%s_normal_centered_z.xlsx"%(type_))
    
    
    pcc_heatmap(normal_df,vmin_=vmin,vmax_=vmax)
    
    
    
    ax9 = fig.add_subplot(gs[20:60,25:26])
    clusterbar(ax9,sub_protein_clusters,rank_gene_list)
    
    plt.savefig("hyper_%s_across_tmt_cluster.pdf"%(type_),dpi=300,bbox_inches="tight")
    plt.savefig("hyper_%s_across_tmt_cluster.png"%(type_),dpi=300,bbox_inches="tight")


rna_f = r"G:\cervical_cancer_multiomics\results\expression_table\tpm_remove_error_sample.xls"

rna_cluster(rna_f,type_="rna_cluster",cluster_n=2) 



    
        
   
    
