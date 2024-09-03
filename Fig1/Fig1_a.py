# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 20:19:03 2022

@author: guixiuqi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Agg') 
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] ='Arial'
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
import os

import collections

os.chdir(r"G:\cervical_cancer_multiomics\results\WES\mutation_number_per_sample")

def rank_gene(gene_list):
    mutation_f = r"G:\cervical_cancer_multiomics\results\WES\mutation_number_per_sample\tumor_gene_mutation_0_1.xlsx"

    mutation_df = pd.read_excel(mutation_f,index_col=0)
    mutation_df = mutation_df.drop(columns=['TB113-1',"TB83"])


    for gene in gene_list:
        if gene not in mutation_df.index:
            print (gene)

    g_mutation_df = mutation_df.loc[gene_list]
    g_mutation_df["#somaticmutation"] = g_mutation_df.sum(axis=1)
    
    g_mutation_df_sort = g_mutation_df.sort_values(by="#somaticmutation",ascending= False)
    gene_rank = list(g_mutation_df_sort.index)
    
    final_mutation_df = g_mutation_df_sort.iloc[:,:-1]
    
    final_mutation_df_sort = final_mutation_df.sort_values(by=gene_rank,ascending=False,axis=1)
    
    return final_mutation_df_sort






mutation_color_dict = {'nonframeshift deletion':'blue',
 'frameshift deletion':'orange',
 'stopgain':'brown',
 'nonsynonymous SNV':'c',
 'stoploss':'r',
 'synonymous SNV':'gray',
 'nonframeshift substitution':'violet',
 'nonframeshift insertion':'y',
 'frameshift insertion':'g',
 'startloss':'lightcoral',
 'unknown':'silver',
 'Splicing':'teal',
 'WT':'white'}

#tumor_maf_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\WES\pQTL\coding_region_mutation_Tumor_annovar_remove_B_113_1_B_75_remove_ncRNA.maf",sep="\t")

def get_gene_mutation_type(gene_list=[]):
    
    tumor_maf_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\WES\remove_error_tumor_mutation_annovar_for_mutsigCV_coding.maf",sep="\t") #only 133samples, no mutation detected in TB154 sample
    
    non_silent_df = tumor_maf_df.loc[tumor_maf_df['Variant_Classification']!='Silent']
    
    collections.Counter(non_silent_df.loc[i,'Tumor_Sample_Barcode']+"_"+non_silent_df.loc[i,'Gene.refGene'] for i in non_silent_df[['Tumor_Sample_Barcode','Gene.refGene']].index)
  #  samples = list(set(tumor_maf_df['Tumor_Sample_Barcode']))
        
    index_list = []
    
    for index in tumor_maf_df.index:
        gene = tumor_maf_df.loc[index,'Hugo_Symbol']
        if gene in gene_list:
            index_list.append(index)
    sub_maf_df = tumor_maf_df.loc[index_list]

    
    gene_list = final_mutation_df_sort.index
    sample_list = final_mutation_df_sort.columns
    
    
    _1d_all_mutation_type=[]
    _2d_mutation_list = []
    for index in gene_list:
        mutation_list = []
        for sample in sample_list:
            gene_sample_df = sub_maf_df.loc[(sub_maf_df['Tumor_Sample_Barcode']==sample)&(sub_maf_df['Hugo_Symbol']==index)]
            if len(gene_sample_df)>0:
                g_mutation_list = list(set(list(gene_sample_df['ExonicFunc.refGene'].values)))
                print(g_mutation_list)
                
                if len(g_mutation_list) >1:
                    if 'synonymous SNV' in g_mutation_list:
                        g_mutation_list.remove('synonymous SNV')
                #else:
                mutation = ";".join(g_mutation_list).replace('.','Splicing')
                mutation_list.append(mutation)
                _1d_all_mutation_type.append(mutation)
            else:
                mutation_list.append("WT")
                _1d_all_mutation_type.append('WT')
        _2d_mutation_list.append(mutation_list)
    
    rank_mutation_df = pd.DataFrame(_2d_mutation_list,columns=sample_list,index = gene_list)
    
    mutation_type_list = list(set(_1d_all_mutation_type))
    
    return rank_mutation_df


def get_mutation_ratio(final_mutation_df_sort):
    
    
    mutation_ratio_dict = {}
    mutation_df_add_ratio = final_mutation_df_sort.copy()
    sample_len = len(mutation_df_add_ratio.columns)
    mutation_df_add_ratio['sum'] = mutation_df_add_ratio.sum(axis=1)
    for index in mutation_df_add_ratio.index:
        ratio = '%.1f'% (mutation_df_add_ratio.loc[index,'sum']/sample_len*100)+"%"
        mutation_ratio_dict[index]=ratio
        
    return mutation_ratio_dict
    
    


    
def mutation_waterfall_plt(rank_mutation_df,ax):
    
    ax.set_xticks([])        
    ax.set_yticks([])
    width=1
    hight = -4
    hight_interval = -1
    bottom=-1
    
    num=0
    pos_list = []
    pos = bottom+0.5*hight
    for gene in rank_mutation_df.index:
        
        i=0
        for mutation_type in rank_mutation_df.loc[gene].values:
            mutation_color = [mutation_color_dict[m] for m in mutation_type.split(";")]

            if len(mutation_color)==1:
                plt.bar(i,hight,bottom = -1 + num*(hight+hight_interval),width=width,
                    color=mutation_color,align='edge')
            else:
                m=0
                for sub_mutation_color in mutation_color:
                    plt.bar(i,hight/len(mutation_color),bottom = -1 + num*(hight+hight_interval) + m * hight/len(mutation_color),
                    width=width,color=sub_mutation_color,align='edge')
                    m +=1
            i += 1
        num+=1
        bottom =-1 + num*(hight+hight_interval)
        pos_list.append(pos)
        pos = bottom+0.5*hight
        #pos_list.append(pos)
    ax.set_xlim([0,len(rank_mutation_df.columns)])
    ax.set_ylim([bottom-hight_interval,-1])
    for i in range(0,len(rank_mutation_df.index)):
        gene = rank_mutation_df.index[i]
        pos = pos_list[i]
        plt.text(-1,pos,gene+"("+mutation_ratio_dict[gene]+")",size=8,horizontalalignment="right",verticalalignment="center")
    
def mutation_type_square(ax):
    
#    fig = plt.figure(figsize=(20,1))
#    ax = fig.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    
    i=0
    start=0
    interval = 4
    bar_width=0.5
    pos=4
    for m_type in ['nonframeshift deletion', 'frameshift deletion',
                   'nonframeshift substitution', 'nonframeshift insertion',
                   'frameshift insertion','synonymous SNV']:
        plt.bar([i*interval+start],[4],bottom=pos-9,width=bar_width,color=mutation_color_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,pos-3,m_type.replace(' ','\n'),size=10,verticalalignment="top")
        i+=1
    plt.xlim([0,i*interval+start+bar_width])
    j=0
    for m_type in ['stopgain', 'nonsynonymous SNV', 'stoploss','startloss', 'Splicing']:
    
        plt.bar([j*interval+start],[4],bottom=pos-15,width=bar_width,color=mutation_color_dict[m_type],align='edge')
        plt.text(j*interval+start+bar_width,pos-12,m_type.replace(' ','\n'),size=10,verticalalignment="top")
        j+=1
    
def mutation_po_num_per_sample(ax):

    sample_list = final_mutation_df_sort.columns
    
    synonymous_num_list = []
    non_synonymous_num_list = []
    for sample in sample_list:
        sample_maf_df = tumor_maf_df.loc[tumor_maf_df['Tumor_Sample_Barcode']==sample]
        
        silent_df = sample_maf_df.loc[sample_maf_df['ExonicFunc.refGene']=='synonymous SNV']
        non_silent_df = sample_maf_df.loc[sample_maf_df['ExonicFunc.refGene']!='synonymous SNV']
        silent_num = len(silent_df)
        non_silent_num=len(non_silent_df)
        
        synonymous_num_list.append(silent_num)
        non_synonymous_num_list.append(non_silent_num)
        
    mutation_num_df = pd.DataFrame([sample_list,synonymous_num_list,non_synonymous_num_list]).T
    mutation_num_df.index = mutation_num_df[0]
    
    mutation_num_df = mutation_num_df.drop(columns=[0])
    
    mutation_num_df.columns = ["Synonymous","Non synonymous"]
    mutation_num_df.to_excel("synonymous_non_synonymous_number.xlsx")
    
        
#    fig = plt.figure(figsize=(22,3))
#    ax=fig.gca()
    ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    
    ax.bar(range(0,len(sample_list)),non_synonymous_num_list,label="Non-synonymous",align='edge')
    ax.bar(range(0,len(sample_list)),synonymous_num_list,bottom=non_synonymous_num_list,
           label="Synonymous",align='edge')
    ax.set_xlim([0,len(sample_list)])

    ax.set_ylabel("Mutations",size=9)
    

import matplotlib as mpl
import matplotlib.cm as cm

def get_color_list(cmap=cm.bwr,value_list=[],vmin=-3,vmax=8):
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    pc1_cmap = cmap
    pc1_m=cm.ScalarMappable(norm=pc1_norm,cmap=pc1_cmap)
    pc1_color_list = []
    for pc in value_list:
        pc1_color_list.append(pc1_m.to_rgba(pc))
    return pc1_m,pc1_color_list

#HPV_type_dict = {"HPV16":"A9",
#                 "HPV18":"A7",
#                 "HPV33":"A9",
#                 'HPV30':'A6',
#                 'HPV31':'A9',
#                 'HPV35':'A9',
#                 'HPV45':'A7',
#                 'HPV52':'A9',
#                 'HPV53':'A6',
#                 'HPV56':'A6',
#                 'HPV58':'A9',
#                 'HPV59':'A7',
#                 'HPV68':'A7',
#                 'HPV73':'A11'}
#
#
#
#def cllinical_add_a7_a9():
#    
#    sample_sub_clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_add_tumor_purity_estimation_2022_6_30_check_hpv.xlsx",
#                                           index_col=0)
#    
#    hpv_clade_list = []
#    for index in sample_sub_clinical_df.index:
#        main_hpv=sample_sub_clinical_df.loc[index,"main_HPV"]
#        if ";" in main_hpv:
#            hpv_clade=[]
#            hpv_list = main_hpv.split(";")
#            for hpv in hpv_list:
#                clade = HPV_type_dict[hpv]
#                hpv_clade.append(clade)
#            hpv_type = ";".join(sorted(hpv_clade))
#        elif main_hpv in HPV_type_dict.keys():
#            hpv_type = HPV_type_dict[main_hpv]
#        else:
#            hpv_type=main_hpv
#        hpv_clade_list.append(hpv_type)
#    
#    sample_sub_clinical_df["Final clade"]=hpv_clade_list
#    sample_sub_clinical_df.to_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_add_tumor_purity_estimation_2022_6_30_check_hpv_add_clade.xlsx")
#                
#    
#cllinical_add_a7_a9()    
    
    


def clinical_bar(ax):
    
    sample_list =  final_mutation_df_sort.columns
    sample_sub_clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_9_21.xlsx",
                                           index_col=0)
    final_sample_sub_clinical_df = sample_sub_clinical_df.loc[sample_list]
    
    stage_colot_dict = {'I':'mistyrose', 
                        'II':'coral',
                        'III':'brown',
                        'IV':'brown'}
    
    cancer_type_color_list = {'Cervical Squamous Cell Carcinoma':'tab:blue',
                        'Cervical Adenocarcinoma':'tab:orange',
                        'Adenosquamous':'tab:green',
                        'Cervical Small Cell Carcinoma ':'tab:red',
                        'Other':'tab:purple'}
        
    diff_degree_dict = {'G1':'g', 'G2':'b', 'G3':'r','Unknown':"white"}
    
    hpv_color_dict = {'A9':"b",'A7':'r',"A6":"c","A11":"m",'Negative':'grey','Unknown':'white'}
    
    
  #  hpv_color_dict = {'HPV16':"lightcoral",'HPV18':'lightseagreen','Negative':'grey','Other':'orange'}
    
    
#    fig=plt.figure(figsize=(24,5))
#    ax=fig.gca()
    y_position_list = []
    
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    hight=-4
    width=1
    hight_interval = -1
    bottom=-1
    
    
    #Histology
    i=0
    for histology_type in final_sample_sub_clinical_df['Histological_Diagnosis'].values:
        cancer_color = [cancer_type_color_list[m] for m in histology_type.split(";")]
        if len(cancer_color)==1:
            plt.bar(i,hight,bottom = bottom,width=width,
                color=cancer_color,align='edge')
        else:
            m=0
            
            for sub_mutation_color in cancer_color:
                plt.bar(i,hight/len(cancer_color),bottom = bottom + m * hight/len(cancer_color),
                width=width,color=sub_mutation_color,align='edge')
                m +=1
        i += 1
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
    
    y_position_list.append(pos)
    
    
    #stage
    stage_list = list(final_sample_sub_clinical_df['Final stage'].values)
    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,align="edge",
         bottom =[bottom]*len(sample_list),color = [stage_colot_dict[i] for i in stage_list])
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)

    y_position_list.append(pos)
    
    #differentiation degree
    
    stage_list = list(final_sample_sub_clinical_df['Grade'].replace(np.nan,"Unknown"))
    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,align="edge",
         bottom =[bottom]*len(sample_list),color = [diff_degree_dict[i] for i in stage_list])
    
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
#    pos= bottom+0.5*hight
    y_position_list.append(pos)
    
    #HPV
    hpv_all_sample = list(final_sample_sub_clinical_df['Final clade'].values)
    i=0
    for hpv_type in hpv_all_sample:
        hpv_color=[]
        hpv_list = hpv_type.split(";")
        for h in hpv_list:
            if h in hpv_color_dict.keys():
                hpv_color.append(hpv_color_dict[h])
            else:
                hpv_color.append(hpv_color_dict['Other'])
        if len(hpv_color)==1:
            plt.bar(i,hight,bottom = bottom,width=width,
                color=hpv_color,align='edge')
        else:
            m=0    
            for sub_hpv_color in hpv_color:
                plt.bar(i,hight/len(hpv_color),bottom = bottom + m * hight/len(hpv_color),
                width=width,color=sub_hpv_color,align='edge')
                m +=1
        i += 1
    
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
    #pos= bottom+0.5*hight
    y_position_list.append(pos)
    #age
#    pc1_m,pc1_color_list = get_color_list(cmap=cm.bwr,value_list=final_sample_sub_clinical_df["Age"].values,vmin=30,vmax=80)  
#    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
#         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
#    
#    pos= bottom+0.5*hight
#    bottom =bottom+(hight+hight_interval)
##    pos= bottom+0.5*hight
#    y_position_list.append(pos)
#    
#    
#    
#    #CA125
#    pc1_m,pc1_color_list = get_color_list(cmap=cm.bwr,value_list=final_sample_sub_clinical_df["CA125(U/ml)"].values,vmin=0,vmax=50)  
#    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
#         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
#    
#    pos= bottom+0.5*hight
#    bottom =bottom+(hight+hight_interval)
##    pos= bottom+0.5*hight
#    y_position_list.append(pos)
#    
#    #CA199
#    pc1_m,pc1_color_list = get_color_list(cmap=cm.bwr,value_list=final_sample_sub_clinical_df["CA199(U/ml)"].values,vmin=0,vmax=40)  
#    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
#         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
#    pos= bottom+0.5*hight
#    bottom = bottom+(hight+hight_interval)
#  #  pos= bottom+0.5*hight
#    y_position_list.append(pos)
#    
#    #SCC
#    pc1_m,pc1_color_list = get_color_list(cmap=cm.bwr,value_list=final_sample_sub_clinical_df["SCC(ng/ml)"].values,vmin=0,vmax=10)  
#    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
#         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
#    pos= bottom+0.5*hight
#    bottom =bottom+(hight+hight_interval)
#  #  pos= bottom+0.5*hight
#    y_position_list.append(pos)
#    
#    #HGB
#    pc1_m,pc1_color_list = get_color_list(cmap=cm.bwr,value_list=final_sample_sub_clinical_df["HGB(g/L)"].values,vmin=80,vmax=160)  
#    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
#         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
#    pos= bottom+0.5*hight
#    bottom =bottom+(hight+hight_interval)
# #   pos= bottom+0.5*hight
#    y_position_list.append(pos)
    
    #size
    pc1_m,pc1_color_list = get_color_list(cmap=cm.YlGn,value_list=final_sample_sub_clinical_df["Final Size(cm)"].values,vmin=1,vmax=8)  
    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
  #  pos= bottom+0.5*hight
    y_position_list.append(pos)
    
    #purity_absolute
    pc1_m,pc1_color_list = get_color_list(cmap=cm.OrRd,value_list=final_sample_sub_clinical_df["purity_absolute"].values,vmin=0,vmax=1)  
    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
  #  pos= bottom+0.5*hight
    y_position_list.append(pos) 
    
    #StromalScore_estimate
#    pc1_m,pc1_color_list = get_color_list(cmap=cm.bwr,value_list=final_sample_sub_clinical_df["StromalScore_estimate"].values,vmin=-2000,vmax=1000)  
#    plt.bar(range(0,len(sample_list)),[hight]*len(sample_list),width=1,
#         bottom =[bottom]*len(sample_list),color = pc1_color_list,align='edge')
#    pos= bottom+0.5*hight
#    bottom =bottom+(hight+hight_interval)
#  #  pos= bottom+0.5*hight
#    y_position_list.append(pos) 
#    
    #
    
        
    
    ax.set_xlim([0,len(sample_list)])
    
#    yticks_list = ['Histology','Stage','Degree','HPV','Age','CA125(U/ml)',
#                                   'CA199(U/ml)','SCC(ng/ml)','HGB(g/L)','Size(cm)',"Absolute purity","Stromal score"]
    
    yticks_list = ['Histology','Stage','Degree','HPV','Size(cm)',"Absolute purity"]
    
    
    for i in range(0,len(yticks_list)):
        plt.text(-1,y_position_list[i],yticks_list[i],horizontalalignment="right",verticalalignment="center",size=8)
    

    
        

def other_square(ax):
    
    stage_colot_dict = {'I':'mistyrose', 
                        'II':'coral',
                        'III+IV':'brown'}
    
    cancer_type_color_list = {'Squamous':'tab:blue',
                        'Adenocarcinoma':'tab:orange',
                        'Adenosquamous':'tab:green',
                        'Small-Cell':'tab:red',
                        'Other':'tab:purple'}
        
    diff_degree_dict = {'G1':'g', 'G2':'b', 'G3':'r'}
    
    
    hpv_color_dict = {'A9':"b",'A7':'r',"A6":"c","A11":"m",'Negative':'grey'}
    
    mutation_type_dict={'Non-synonymous':'tab:blue','Synonymous':'tab:orange'}

    
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    
    i=0
    start=0
    interval = 0.9
    height=-4
    bar_width=0.1
    bottom =-1
    
    text_size=6
  #  pos=4
    for m_type in mutation_type_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=mutation_type_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom,m_type.replace(' ','\n'),size=text_size,verticalalignment="top")
        i+=1
    bottom = bottom+(-interval+height)
    
    
    i=0
    m=0
    interval=0.7
    for m_type in cancer_type_color_list.keys():
        if i >0 and i%2==0:
            bottom = bottom+(-interval+height)
            
            m=0
            plt.bar([m*interval+start],[height],bottom=bottom,width=bar_width,color=cancer_type_color_list[m_type],align='edge')
            plt.text(m*interval+start+bar_width,bottom,m_type.replace(' ','\n'),size=text_size,verticalalignment="top")
            m+=1
            
        else:
            plt.bar([m*interval+start],[height],bottom=bottom,width=bar_width,color=cancer_type_color_list[m_type],align='edge')
            plt.text(m*interval+start+bar_width,bottom,m_type.replace(' ','\n'),size=text_size,verticalalignment="top")
            m+=1
        i+=1
    bottom = bottom+(-interval+height)
           
    i=0
    interval=0.3
    for m_type in hpv_color_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=hpv_color_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom,m_type.replace(' ','\n'),size=text_size,verticalalignment="top")
        i+=1
    bottom = bottom+(-interval+height)
    
    i=0
    for m_type in stage_colot_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=stage_colot_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom,m_type.replace(' ','\n'),size=text_size,verticalalignment="top")
        i+=1
    bottom = bottom+(-interval+height)    
    
    i=0
    interval=0.3
    for m_type in diff_degree_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=diff_degree_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom,m_type.replace(' ','\n'),size=text_size,verticalalignment="top")
        i+=1    
 #   xlim=0
 #   xmax = (i-1)*interval+start
  #  plt.text(0,bottom+height-interval,'Well                 Poor',size=text_size,verticalalignment="top")
    
    bottom = bottom+(-interval+height)  
    
def drew_bar_size(ax):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = 1,vmax = 8)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.YlGn,norm=pc1_norm,ticks=[1,8],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("Size(cm)",size=6)
    cb1.ax.set_xticklabels(['1','8'],rotation=0,fontdict={'fontsize':5})


def drew_bar_tumor_purity_absolute(ax):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = 0,vmax = 1)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.OrRd,norm=pc1_norm,ticks=[0,1],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("Tumor purity",size=6)
    cb1.ax.set_xticklabels(['0%','100%'],rotation=0,fontdict={'fontsize':5})       


   
    
#tumor_maf_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\WES\all_sample_all_mutation_annovar_to_maf_add_ref_alt_count_remove_11_sample_coding.maf",sep="\t")    
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

gene_list =  ["PIK3CA","HLA-A","RHOA","ARID1A","FBXW7","LATS1","PTEN","FAT1","FAS","RB1","TP53","EP300","ERBB3","MAPK1","CASP8","KMT2C","KMT2D","KRAS"]

final_mutation_df_sort = rank_gene(gene_list)

mutation_ratio_dict = get_mutation_ratio(final_mutation_df_sort)
rank_mutation_df = get_gene_mutation_type(gene_list=gene_list)

rank_mutation_df.to_excel(r"fig1a_gene_mutation.xlsx")

fig = plt.figure(figsize=(12,8),dpi=300)
gs = GridSpec(30,30,figure=fig)

ax1 = fig.add_subplot(gs[0:6,0:20])
mutation_po_num_per_sample(ax1)
ax2 = fig.add_subplot(gs[6:13,0:20])

clinical_bar(ax2)

ax3 = fig.add_subplot(gs[13:26,0:20])

mutation_waterfall_plt(rank_mutation_df,ax3)  

ax4 = fig.add_subplot(gs[27:29,0:20])
mutation_type_square(ax4)  

ax5 = fig.add_subplot(gs[5:10,20:25])
other_square(ax5)

ax6 = fig.add_subplot(gs[11:12,20:22])
drew_bar_size(ax6)

ax6 = fig.add_subplot(gs[11:12,23:25])
drew_bar_tumor_purity_absolute(ax6)

plt.savefig("mutation_pattern_all_134_remove_TB113-TB83_tumor_sample.pdf",dpi=300,bbox_inches="tight")

    
    
    
    
    
    
    
    
    
    
    
    
    
    

        
            
        
    
    
    
    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



