# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 14:00:10 2022

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
#import scipy
import matplotlib
import argparse
import os
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
from statsmodels.stats.multitest import multipletests

def get_color_list(cmap=cm.bwr,value_list=[],vmin=-3,vmax=8):
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    pc1_cmap = cmap
    pc1_m=cm.ScalarMappable(norm=pc1_norm,cmap=pc1_cmap)
    pc1_color_list = []
    for pc in value_list:
        pc1_color_list.append(pc1_m.to_rgba(pc))
    return pc1_m,pc1_color_list



# def marker_gene_heatmap(cluster_2d,ax,ax2,vmin_=-2,vmax_=2):
#     gene_list = []
#     cluster_gene_number_dict = {}
#     for key in marker_gene_dict.keys():
#         values = marker_gene_dict[key]
#         n=0
#         for v in values:
#             if v in normalized_df.index:
#                 gene_list.append(v)
#                 n+=1
#             else:
#                 pass
#         cluster_gene_number_dict[key] = n
    
#     df = pd.DataFrame(index=gene_list,columns=cluster_2d.columns)
        
#     for column in cluster_2d.columns:
#         if column in normalized_df.columns:
#             df[column] = normalized_df.loc[gene_list,column].values
#         else:
#             pass
#     df_normalized = df.sub(df.median(axis=1), axis=0)    
#     tumor_normalized_df = df_normalized.replace(np.nan,0)
    
#     normal_sample_list =[]
#     for column in normalized_df:
#         if column.startswith("N"):
#             normal_sample_list.append(column)
#         else:
#             pass
#     normal_df = normalized_df.loc[gene_list,normal_sample_list]
    
#     normal_normalized_df = normal_df.sub(df.median(axis=1), axis=0)  
        
    
#     g = sns.heatmap(tumor_normalized_df,vmin=vmin_,vmax=vmax_,cmap='bwr',cbar=False,xticklabels=False,yticklabels=False,ax=ax)
#     #plt.yticks(fontsize=4)
#     g = sns.heatmap(normal_normalized_df,vmin=vmin_,vmax=vmax_,cmap='bwr',cbar=False,xticklabels=False,yticklabels=True,ax=ax2)
#     plt.yticks(fontsize=6)
    
    
#     return tumor_normalized_df,normal_normalized_df
        

def data_for_feature(df):

    subgroup_list = set(df.iloc[:,1].values)
    color_list = list(set(df.iloc[:,0].values))
 #   i=0
 #   table_list = [[]]*len(subgroup_list)
    all_table_list = []
    for subgroup in subgroup_list:
        new_df = df.loc[df.iloc[:,1]==subgroup]
        sub_color_list = list(new_df.iloc[:,0].values)
   
        sub_color_dict = collections.Counter(sub_color_list)
         
        print(sub_color_dict)
        table_list = []
        for color in color_list:
            table_list.append(sub_color_dict[color])
        all_table_list.append(table_list)

#        i += 1
    return all_table_list


def chi_pvalue(add_feature_cluster_df,groups_):
    
    pvalue_list = []
    for column in ['Histological_Diagnosis', 'Stage', 'Grade', 'HPV clade','Lymph Node Status(0.Negative; 1.Positive)']:
  #   column = "TNM stage"
        if column !="HPV clade":
                
            df = add_feature_cluster_df[[column,"%s"%groups_]]
 #    df = df
            a = data_for_feature(df)
            g, pvalue, dof, expctd = chi2_contingency(a, lambda_="log-likelihood")
            pvalue_list.append(pvalue)
            
        else:
            
            df = add_feature_cluster_df[[column,"%s"%groups_]]
            
            df = df.loc[(df[column]=="A7")|(df[column]=="A9")]
            
            a = data_for_feature(df)
            g, pvalue, dof, expctd = chi2_contingency(a, lambda_="log-likelihood")
            pvalue_list.append(pvalue)
            

    for column in ['Tumor_size',
       'Age', 'StromalScore_estimate', 'ImmuneScore', ]:
        
        
        
        a='%s ~ C(%s)'%(column,groups_)
        model = ols(a, data=add_feature_cluster_df).fit()
        anova_table = sm.stats.anova_lm(model,typ=2)
        pc1_anova_pvalue = anova_table.loc['C(%s)'%(groups_),'PR(>F)']
        
        pvalue_list.append(pc1_anova_pvalue)
        
      
        
    return pvalue_list
        
    

def cluster_map_sample_type(cluster_2d,ax,add_feature_cluster_df,groups_):  
    
    sample_sub_clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",
                                           index_col=0)       
    final_sample_sub_clinical_df = sample_sub_clinical_df.loc[cluster_2d.columns]
    
    #final_sample_sub_clinical_df[]
    stage_colot_dict = {'I':'mistyrose', 
                        'II':'coral',
                        'III':'brown',
                        'IV':'brown'}
    
    cancer_type_color_list = {'Cervical Squamous Cell Carcinoma':'tab:blue',
                        'Cervical Adenocarcinoma':'tab:orange',
                        'Adenosquamous':'tab:green',
                        'Cervical Small Cell Carcinoma ':'tab:red',
                        'Other':'tab:purple'}
        
    diff_degree_dict = {'G1':'b', 'G2':'b', 'G3':'r',np.nan:"white"}
    
    lymph_node_status_dict = {0:"grey",1:"black"}
    
    hpv_color_dict = {'A9':"b",'A7':'r',"A7;A9":"g","A6":"c","A11":"m",'Negative':'grey','Unknown':'white',np.nan:"white"}
    
    tmt_color_dict = {1:"tab:green",2:"tab:blue",3:"tab:red",4:"tab:orange"}
    
    y_position_list = []
#    fig= plt.figure(figsize=(5,1))
#    ax=fig.gca()
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
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,align="edge",
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = sample_cluster_color_list)
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)

    y_position_list.append(pos)   
    
    
    
    # tmt_list = list(add_feature_cluster_df_sort['subgroup'].values)
    # plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,align="edge",
    #      bottom =[bottom]*len(final_sample_sub_clinical_df),color = [tmt_color_dict[i] for i in tmt_list])
    # pos= bottom+0.5*hight
    # bottom =bottom+(hight+hight_interval)

    # y_position_list.append(pos)
    
    
    
    
    
    
    add_feature_cluster_df['Histological_Diagnosis']=final_sample_sub_clinical_df['Histological_Diagnosis']
    
    i=0
    for histology_type in final_sample_sub_clinical_df['Histological_Diagnosis'].values:
       # histology_type = str(int(histology_type))
        cancer_color = [cancer_type_color_list[histology_type]]
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
    

    stage_combin_dict = {'I':'I', 
                        'II':'II',
                        'III':'III/IV',
                        'IV':'III/IV'}
    
    stage_list = list(final_sample_sub_clinical_df['Final stage'].values)
    add_feature_cluster_df['Stage'] = [stage_combin_dict[i] for i in stage_list]
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,align="edge",
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = [stage_colot_dict[i] for i in stage_list])
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)

    y_position_list.append(pos)
    
    #differentiation degree
    
    stage_list = list(final_sample_sub_clinical_df['Grade'])
    new_list = []
    
    for i in stage_list:
        if i=="G1":
            new_list.append("G1/2")
        elif i=="G2":
            new_list.append("G1/2")
        else:
            new_list.append("G3")
    
    
    add_feature_cluster_df['Grade']= new_list
    
    
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,align="edge",
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = [diff_degree_dict[i] for i in stage_list])
    
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
#    pos= bottom+0.5*hight
    y_position_list.append(pos)
    
    #HPV
    hpv_all_sample = list(final_sample_sub_clinical_df['Final clade'].values)
    
    add_feature_cluster_df["HPV clade"] = hpv_all_sample
    
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,align="edge",
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = [hpv_color_dict[i] for i in hpv_all_sample])
    
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
#    pos= bottom+0.5*hight
    y_position_list.append(pos)
    
    cl_value_dict = {'Lymph Node Status(0.Negative; 1.Positive)':0}
    
    for c in cl_value_dict.keys():
        if c=='HGB(g/L)' :
            color_list = []
            for sample in cluster_2d.columns:
                value_ = final_sample_sub_clinical_df.loc[sample,c]
                if np.isnan(value_):
                    color_list.append("white") 
                elif value_ < cl_value_dict[c]:
                    color_list.append("black")
                else:
                    color_list.append("grey")
            plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,
            bottom =[bottom]*len(final_sample_sub_clinical_df),color = color_list,align='edge')
    
            pos= bottom+0.5*hight
            bottom =bottom+(hight+hight_interval)
            y_position_list.append(pos)
            add_feature_cluster_df[c] = color_list
            
        else:
            color_list = []
            for sample in cluster_2d.columns:
                value_ = final_sample_sub_clinical_df.loc[sample,c]
                if np.isnan(value_):
                    color_list.append("white")
                elif value_ > cl_value_dict[c]:
                    color_list.append("black")
                else:
                    color_list.append("grey")
            plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,
            bottom =[bottom]*len(final_sample_sub_clinical_df),color = color_list,align='edge')
            add_feature_cluster_df[c] = color_list
    
            pos= bottom+0.5*hight
            bottom =bottom+(hight+hight_interval)
            y_position_list.append(pos)
            
                    
    #size
    
    
    pc1_m,pc1_color_list = get_color_list(cmap=cm.YlGn,value_list=final_sample_sub_clinical_df["Final Size(cm)"].values,vmin=1,vmax=6)  
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = pc1_color_list,align='edge')
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
  #  pos= bottom+0.5*hight
    y_position_list.append(pos)
  
    #age
    pc1_m,pc1_color_list = get_color_list(cmap=cm.OrRd,value_list=final_sample_sub_clinical_df["Age"].values,vmin=30,vmax=70)  
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = pc1_color_list,align='edge')
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
  #  pos= bottom+0.5*hight
    y_position_list.append(pos)


      
    pc1_m,pc1_color_list = get_color_list(cmap=cm.OrRd,value_list=final_sample_sub_clinical_df["StromalScore_estimate"].values,vmin=-2000,vmax=2000)  
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = pc1_color_list,align='edge')
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
  #  pos= bottom+0.5*hight
    y_position_list.append(pos)
    
    ###ImmuneScore
    pc1_m,pc1_color_list = get_color_list(cmap=cm.OrRd,value_list=final_sample_sub_clinical_df["ImmuneScore"].values,vmin=-1000,vmax=3000)  
    plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,
         bottom =[bottom]*len(final_sample_sub_clinical_df),color = pc1_color_list,align='edge')
    pos= bottom+0.5*hight
    bottom =bottom+(hight+hight_interval)
  #  pos= bottom+0.5*hight
    y_position_list.append(pos)    
    
    
  #   pc1_m,pc1_color_list = get_color_list(cmap=cm.OrRd,value_list=final_sample_sub_clinical_df["purity_absolute"].values,vmin=0,vmax=1)  
  #   plt.bar(range(0,len(final_sample_sub_clinical_df)),[hight]*len(final_sample_sub_clinical_df),width=1,
  #        bottom =[bottom]*len(final_sample_sub_clinical_df),color = pc1_color_list,align='edge')
  #   pos= bottom+0.5*hight
  #   bottom =bottom+(hight+hight_interval)
  # #  pos= bottom+0.5*hight
  #   y_position_list.append(pos)
    

    
    add_feature_cluster_df["Tumor_size"] = final_sample_sub_clinical_df["Final Size(cm)"]
    add_feature_cluster_df['Age']=final_sample_sub_clinical_df["Age"]
    add_feature_cluster_df['StromalScore_estimate']=final_sample_sub_clinical_df["StromalScore_estimate"]
    add_feature_cluster_df['ImmuneScore']=final_sample_sub_clinical_df["ImmuneScore"]
   # add_feature_cluster_df['purity_absolute']=final_sample_sub_clinical_df["purity_absolute"]
    
    
    
    ax.set_xlim([0,len(final_sample_sub_clinical_df)])
    

    
#    cl_value_dict = {'CA125(U/ml)':35,'CA199(U/ml)':27,'SCC(ng/ml)':1.5,'HGB(g/L)':115}   
    pvalue_list = chi_pvalue(add_feature_cluster_df,groups_)
    
    corrected_p_values = multipletests(pvalue_list, method='fdr_bh')[1]
    
    yticks_list = ['HPV_Group','Histology','Stage','Degree','HPV','Lymph node','Size(cm)','Age',"StromalScore","ImmuneScore"]#+gene_list
    
    HPV_cluster_chi_df = pd.DataFrame([yticks_list[1:],pvalue_list,corrected_p_values]).T
    HPV_cluster_chi_df.columns = ["Feature","chi-p value","BH-Pvalue"]
    HPV_cluster_chi_df.to_excel(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\hpv_cluster_chi_bh_pvalue.xlsx")
    
    for i in range(0,len(yticks_list)):
        plt.text(-1,y_position_list[i],yticks_list[i],horizontalalignment="right",verticalalignment="center",size=8)    
        if i < len(yticks_list)-1:
        
            plt.text(len(final_sample_sub_clinical_df),y_position_list[i+1],"%.2e"%(corrected_p_values[i]),horizontalalignment="left",verticalalignment="center",size=8) 




def pcc_heatmap(cluster_2d,vmin_=-0.5,vmax_=1):
    
    sns.heatmap(cluster_2d,vmin=vmin_,vmax=vmax_,cmap='bwr',cbar=False,yticklabels=True,xticklabels=False)

def other_square(ax):
    
    # group_dict = {"Group1":"tab:green","Group2":"tab:blue","Group3":"tab:red"}
    
    stage_colot_dict = {'I':'mistyrose', 
                        'II':'coral',
                        'III+IV':'brown'}
    
    cancer_type_color_list = {'Squamous':'tab:blue',
                        'Adenocarcinoma':'tab:orange',
                        'Adenosquamous':'tab:green',
                        'Small cell carcinoma':'tab:red',
                        'Other':'tab:purple'}
        
    diff_degree_dict = {'Grade1/2':'b', 'Grade3':'r'}
    
    hpv_color_dict = {'A9':"b",'A7':'r',"A7;A9":"g","A6":"c","A11":"m",'Negative':'grey'}
    
    cutoff_dict = {'Yes':"black",'No':'grey'}
    
  #  sample_type_dict = {'Tumor':'lightcoral','Normal':'lightseagreen'}
        
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    
    i=0
    start=0
    interval = 1
    height=-4
    bar_width=0.1
    bottom =-1
    
    text_size=8
    
    # i=0
    # for m_type in group_dict.keys():
    #     plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=group_dict[m_type],align='edge')
    #     plt.text(i*interval+start+bar_width,bottom+height/2,m_type,size=text_size,verticalalignment="center")
    #     bottom = bottom+(-interval+height)
        
    # bottom = bottom+(-interval+height)
    
    i=0
    m=0
    
    interval = 1
    for m_type in cancer_type_color_list.keys():
        if i >=0:
            
            m=0
            plt.bar([m*interval+start],[height],bottom=bottom,width=bar_width,color=cancer_type_color_list[m_type],align='edge')
            plt.text(m*interval+start+bar_width,bottom+height/2,m_type,size=text_size,verticalalignment="center")
            bottom = bottom+(-interval+height)
        else:
            plt.bar([m*interval+start],[height],bottom=bottom,width=bar_width,color=cancer_type_color_list[m_type],align='edge')
            plt.text(m*interval+start+bar_width,bottom+height/2,m_type,size=text_size,verticalalignment="center")
            bottom = bottom+(-interval+height)

    bottom = bottom+(-interval+height)
    
    interval=1
#    i=0
#    for m_type in sample_type_dict.keys():
#        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=sample_type_dict[m_type],align='edge')
#        plt.text(i*interval+start+bar_width,bottom+height/2,m_type.replace(' ',' '),size=text_size,verticalalignment="center")
#        bottom = bottom+(-interval+height)
#       # i+=1
#    bottom = bottom+(-interval+height)
           
    i=0
    for m_type in hpv_color_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=hpv_color_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom+height/2,m_type,size=text_size,verticalalignment="center")
        bottom = bottom+(-interval+height)
     #   i+=1
    bottom = bottom+(-interval+height)
    
    i=0
    for m_type in stage_colot_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=stage_colot_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom+height/2,m_type.replace(' ','\n'),size=text_size,verticalalignment="center")
        bottom = bottom+(-interval+height)
      #  i+=1
    bottom = bottom+(-interval+height)    
    
    
    i=0
    for m_type in cutoff_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=cutoff_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom+height/2,m_type.replace(' ','\n'),size=text_size,verticalalignment="center")
        bottom = bottom+(-interval+height)
       # i+=1
    bottom = bottom+(-interval+height)   
    
    
    
    
    i=0
   # interval=0.2
    for m_type in diff_degree_dict.keys():
        plt.bar([i*interval+start],[height],bottom=bottom,width=bar_width,color=diff_degree_dict[m_type],align='edge')
        plt.text(i*interval+start+bar_width,bottom+height/2,m_type.replace(' ','\n'),size=text_size,verticalalignment="center")
        bottom = bottom+(-interval+height)
      #  i+=1
    
 #   xlim=0
 #   xmax = (i-1)*interval+start
#    plt.text(0,bottom+height-interval,'Well                 Poor',size=text_size,verticalalignment="center")
    
    bottom = bottom+(-interval+height)  


def drew_bar_pcc(ax):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = -0.5,vmax = 1)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.bwr,norm=pc1_norm,ticks=[-0.5,1],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("PCC",size=6)
    cb1.ax.set_xticklabels(['-0.5','1'],rotation=0,fontdict={'fontsize':5})

def drew_bar_size(ax):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = 1,vmax = 6)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.YlGn,norm=pc1_norm,ticks=[1,6],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("Size(cm)",size=6)
    cb1.ax.set_xticklabels(['1','6'],rotation=0,fontdict={'fontsize':5})
    
    
def drew_bar_age(ax):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = 30,vmax = 70)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.OrRd,norm=pc1_norm,ticks=[30,70],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("Age(Year)",size=6)
    cb1.ax.set_xticklabels(['30','70'],rotation=0,fontdict={'fontsize':5})
    
# def drew_bar_tumor_purity_absolute(ax):
#     #gs = GridSpec(30,10,figure=fig)
#     #ax = fig.add_subplot(gs[26:27,7:9])
#     pc1_norm = mpl.colors.Normalize(vmin = 0,vmax = 1)
#     cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.OrRd,norm=pc1_norm,ticks=[0,1],orientation='horizontal')
#     #cb1.set_label('PC1',size=15)
#     ax.set_title("Tumor purity",size=6)
#     cb1.ax.set_xticklabels(['0%','100%'],rotation=0,fontdict={'fontsize':5})    
    
    
def drew_bar_tumor_immnuescore(ax):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = -1000,vmax = 3000)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.OrRd,norm=pc1_norm,ticks=[-1000,3000],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("ImmuneScore",size=6)
    cb1.ax.set_xticklabels(['-1000','3000'],rotation=0,fontdict={'fontsize':5})        
    

    
def drew_bar_tumor_stromalscore(ax):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = -2000,vmax = 2000)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.OrRd,norm=pc1_norm,ticks=[-2000,2000],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("StromalScore",size=6)
    cb1.ax.set_xticklabels(['-2000','2000'],rotation=0,fontdict={'fontsize':5})      
    
    


def drew_bar_z(ax,vmin = -3,vmax = 3):
    #gs = GridSpec(30,10,figure=fig)
    #ax = fig.add_subplot(gs[26:27,7:9])
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    cb1=mpl.colorbar.ColorbarBase(ax,cmap=cm.bwr,norm=pc1_norm,ticks=[vmin,vmax],orientation='horizontal')
    #cb1.set_label('PC1',size=15)
    ax.set_title("z-statistic",size=6)
    cb1.ax.set_xticklabels([str(vmin),str(vmax)],rotation=0,fontdict={'fontsize':5})

def pca_analysis(pearsonr_cor):
    from sklearn.decomposition import PCA
    X_new =  pearsonr_cor.T
    pca = PCA(n_components=2)
    X_r = pca.fit(X_new).transform(X_new)
    print('explained variance ratio (first three components): %s'
      % str(pca.explained_variance_ratio_))
    pc_df = pd.DataFrame(X_r,index = X_new.index,columns = ["PC1","PC2"])
    return pca,pc_df


        




if __name__=="__main__":
    

    
    tmt_cluster_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\cluster_3.csv",
                                 sep=',',index_col=0)
    
    tmt_cluster_df = tmt_cluster_df.sort_values(by='x')
    collections.Counter(tmt_cluster_df['x'])
    #{1: 21, 2: 52, 3: 63}
    clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)
    
    tmt_cluster_df['type'] = clinical_df.loc[tmt_cluster_df.index,'Histological_Diagnosis']
    tmt_cluster_df = tmt_cluster_df.sort_values(by=['x','type'])

    
    os.chdir(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression")

    
    dia_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\dia_z_statistic_table_removed_11_error_sample.xls",
                             sep="\t",index_col=0)
    
    tmt_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\tmt_z_statistic_table_remove_error_sample.csv",sep="\t",index_col=0)
    
    phos_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\tmt_phos_half_quantified_dreamai_imputation_res.csv",sep=",",index_col=0)
    
    
    hpv_gene_expression_tpm_raw = pd.read_csv("hpv_gene_expression_tpm.csv",sep="\t",index_col=0)

    hpv_gene_expression_tpm = np.log2(hpv_gene_expression_tpm_raw+1)
    hpv_gene_expression_tpm = hpv_gene_expression_tpm.drop(index=["E10"])
    
    
    hpv_add_clade = hpv_gene_expression_tpm.copy().T
    
    hpv_add_clade["sample_type"] = clinical_df.loc[hpv_add_clade.index,"sample type"]
    hpv_add_clade["hpv clade"] = clinical_df.loc[hpv_add_clade.index,"Final clade"]
    
    hpv_add_clade = hpv_add_clade.loc[(hpv_add_clade["hpv clade"]=="A7")|(hpv_add_clade["hpv clade"]=="A9")]
    
    
    #HPV gene expression accross gene type
#    hpv_color_dict = {'A9':"b",'A7':'r'}
#    for gene in hpv_add_clade.columns[0:-2]:
#        plt.figure(figsize=(3,5))
#        sns.violinplot(x="hpv clade",y=gene,data = hpv_add_clade,palette=hpv_color_dict)
#        sns.stripplot(x="hpv clade",y=gene,data = hpv_add_clade,color="black",order=["A9", "A7"])
#        
#        stats_,pvalue = ranksums(hpv_add_clade.loc[hpv_add_clade['hpv clade']=="A9",gene],hpv_add_clade.loc[hpv_add_clade['hpv clade']=="A7",gene])
#        
#        
#        plt.title("Wilcoxon rank-sum test\n%s"%("%.2e"%pvalue))
#        plt.xticks([0,1],["A9(95)","A7(35)"])
#        plt.ylabel("Log2(TPM) of %s"%gene)
#     #   plt.show()
#     #   plt.close("all")
#
#        plt.savefig(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\tpm_across_a7_a9\%s_tpm_across_a7_a9.png"%(gene),dpi=300,bbox_inches="tight")
#    
    ### scatterplot of HPV gene expression
    plt.figure(figsize=(7,7),dpi=300)
    plt.tight_layout()
    num=1
    index_list = hpv_gene_expression_tpm.index
    for gene in range(0,len(hpv_gene_expression_tpm.index)):
        for index in range(gene+1,len(hpv_gene_expression_tpm.index)):
            
            num= 7*(gene)+index
            plt.subplot(7,7,num)
            
            plt.scatter(hpv_gene_expression_tpm.iloc[index],hpv_gene_expression_tpm.iloc[gene],s=0.2)
            
            stats_,pvalue = spearmanr(hpv_gene_expression_tpm.iloc[index],hpv_gene_expression_tpm.iloc[gene])
            plt.text(5,2,"r=%s\np=%s"%("%.2f"%stats_,"%.2e"%(pvalue)),size=5)
            
            if index==gene+1:
                plt.xlabel(index_list[index],size=6)
                plt.ylabel(index_list[gene],size=6)
           # num+=1
                plt.xticks([0,5,10,15],size=6)
                plt.yticks([0,5,10,15],size=6)
                plt.plot([0,15],[0,15],linestyle="--",color="grey")
            else:
             #   plt.xlabel(index_list[index],size=6)
             #   plt.ylabel(index_list[gene],size=6)
           # num+=1
                plt.xticks([],size=6)
                plt.yticks([],size=6)
                plt.plot([0,15],[0,15],linestyle="--",color="grey")
                
#    plt.tight_layout()
    plt.savefig(r"gene_expression_scatterplot.pdf",dpi=300)
    
 
    ###densityplot of hpv expression
    plt.figure(figsize=(1,4),dpi=300)
    
    df =  hpv_gene_expression_tpm.T
    num=1
    for gene in hpv_gene_expression_tpm.index:
       
        plt.subplot(7,1,num)
      #  plt.hist(tmp_list,bins=50)
        #df.gene.plot.density(color='green')
        sns.kdeplot(data=df, x=gene)
        
        if num==7:
            plt.xlabel("Log$_2$(TPM)",size=6)
            plt.xticks([0,5,10,15],size=6)
            plt.ylabel(gene,size=6)
            plt.yticks([])
        else:
            plt.xlabel("")
            plt.xticks([])
            plt.yticks([])
            plt.ylabel(gene,size=6)
            
            
        
       # plt.yticks([0,25,50],size=6)
            
        num+=1
        
    plt.savefig("hpv_gene_expression_densityplot.pdf",bbox_inches="tight",dpi=300)
     
    
    ####HPV cluster
    
    hpv_gene_expression_tpm_raw = pd.read_csv("hpv_gene_expression_tpm.csv",sep="\t",index_col=0)

    hpv_gene_expression_tpm = np.log2(hpv_gene_expression_tpm_raw+1)
    hpv_gene_expression_tpm = hpv_gene_expression_tpm.drop(index=["E10"])
 
  #  hpv_gene_expression_tpm = hpv_gene_expression_tpm.drop(index=["E7","E1","L1"])
    
    hpv_gene_expression_tpm = hpv_gene_expression_tpm.drop(index=["E7","E1","E6","L1","L2"])
    
    tumor_sample_list = []
    normal_sample_list = []
    for column in hpv_gene_expression_tpm.columns:
        if column.startswith("T"):
            tumor_sample_list.append(column)
        else:
            normal_sample_list.append(column)
            
            
    final_tumor_sample = list(set(tumor_sample_list)&set(tmt_cluster_df.index))
    
    tmt_hyper_z_cluster_df = hpv_gene_expression_tpm[final_tumor_sample]
    
  #  g = sns.clustermap(tmt_hyper_z_cluster_df,col_cluster=True,row_cluster= False,cmap="bwr",method="ward",metric="correlation")
  
    g = sns.clustermap(tmt_hyper_z_cluster_df,col_cluster=True,row_cluster= False,cmap="bwr",method="ward")
    
    
    
    plt.savefig("hpv_expression_cluster.pdf",dpi=300,bbox_inches="tight")
    
 #   cluster_2d = g.data2d
    
    
    
    cluster_color_dict = {2:"#F8B62D",1:"#601986",3:"tab:red",4:"tab:orange"}

        
    add_feature_cluster_df = tmt_cluster_df.loc[final_tumor_sample]
    add_feature_cluster_df.columns=["subgroup","type"]
        
    for cluster_num in [2]:
        samples_indexs=scipy.cluster.hierarchy.fcluster(g.dendrogram_col.calculated_linkage,t=cluster_num,criterion='maxclust')
        
        sub_sample_clusters={}
        for sample ,cluster in zip(tmt_hyper_z_cluster_df,samples_indexs):
            if cluster in sub_sample_clusters.keys():
                sub_sample_clusters[cluster].append(sample)
            else:
                sub_sample_clusters[cluster] = [sample]
                
        add_feature_cluster_df["hpv_cluster_%s"%(cluster_num)] = [0]*len(add_feature_cluster_df)
        
        for index in range(1,cluster_num+1):
            samples = sub_sample_clusters[index]
            for sample in samples:
            
                add_feature_cluster_df.loc[sample,"hpv_cluster_%s"%(cluster_num)] = index
                
        add_feature_cluster_df_sort = add_feature_cluster_df.sort_values(by="hpv_cluster_%s"%(cluster_num))
            
    
        cluster_2d = g.data2d
        cluster_2d = cluster_2d[add_feature_cluster_df_sort.index]
        
        cluster_color_dict = {2:"#F8B62D",1:"#601986",3:"tab:red",4:"tab:orange"}
        
        
        sample_cluster_color_list = []
        for sample in cluster_2d.columns:   
            color = cluster_color_dict[add_feature_cluster_df_sort.loc[sample,"hpv_cluster_%s"%(cluster_num)]]
            sample_cluster_color_list.append(color)
            
        import argparse
        import matplotlib.gridspec as gridspec
        from matplotlib.gridspec import GridSpec
        
        fig = plt.figure(figsize=(10,10),dpi=300)
        gs = GridSpec(60,30,figure=fig)
        
        
        ax2 = fig.add_subplot(gs[0:30,8:25])
        
        cluster_map_sample_type(cluster_2d,ax=ax2,add_feature_cluster_df=add_feature_cluster_df_sort,groups_="hpv_cluster_%s"%(cluster_num))
        
        ax3 = fig.add_subplot(gs[30:40,8:25])    
        pcc_heatmap(cluster_2d,vmin_=0,vmax_=15)
        
        
        
        ax5 = fig.add_subplot(gs[0:30,0:1])
        other_square(ax5)
        
        
        x=0
        y=2
        ax6 = fig.add_subplot(gs[32:33,x:y])
        drew_bar_size(ax6)
        
        ax7 = fig.add_subplot(gs[36:37,x:y])
        drew_bar_age(ax7)
        
        ax9 = fig.add_subplot(gs[40:41,x:y])
        drew_bar_tumor_stromalscore(ax9)
        
        
        ax10 = fig.add_subplot(gs[44:45,x:y])
        drew_bar_tumor_immnuescore(ax10)
        
        
        # ax10 = fig.add_subplot(gs[48:49,x:y])
        # drew_bar_tumor_purity_absolute(ax10)  
        
    
        os.chdir(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression")
    
        plt.savefig("hpv_cluster_%s_heatmap_new_BH_adjusted_P.pdf"%(cluster_num),dpi=300,bbox_inches="tight")
        
        add_feature_cluster_df_sort.to_csv("hpv_gene_cluster_df.csv",sep="\t")



for i in ["Age"]:

    plt.figure(figsize=(1,2),dpi=300)
    palette = cluster_color_dict = {1:"tab:green",2:"tab:blue"}

#    sns.violinplot(x='hpv_cluster_2',y=i,data = add_feature_cluster_df_sort,palette=palette,fliersize=1,linewidth=0.5)  
    sns.boxplot(x='hpv_cluster_2',y=i,data = add_feature_cluster_df_sort,color="white",fliersize=1,linewidth=0.5)  
    sns.stripplot(x='hpv_cluster_2',y=i,data = add_feature_cluster_df_sort,palette=palette,s=2)
    plt.xticks([0,1],["HPV_C1","HPV_C2"],size=6)
    #plt.xlabel("")
    plt.ylim([20,100])
    
    plt.savefig(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\HPV_e2_e5_cluster_2\%s_across_hpv_cluster_box.pdf"%(i),
                dpi=300,bbox_inches="tight")
    


#immnuescore across HPV clade

hpv_h7_h9_df = add_feature_cluster_df_sort.loc[(add_feature_cluster_df_sort["HPV clade"]=="A7")|(add_feature_cluster_df_sort["HPV clade"]=="A9")]







    
    
    
add_feature_cluster_df["E2"] = cluster_2d.loc["E2",add_feature_cluster_df.index]

add_feature_cluster_df["E5"] = cluster_2d.loc["E5_ALPHA",add_feature_cluster_df.index]


scipy.stats.ttest_ind(add_feature_cluster_df.loc[add_feature_cluster_df["hpv_cluster_2"]==1,"E2"],
                      add_feature_cluster_df.loc[add_feature_cluster_df["hpv_cluster_2"]==2,"E2"]) #P = 8.95e-24

scipy.stats.ttest_ind(add_feature_cluster_df.loc[add_feature_cluster_df["hpv_cluster_2"]==1,"E5"],
                      add_feature_cluster_df.loc[add_feature_cluster_df["hpv_cluster_2"]==2,"E5"]) #P = 4.16e-52

    
    
    
    
    
    
    
    
    




