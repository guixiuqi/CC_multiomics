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

from scipy import stats
from statsmodels.stats.multitest import multipletests




def get_color_list(cmap=cm.bwr,value_list=[],vmin=-3,vmax=8):
    pc1_norm = mpl.colors.Normalize(vmin = vmin,vmax = vmax)
    pc1_cmap = cmap
    pc1_m=cm.ScalarMappable(norm=pc1_norm,cmap=pc1_cmap)
    pc1_color_list = []
    for pc in value_list:
        pc1_color_list.append(pc1_m.to_rgba(pc))
    return pc1_m,pc1_color_list




        
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
        table_list = []
        for color in color_list:
            table_list.append(sub_color_dict[color])
        all_table_list.append(table_list)

#        i += 1
    return all_table_list


def chi_pvalue(add_feature_cluster_df,cell_list):
    
    pvalue_list = []
    
    for column in cell_list:
        
        
        
        a='%s ~ C(subgroup)'%(column)
        model = ols(a, data=add_feature_cluster_df).fit()
        anova_table = sm.stats.anova_lm(model,typ=2)
        pc1_anova_pvalue = anova_table.loc['C(subgroup)','PR(>F)']
        
        pvalue_list.append(pc1_anova_pvalue)


    return pvalue_list
        


def cluster_map_sample_type(add_feature_cluster_df):  
    
    sample_sub_clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",
                                           index_col=0)       
    final_sample_sub_clinical_df = sample_sub_clinical_df.loc[tmt_cluster_df.index]
    stage_colot_dict = {'I':'mistyrose', 
                        'II':'coral',
                        'III':'brown',
                        'IV':'brown'}
    
    cancer_type_color_list = {'Cervical Squamous Cell Carcinoma':'tab:blue',
                        'Cervical Adenocarcinoma':'tab:orange',
                        'Adenosquamous':'tab:green',
                        'Cervical Small Cell Carcinoma ':'tab:red',
                        'Other':'tab:purple'}
        
    diff_degree_dict = {'G1':'g', 'G2':'b', 'G3':'r',np.nan:"white"}
    
    hpv_color_dict = {'A9':"b",'A7':'r',"A7;A9":"g","A6":"c","A11":"m",'Negative':'grey','Unknown':'white',np.nan:"white"}
    
#    y_position_list = []

    cell_list = []
    
    for column in sample_sub_clinical_df.columns:
        if "xcell" in column:
            cell_list.append(column)
            
            add_feature_cluster_df[column] = final_sample_sub_clinical_df[column]

        
        
#    cl_value_dict = {'CA125(U/ml)':35,'CA199(U/ml)':27,'SCC(ng/ml)':1.5,'HGB(g/L)':115}   
    pvalue_list = chi_pvalue(add_feature_cluster_df,cell_list)

    df_ = pd.DataFrame([cell_list,pvalue_list]).T
    df_.columns = ["cell_type","ANOVA_pvlaue"]
    corrected_p_values = multipletests(df_["ANOVA_pvlaue"], method='fdr_bh')[1]
    df_["BH p_value"] = corrected_p_values
    
    df_.to_excel(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\cell_type_xcell\tmt_cluster_x_cell_anova_bh_p.xlsx")
    return df_
        
 






if __name__=="__main__":
    os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance")
    
    tmt_hyper_z_df = pd.read_csv("top_variance_all_z_df_dropna_3000_1674.csv",sep=",",index_col=0)
    
    tmt_cluster_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\cluster_3.csv",
                                 sep=',',index_col=0)
#    tmt_cluster_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\consensusclusterplus\euclidean_km\cluster_5.csv",
#                                 sep=',',index_col=0)
    
    os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance")
    
    tmt_cluster_df = tmt_cluster_df.sort_values(by='x')
    collections.Counter(tmt_cluster_df['x'])
    #{1: 21, 2: 52, 3: 63}
    clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)
    
    tmt_cluster_df['type'] = clinical_df.loc[tmt_cluster_df.index,'Histological_Diagnosis']
    tmt_cluster_df = tmt_cluster_df.sort_values(by=['x','type'])
    
    final_sample = []

    for index in tmt_cluster_df.index:
        rna= clinical_df.loc[index,'StromaScore_xcell']
        if np.isnan(rna):
            pass
        else:
            final_sample.append(index)
            
    tmt_cluster_df = tmt_cluster_df.loc[final_sample]        
            
   
    tmt_hyper_z_cluster_df = tmt_hyper_z_df[tmt_cluster_df.index]
    
    tmt_z_hyper_normal_sample_df = tmt_hyper_z_df.drop(columns=tmt_cluster_df.index)
    
    
    cluster_color_dict = {1:"tab:green",2:"tab:blue",3:"tab:red"}
    
    sample_cluster_color_list = []
    for index in tmt_hyper_z_cluster_df.columns:
        color = cluster_color_dict[tmt_cluster_df.loc[index,"x"]]
        sample_cluster_color_list.append(color)
        
    add_feature_cluster_df = tmt_cluster_df.copy()
    add_feature_cluster_df.columns=["subgroup","type"]
    

    os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\cell_type_xcell")

    
    cluster_map_sample_type(add_feature_cluster_df)
    
    




    
    
    
    
    
    
    




