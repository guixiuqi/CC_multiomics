# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 20:31:45 2022

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

import seaborn as sns
import os
os.chdir(r"G:\cervical_cancer_multiomics\results\expression_table\\")

dia_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\dia_z_statistic_table_removed_11_error_sample.xls",sep="\t",index_col=0)


tmt_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\tmt_z_statistic_table_remove_error_sample.csv",sep="\t",index_col=0)

tmt_p_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_p_protein\tmt_p_normalized_log2_statistic_table_remove_error_sample_non_uniq_index.csv",sep="\t",index_col=0)

rna_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\rna_seq_normalized_read_counts_final_sample_name_removed_11_error_sample.xls",
                     sep="\t",index_col=0)

tmt_ace_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_ace_protein\tmt_ace_normalized_log2_statistic_table_remove_error_sample_non_uniq_gene.csv",sep="\t",index_col=0)

pair_WES_df =  pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\all_pair.txt",sep="\t")
no_blood_WES_sample = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\final_without_blood_sample.txt",sep="\t",index_col=0,header=None)

def get_tumor_normal_sample(dia_df):
    
    tumor_sample = []
    normal_sample = []
    
    for column in dia_df.columns:
        if column.startswith("T") and '.' not in column:
            tumor_sample.append(column)
        elif column.startswith("N"):
            normal_sample.append(column)
    print(len(tumor_sample),len(normal_sample))
    return tumor_sample,normal_sample

def rename(index):
    sets = index.split("-")
    if len(sets) == 3:
        new_index = sets[0]+sets[1]+"-"+sets[2]
    elif len(sets)==2:
        if len(sets[0])==1:
            new_index = sets[0]+sets[1]
        else:
            new_index = sets[0]+"-"+sets[1]
    else:
        new_index = sets[0]
    if new_index.endswith("-2"):
        final_index="N"+new_index
    #rna_new_id_dict[index] = new_index
 #       final_index_list.append(final_index)
    else:
        final_index = "T"+ new_index
    return final_index

def WES_tumor_dict(pair_WES_df,no_blood_WES_sample):
    
    WES_tumor_list = []
    wes_tumor_blood_list = []
    WES_normal_list = []
    
    for index in pair_WES_df.index:
        new_name = rename(index)
        if new_name.startswith("T"):
            WES_tumor_list.append(new_name)
            wes_tumor_blood_list.append(new_name)
        else:
            WES_normal_list.append(new_name)
    for index in no_blood_WES_sample.index:
        new_name = rename(index)
        if new_name.startswith("T"):
            WES_tumor_list.append(new_name)
            #wes_tumor_blood_dict[new_name]=0
        else:
            WES_normal_list.append(new_name)
    print(len(WES_tumor_list),len(wes_tumor_blood_list),len(WES_normal_list))
    return WES_tumor_list,wes_tumor_blood_list,WES_normal_list
        
    


rna_tumor_sample,rna_normal_sample = get_tumor_normal_sample(rna_df)

tmt_tumor_sample,tmt_normal_sample = get_tumor_normal_sample(tmt_df)
dia_tumor_sample,dia_normal_sample = get_tumor_normal_sample(dia_df)

tmt_p_tumor_sample,tmt_p_normal_sample = get_tumor_normal_sample(tmt_p_df)

tmt_ace_tumor_sample,tmt_ace_normal_sample = get_tumor_normal_sample(tmt_ace_df)



wes_tumor_list,wes_tumor_with_blood,wes_normal_list =  WES_tumor_dict(pair_WES_df,no_blood_WES_sample)



all_tumor_sample = list(set(rna_tumor_sample)|set(tmt_tumor_sample)|set(dia_tumor_sample)|set(tmt_p_tumor_sample)|set(wes_tumor_list)|set(tmt_ace_tumor_sample))
all_normal_sample = list(set(rna_normal_sample)|set(tmt_normal_sample)|set(dia_normal_sample)|set(tmt_p_normal_sample)|set(wes_normal_list)|set(tmt_ace_normal_sample))


all_data_2d = []

for data in [wes_tumor_list,wes_tumor_with_blood,wes_normal_list,rna_tumor_sample,rna_normal_sample,tmt_tumor_sample,
             tmt_normal_sample,tmt_p_tumor_sample,tmt_p_normal_sample,tmt_ace_tumor_sample,tmt_ace_normal_sample,dia_tumor_sample,dia_normal_sample]:
    data_1d_list = []
    for pati in all_tumor_sample+["NB112-2","NB139-2","NB142-2"]:
        if pati in data or "N"+pati[1:-1]+"2" in data:
            data_1d_list.append(1)
        else:
            data_1d_list.append(0)
    all_data_2d.append(data_1d_list)
    
data_df = pd.DataFrame(all_data_2d,index=["Tumor tissue WES(%s)"%(len(wes_tumor_list)),
                                          "Blood WES(%s)"%(len(wes_tumor_with_blood)),
                                          "Adjacent non-tumor tissue WES(%s)"%(len(wes_normal_list)),
                                          "Tumor tissue RNA(%s)"%(len(rna_tumor_sample)),
                                          "Adjacent non-tumor tissue RNA(%s)"%(len(rna_normal_sample)),
                                          "Tumor tissue proteome(TMT:%s)"%(len(tmt_tumor_sample)),
                                          "Adjacent non-tumor tissue proteome(TMT:%s)"%(len(tmt_normal_sample)),                                         
                                          "Tumor tissue phosphoproteome(TMT:%s)"%(len(tmt_p_tumor_sample)),
                                          "Adjacent non-tumor tissue phosphoproteome(TMT:%s)"%(len(tmt_p_normal_sample)),
                                            "Tumor tissue acetylproteome(TMT:%s)"%(len(tmt_ace_tumor_sample)),
                                          "Adjacent non-tumor tissue acetylproteome(TMT:%s)"%(len(tmt_ace_normal_sample)),
                                           "Tumor tissue proteome(DIA:%s)"%(len(dia_tumor_sample)),
                                          "Adjacent non-tumor tissue proteome(DIA:%s)"%(len(dia_normal_sample)),],columns=all_tumor_sample+["NB112-2","NB139-2","NB142-2"])


clinical_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_add_tumor_purity_estimation_check_hpv_add_clade_follow_up_2022_8_31.csv",sep=",",index_col=0)

data_df_sort = data_df.sort_values(by=list(data_df.index),axis=1,ascending = False)

data_df_sort.to_excel("Figs1b.xlsx")

clinical_df["sample_type_bar"] = clinical_df["Histological_Diagnosis"]
for sample in ["NB112-2","NB139-2","NB142-2"]:
    clinical_df.loc[sample,"sample_type_bar"] = "Patients with only\n adjacent non-tumor tissue"


all_df_list = []
sample_type_list = ['Cervical Squamous Cell Carcinoma','Cervical Adenocarcinoma','Adenosquamous',
                    'Cervical Small Cell Carcinoma ','Other',"Patients with only\n adjacent non-tumor tissue"]
sample_type_len_dict = {}
for sample_type in sample_type_list:
    sample = list(clinical_df.loc[clinical_df["sample_type_bar"]==sample_type].index)
    common_sample = set(sample) & set(all_tumor_sample+["NB112-2","NB139-2","NB142-2"])
    sample_type_len_dict[sample_type]=len(common_sample)
    
    sample_df = data_df_sort[common_sample]
    sample_df = sample_df.sort_values(by=list(sample_df.index),axis=1,ascending = False)
    
    
    all_df_list.append(sample_df)


all_concat_df = pd.concat(all_df_list,axis=1)

all_concat_df.to_excel(r"G:\cervical_cancer_multiomics\results\expression_table\Figs1b.xlsx")

color_list = ["#F8766D", "#E88526", "#D39200", "#B79F00", "#93AA00", "#5EB300", "#00BA38", "#00BF74", "#00C19F",
"#00BFC4", "#00B9E3", "#00ADFA", "#619CFF", "#AE87FF", "#DB72FB", "#F564E3", "#FF61C3", "#FF699C"]

cancer_type_color_dict = {'Cervical Squamous Cell Carcinoma':'tab:blue',
                        'Cervical Adenocarcinoma':'tab:orange',
                        'Adenosquamous':'tab:green',
                        'Cervical Small Cell Carcinoma ':'tab:red',
                        'Other':'tab:purple',
                        "Patients with only\n adjacent non-tumor tissue":"grey"}

sample_cancer_type_color_list = []

for sample in all_concat_df.columns:
    sample_type = clinical_df.loc[sample,"sample_type_bar"]
    sample_cancer_type_color_list.append(cancer_type_color_dict[sample_type])
    




fig = plt.figure(figsize=(10,4))   

import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(10,4),dpi=300)
plt.rcParams['xtick.labelsize']=15
plt.rcParams['ytick.labelsize']=15

gs = GridSpec(10,20,figure=fig)

ax = fig.add_subplot(gs[0:10,0:15])

#ax=fig.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

y_loc_list = []
i=0
height = -1
interval=-0.1
bottom=0
sn = len(all_tumor_sample+["NB112-2","NB139-2","NB142-2"])

plt.bar([i for i in range(0,sn)],[0.5]*sn,bottom=[0.1]*sn,color=sample_cancer_type_color_list,align="edge")

for index in all_concat_df.index:
    values = all_concat_df.loc[index].values
    index_color_list = []
    for value in values:
        if value==1:
            index_color_list.append(color_list[i])
        else:
            index_color_list.append("whitesmoke")
    plt.bar([i for i in range(0,sn)],[height]*sn,bottom=[bottom]*sn,color=index_color_list,align="edge")
    
    y_loc_list.append(bottom+height/2)
    bottom=bottom+height+interval
    i+=1
plt.xlim([0,sn])
plt.xticks([])
plt.yticks(y_loc_list,list(all_concat_df.index))
ax.tick_params(axis='both', which='both', length=0)


cancer_type_color_dict = {'Squamous(112)':'tab:blue',
                        'Adenocarcinoma(19)':'tab:orange',
                        'Adenosquamous(2)':'tab:green',
                        'Small Cell Carcinoma(4)':'tab:red',
                        'Other(2)':'tab:purple',
                        "Patients with only\n adjacent non-tumor tissue(3)":"grey"}


ax1 = fig.add_subplot(gs[1:5,15:16])
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.set_xticks([])
ax1.set_yticks([])
i=0
m=0
i=0
start=0
interval = 1
height=-4
bar_width=0.1
bottom =-1

text_size=8
interval = 1
for m_type in cancer_type_color_dict.keys():
    if i >=0:
        
        m=0
        plt.bar([m*interval+start],[height],bottom=bottom,width=bar_width,color=cancer_type_color_dict[m_type],align='edge')
        plt.text(m*interval+start+bar_width,bottom+height/2,m_type,size=text_size,verticalalignment="center")
        bottom = bottom+(-interval+height)
    else:
        plt.bar([m*interval+start],[height],bottom=bottom,width=bar_width,color=cancer_type_color_dict[m_type],align='edge')
        plt.text(m*interval+start+bar_width,bottom+height/2,m_type,size=text_size,verticalalignment="center")
        bottom = bottom+(-interval+height)

plt.savefig(r"G:\cervical_cancer_multiomics\results\expression_table\all_data_no_remove_sample_bar.pdf",dpi=300,bbox_inches="tight")
        
    
plt.savefig(r"G:\cervical_cancer_multiomics\results\expression_table\all_data_no_remove_sample_bar.png",dpi=300,bbox_inches="tight")
        
    

















   


