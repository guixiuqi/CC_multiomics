# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 14:15:52 2022

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

import collections
import os

sample_info_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)

sample_type = []

for index in sample_info_df.index:
    if index.startswith("T"):
        sample_type.append("Tumor")
    else:
        sample_type.append("Normal")
        
sample_info_df["type"] = sample_type

squ_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["type"]=="Tumor")&(sample_info_df["Histological_Diagnosis"]=="Cervical Squamous Cell Carcinoma"),"Final clade"].values)

squ_hpv_type_df = pd.DataFrame.from_dict(squ_hpv_type_dict,orient ='index')

#squ_hpv_type_df = squ_hpv_type_df.drop(index="Unknown")

ratio_list = []
sum_=sum(squ_hpv_type_df[0])

for index in squ_hpv_type_df.index:
    num = squ_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
squ_hpv_type_df["ratio"] = ratio_list

squ_hpv_type_df.columns= ["squ_num","squ_ratio"]



ade_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["type"]=="Tumor")&(sample_info_df["Histological_Diagnosis"]=="Cervical Adenocarcinoma"),"Final clade"].values)

ade_hpv_type_df = pd.DataFrame.from_dict(ade_hpv_type_dict,orient ='index')
ade_hpv_type_df = ade_hpv_type_df.drop(index="Unknown")
ratio_list = []
sum_=sum(ade_hpv_type_df[0])

for index in ade_hpv_type_df.index:
    num = ade_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
ade_hpv_type_df["ratio"] = ratio_list

ade_hpv_type_df.columns= ["ade_num","ade_ratio"]



small_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["type"]=="Tumor")&(sample_info_df["Histological_Diagnosis"]=="Cervical Small Cell Carcinoma "),"Final clade"].values)

small_hpv_type_df = pd.DataFrame.from_dict(small_hpv_type_dict,orient ='index')
small_hpv_type_df = small_hpv_type_df
ratio_list = []
sum_=sum(small_hpv_type_df[0])

for index in small_hpv_type_df.index:
    num = small_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
small_hpv_type_df["ratio"] = ratio_list

small_hpv_type_df.columns= ["small_num","small_ratio"]







Adenosquamous_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["type"]=="Tumor")&(sample_info_df["Histological_Diagnosis"]=="Adenosquamous"),"Final clade"].values)

Adenosquamous_hpv_type_df = pd.DataFrame.from_dict(Adenosquamous_hpv_type_dict,orient ='index')
Adenosquamous_hpv_type_df = Adenosquamous_hpv_type_df
ratio_list = []
sum_=sum(Adenosquamous_hpv_type_df[0])

for index in Adenosquamous_hpv_type_df.index:
    num = Adenosquamous_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
Adenosquamous_hpv_type_df["ratio"] = ratio_list

Adenosquamous_hpv_type_df.columns= ["Adenosquamous_num","Adenosquamous_ratio"]




normal_squ_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["Histological_Diagnosis"]=="Cervical Squamous Cell Carcinoma")&(sample_info_df["type"]=="Normal"),"Final clade"].values)    
   
normal_squ_hpv_type_df = pd.DataFrame.from_dict(normal_squ_hpv_type_dict,orient ='index')
normal_squ_hpv_type_df = normal_squ_hpv_type_df.drop(index="Unknown")

ratio_list = []
sum_=sum(normal_squ_hpv_type_df[0])

for index in normal_squ_hpv_type_df.index:
    num = normal_squ_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
normal_squ_hpv_type_df["ratio"] = ratio_list

normal_squ_hpv_type_df.columns= ["normal_squ_num","normal_squ_ratio"]



normal_ade_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["Histological_Diagnosis"]=="Cervical Adenocarcinoma")&(sample_info_df["type"]=="Normal"),"Final clade"].values)    
   
normal_ade_hpv_type_df = pd.DataFrame.from_dict(normal_ade_hpv_type_dict,orient ='index')
#normal_ade_hpv_type_df = normal_ade_hpv_type_df.drop(index="Unknown")

ratio_list = []
sum_=sum(normal_ade_hpv_type_df[0])

for index in normal_ade_hpv_type_df.index:
    num = normal_ade_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
normal_ade_hpv_type_df["ratio"] = ratio_list

normal_ade_hpv_type_df.columns= ["normal_ade_num","normal_ade_ratio"]




normal_small_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["Histological_Diagnosis"]=="Cervical Small Cell Carcinoma ")&(sample_info_df["type"]=="Normal"),"Final clade"].values)    
   
normal_small_hpv_type_df = pd.DataFrame.from_dict(normal_small_hpv_type_dict,orient ='index')
normal_small_hpv_type_df = normal_small_hpv_type_df.drop(index="Unknown")

ratio_list = []
sum_=sum(normal_small_hpv_type_df[0])

for index in normal_small_hpv_type_df.index:
    num = normal_small_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
normal_small_hpv_type_df["ratio"] = ratio_list

normal_small_hpv_type_df.columns= ["normal_small_num","normal_small_ratio"]





normal_Adenosquamous_hpv_type_dict = collections.Counter(sample_info_df.loc[(sample_info_df["Histological_Diagnosis"]=="Adenosquamous")&(sample_info_df["type"]=="Normal"),"Final clade"].values)    
   
normal_Adenosquamous_hpv_type_df = pd.DataFrame.from_dict(normal_Adenosquamous_hpv_type_dict,orient ='index')
normal_Adenosquamous_hpv_type_df = normal_Adenosquamous_hpv_type_df

ratio_list = []
sum_=sum(normal_Adenosquamous_hpv_type_df[0])

for index in normal_Adenosquamous_hpv_type_df.index:
    num = normal_Adenosquamous_hpv_type_df.loc[index,0]
    ratio= num/sum_
    ratio_list.append(ratio)
normal_Adenosquamous_hpv_type_df["ratio"] = ratio_list

normal_Adenosquamous_hpv_type_df.columns= ["normal_Adenosquamous_num","normal_Adenosquamous_ratio"]







all_df = pd.concat([squ_hpv_type_df,ade_hpv_type_df,small_hpv_type_df,Adenosquamous_hpv_type_df,normal_squ_hpv_type_df,normal_ade_hpv_type_df,normal_small_hpv_type_df,normal_Adenosquamous_hpv_type_df],axis=1)
all_df = all_df.replace(np.nan,0)

hpv_color_dict = {'A9':"b",'A7':'r',"A6":"c","A11":"m",'Negative':'grey','A7;A9':'g'}

all_df = all_df.loc[['A9', 'A7', 'A7;A9', 'A6', 'A11','Negative']]
os.chdir(r"G:\cervical_cancer_multiomics\results\expression_table\\")
all_df.to_excel("hpv_clade_num_new.xlsx")



import numpy as np
import matplotlib.pyplot as plt


category_names = list(all_df.index)
ratio_results = {
    'Squamous': [i for i in list(all_df['squ_ratio'].values)],
    'Adenocarcinoma':[i for i in list(all_df['ade_ratio'].values)],
    'Adenosquamous':[i for i in list(all_df['Adenosquamous_ratio'].values)],
    'Small Cell':[i for i in list(all_df['small_ratio'].values)],    
    'Normal_Squamous':[i for i in  list(all_df['normal_squ_ratio'].values)],
    'Normal_Adenocarcinoma':[i for i in  list(all_df['normal_ade_ratio'].values)],
    'Normal_Adenosquamous':[i for i in  list(all_df['normal_Adenosquamous_ratio'].values)],
    'Normal_Small_Cell':[i for i in  list(all_df['normal_small_ratio'].values)]
    
}


results = {
    'Squamous': [i for i in list(all_df['squ_num'].values)],
    'Adenocarcinoma':[i for i in list(all_df['ade_num'].values)],
    'Adenosquamous':[i for i in list(all_df['Adenosquamous_num'].values)],
    'Small Cell':[i for i in list(all_df['small_num'].values)],    
    'Normal_Squamous':[i for i in  list(all_df['normal_squ_num'].values)],
    'Normal_Adenocarcinoma':[i for i in  list(all_df['normal_ade_num'].values)],
    'Normal_Adenosquamous':[i for i in  list(all_df['normal_Adenosquamous_num'].values)],
    'Normal_Small_Cell':[i for i in  list(all_df['normal_small_num'].values)]
    
}


def survey(results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(ratio_results.values()))
    data1 = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = [hpv_color_dict[i] for i in category_names]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        print(i)
        widths = data[:, i]
        widths1 = data1[:,i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        xcenters = starts + widths / 2

       # r, g, b, _ = color
        text_color = 'white'
        for y, (x, c) in enumerate(zip(xcenters, widths1)):
            
            if c>0:
                ax.text(x, y, str(int(c)), ha='center', va='center',
                    color=text_color,size=15)
    ax.legend(ncol=int(len(category_names)/2), bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='large')

    return fig, ax


survey(results, category_names)

os.chdir(r"G:\cervical_cancer_multiomics\results\expression_table\\")
plt.savefig(r"HPV_barh.pdf",dpi=300,bbox_inches="tight")






























