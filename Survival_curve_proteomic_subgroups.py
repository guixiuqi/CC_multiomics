# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 14:43:21 2022

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
import os
import matplotlib
import argparse
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec

#matplotlib.use('Agg')
matplotlib.rcParams['font.family'] ='Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import glob
import collections
import seaborn as sns

os.chdir(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km")



response_df = pd.read_excel(r"G:\余竟_ppt\宫颈癌随访数据更新至2023年7月24日-to gxq.xlsx",index_col=0)

response_df = response_df[['是否进展（0，否；1，是）（截止到2023年7月24日）','PFS（截止到2023年7月24日）']]
response_df.columns = ["Progress","PFS"]

response_df.to_excel("survival_time.xlsx")

response_df


clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)


tmt_cluster_info = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\cluster_3.csv",
                               sep=",",index_col=0)

tumor_response_df = response_df.loc[tmt_cluster_info.index]

tumor_response_df["tmt_cluster_top_valance"] = tmt_cluster_info.loc[tumor_response_df.index,"x"]


chemo_radio_df =tumor_response_df.loc[(clinical_df['Treatment']!="Surgery")&(clinical_df["sample type"]!="Normal")&(clinical_df['Treatment']!='Other')]

chemo_radio_df = chemo_radio_df[['Progress', 'PFS','tmt_cluster_top_valance']]
chemo_radio_df = chemo_radio_df.dropna(how="any")
chemo_radio_df.to_excel("tmt_radio_km_data.xlsx")

cluster_name = "tmt_cluster_top_valance"
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test
kmf = KaplanMeierFitter()
fig = plt.figure(figsize=(3,3),dpi=300)
color_dict={1:"green",2:"blue",3:"red"}
ax = fig.gca()
  #  survival_df[gene] = snp_df.loc[gene,survival_df.index]
for name, grouped_df in chemo_radio_df.groupby(cluster_name):
    kmf.fit(grouped_df['PFS'], grouped_df['Progress'], label="Group"+str(name)+"(%s)"%str(len(grouped_df)),
                alpha =0.1)

    kmf.plot(ax=ax,show_censors=True,ci_show=False,color=color_dict[name],
             censor_styles={'ms': 6},linewidth=1)
    
results = multivariate_logrank_test(chemo_radio_df['PFS'],
                                   chemo_radio_df[cluster_name], 
                                   chemo_radio_df['Progress'])


T1 = chemo_radio_df.loc[chemo_radio_df['tmt_cluster_top_valance']==1,'PFS']#[1, 4, 10, 12, 12, 3, 5.4]
E1 = chemo_radio_df.loc[chemo_radio_df['tmt_cluster_top_valance']==1,'Progress'] #[1, 0, 1,  0,  1,  1, 1]

T2 = chemo_radio_df.loc[chemo_radio_df['tmt_cluster_top_valance']==2,'PFS']
E2 = chemo_radio_df.loc[chemo_radio_df['tmt_cluster_top_valance']==2,'Progress']


T3 = chemo_radio_df.loc[chemo_radio_df['tmt_cluster_top_valance']==3,'PFS']
E3 = chemo_radio_df.loc[chemo_radio_df['tmt_cluster_top_valance']==3,'Progress']

from lifelines.statistics import logrank_test
results = logrank_test(T1, T3, event_observed_A=E1, event_observed_B=E3)
g1_g3_p = results.p_value

results = logrank_test(T2, T3, event_observed_A=E2, event_observed_B=E3)
g2_g3_p = results.p_value

plt.rcParams['xtick.labelsize']=6
plt.rcParams['ytick.labelsize']=6
plt.ylim([0,1.1])
plt.text(5,0.5,"Group1 vs Group3 P-value: "+'{:.2e}'.format(g1_g3_p)+"\n"+"Group2 vs Group3 P-value: "+'{:.2e}'.format(g2_g3_p),size=6)
ax.set_xlabel("Progress free time(Month)",size=6)
ax.set_ylabel("Non-progress ratio",size=6)

plt.savefig("PFS_TMT_top_variance_3_cluster_survival_analysis.pdf",dpi=300,bbox_inches="tight")
