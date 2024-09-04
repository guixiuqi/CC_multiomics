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

os.chdir(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\HPV_e2_e5_cluster_2")

response_df = pd.read_excel(r"G:\余竟_ppt\宫颈癌随访数据更新至2023年7月24日-to gxq.xlsx",index_col=0)

response_df = response_df[['是否进展（0，否；1，是）（截止到2023年7月24日）','PFS（截止到2023年7月24日）']]
response_df.columns = ["Progress","PFS"]

clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)


tmt_cluster_info = pd.read_csv(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\hpv_gene_cluster_df.csv",
                               sep="\t",index_col=0)

tumor_response_df = response_df.loc[tmt_cluster_info.index]

tumor_response_df["hpc_cluster_2"] = tmt_cluster_info.loc[tumor_response_df.index,"hpv_cluster_2"]


chemo_radio_df =tumor_response_df.loc[(clinical_df['Treatment']!="Surgery")&(clinical_df["sample type"]!="Normal")&(clinical_df['Treatment']!='Other')]

chemo_radio_df = chemo_radio_df[['Progress', 'PFS','hpc_cluster_2']]
chemo_radio_df = chemo_radio_df.dropna(how="any")


chemo_radio_df.to_excel(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\HPV_gene_expression\HPV_e2_e5_cluster_2\hpv_cluster_survival.xlsx")



cluster_name = "hpc_cluster_2"
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
    


T1 = chemo_radio_df.loc[chemo_radio_df['hpc_cluster_2']==1,'PFS']#[1, 4, 10, 12, 12, 3, 5.4]
E1 = chemo_radio_df.loc[chemo_radio_df['hpc_cluster_2']==1,'Progress'] #[1, 0, 1,  0,  1,  1, 1]

T2 = chemo_radio_df.loc[chemo_radio_df['hpc_cluster_2']==2,'PFS']
E2 = chemo_radio_df.loc[chemo_radio_df['hpc_cluster_2']==2,'Progress']




from lifelines.statistics import logrank_test
results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
g1_g2_p = results.p_value


plt.rcParams['xtick.labelsize']=6
plt.rcParams['ytick.labelsize']=6
plt.ylim([0,1.1])
plt.text(5,0.5,"hpv_c1 vs hpv_c2 P-value: "+'{:.2e}'.format(g1_g2_p),size=6)
ax.set_xlabel("Progress free time",size=6)
ax.set_ylabel("Non-progress ratio",size=6)

plt.savefig("PFS_hpv_e2_e5_7月份随访结果.pdf",dpi=300,bbox_inches="tight")
