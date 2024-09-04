# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 10:47:28 2023

@author: guixiuqi
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
from scipy.stats import spearmanr
from scipy.stats import pearsonr
#from matplotlib_venn import venn3
import os
#os.chdir(r"D:\python_module")
#import venn
#import mkdir_
#import get_rgb_color_list

import scipy.stats as stats
from scipy.stats import chi2_contingency
import scipy
import matplotlib
import argparse

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
from statannotations.Annotator import Annotator
import argparse
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import matplotlib.cm as cm
plt.rcParams['xtick.labelsize']=6
plt.rcParams['ytick.labelsize']=6
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'


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
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
import os
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.statistics import multivariate_logrank_test
import seaborn as sns

filedir = r"G:\cervical_cancer_multiomics\results\TMT_protein\feature_PFS_cox\\"
os.chdir(filedir)

#subgroup_df = pd.read_csv(filedir + "subgroup_by_hyper_across_all_sample_pc1_survival_df.csv",sep="\t",index_col=0)

clinical_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\expression_table\sample_clinical_info_estimate_xcell_absolute_2022_10_12.xlsx",index_col=0)

tmt_group_df = pd.read_excel(r"G:\cervical_cancer_multiomics\results\TMT_protein\top_variance\euclidean_km\multi_omics_cluster.xlsx",index_col=0)
tmt_group_df = tmt_group_df["TMT_proteomic_subgroup"].to_frame()
tmt_group_df = tmt_group_df.dropna(how="any")

PFS_df = pd.read_excel(r"G:\余竟_ppt\宫颈癌随访数据更新至2023年7月24日-to gxq.xlsx",index_col=0)


chemoradio_df = clinical_df.loc[(clinical_df['Treatment']!="Surgery")&(clinical_df["sample type"]!="Normal")&(clinical_df['Treatment']!="Other")]
sample_list = list(set(chemoradio_df.index)&set(tmt_group_df.index))

sub_group_df = tmt_group_df.loc[sample_list]
sub_group_df.columns = ['tmt_group']


sub_group_df["Progress"] = PFS_df.loc[sub_group_df.index,'是否进展（0，否；1，是）（截止到2023年7月24日）']
sub_group_df["PFS"] = PFS_df.loc[sub_group_df.index,'PFS（截止到2023年7月24日）']

sub_group_df["death_or_not"] = PFS_df.loc[sub_group_df.index,'是否死亡（0，否；1，是）（截止至2023年7月24日）']
sub_group_df["OS"] = PFS_df.loc[sub_group_df.index,'OS（月，截止到2023年7月24日）']




sub_group_df1 = sub_group_df.dropna(how="any")
sub_group_df = sub_group_df.dropna(how="any")


        

def gene_survival(chemo_radio_tmt_df,gene,type_=""):
    

    tmt_gene_df = chemo_radio_tmt_df.loc[gene].to_frame()
    tmt_gene_df["PFS"] = sub_group_df.loc[tmt_gene_df.index,"PFS"]
    tmt_gene_df['Progress'] = sub_group_df.loc[tmt_gene_df.index,"Progress"]
    tmt_gene_df = tmt_gene_df.sort_values(by=gene)
    
    min_half = 10
    
#    for  half_  in range(10,50):
    
    
    half_ =int(len(tmt_gene_df)/2)
    
    tmt_gene_df['group']=["Low"]*half_ +['High']*(len(tmt_gene_df)-half_)
    
    survival_df = tmt_gene_df
    
#    survival_df["tmt_proteomic_group"] = tmt_group_df.loc[survival_df.index,"TMT_proteomic_subgroup"]
    print(survival_df)
    color_dict={"High":"red","Low":"blue"}
    
    plt.figure(figsize=(8,4))
    ax = plt.subplot(121)
    sns.violinplot(x='group',y=gene,data=tmt_gene_df)
    
    
    ax = plt.subplot(122)
    kmf1 = KaplanMeierFitter()
    for name, grouped_df in survival_df.groupby("group"):
        kmf1.fit(grouped_df["PFS"], grouped_df["Progress"], label=name+"(%s)"%(len(grouped_df)),
                alpha =0.1)

        kmf1.plot(ax=ax,show_censors=True,ci_show=False,color=color_dict[name],
              censor_styles={'ms': 6},linewidth=0.5)
    plt.rcParams['xtick.labelsize']=12
    plt.rcParams['ytick.labelsize']=12
    plt.xlabel("PFS(month)",size=12)
    plt.ylabel("Progress free ratio",size=12)
 #   ax.yaxis.set_minor_locator(AutoMinorLocator())
#    ax.xaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim([0.1,1])
    
    #from lifelines.statistics import multivariate_logrank_test
    T1 =  survival_df.loc[survival_df["group"]=="High"]["PFS"].values
    E1 = survival_df.loc[survival_df["group"]=="High"]["Progress"].values
    T2 = survival_df.loc[survival_df["group"]=="Low"]["PFS"].values
    E2 = survival_df.loc[survival_df["group"]=="Low"]["Progress"].values
    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    
    pvalue = results.p_value
        
    ax.text(10,0.6,"pvalue:%s"%(str("{:.2e}".format(pvalue))),size=12)
    plt.tight_layout()
    os.chdir(r"G:\cervical_cancer_multiomics\results\chemoradio_survival")
    
    survival_df.to_excel("%s_%s_exp.xlsx"%(gene,type_))
    
    
    
    
#    plt.savefig("%s_z_boxplot_logrank.pdf"%(gene),dpi=300)
#    plt.savefig("%s_z_boxplot_logrank_%s.pdf"%(gene.replace(":","_"),type_),dpi=300,bbox_inches="tight")
    plt.close()
    
    return gene,pvalue


tmt_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_protein\tmt_z_statistic_table_remove_error_sample.csv",index_col=0,sep="\t")
chemo_radio_tmt_df = tmt_df[sub_group_df.index]
chemo_radio_tmt_df = chemo_radio_tmt_df.dropna(thresh=32)
gene_survival(chemo_radio_tmt_df,"PRKCB",type_="tmt")
gene_survival(chemo_radio_tmt_df,"FOSL2",type_="tmt")


RNA_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\tissue_RNA\summary\tpkm_remove_error_sample.xls",index_col=0,sep="\t")

chemo_radio_tmt_df = RNA_df[sub_group_df.index]
chemo_radio_tmt_df = chemo_radio_tmt_df.dropna(thresh=32)
gene_survival(chemo_radio_tmt_df,"PRKCB",type_="RNA")
gene_survival(chemo_radio_tmt_df,"FOSL2",type_="RNA")




DIA_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\expression_table\dia_z_statistic_table_removed_11_error_sample.xls",index_col=0,sep="\t")

chemo_radio_tmt_df = DIA_df[sub_group_df.index]
chemo_radio_tmt_df = chemo_radio_tmt_df.dropna(thresh=32)
gene_survival(chemo_radio_tmt_df,"PRKCB",type_="DIA")
gene_survival(chemo_radio_tmt_df,"FOSL2",type_="DIA")





tmt_phos_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_p_protein\tmt_phosprotein_normalized_log2_statistic_table_remove_error_sample_tumor.csv",index_col=0,sep=",")

chemo_radio_tmt_df = tmt_phos_df[sub_group_df.index]
chemo_radio_tmt_df = chemo_radio_tmt_df.dropna(thresh=32)
gene_survival(chemo_radio_tmt_df,"PRKCB",type_="tmt_phos")
gene_survival(chemo_radio_tmt_df,"FOSL2",type_="tmt_phos")




tmt_phos_site_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_p_protein\tmt_p_normalized_log2_statistic_table_remove_error_sample.csv",index_col=0,sep=",")

chemo_radio_tmt_df = tmt_phos_site_df[sub_group_df.index]
chemo_radio_tmt_df = chemo_radio_tmt_df.dropna(thresh=32)
#gene_survival(chemo_radio_tmt_df,"PRKCB:S11",type_="tmt_phos_site")
#gene_survival(chemo_radio_tmt_df,"PRKCB:T642",type_="tmt_phos_site")
gene_survival(chemo_radio_tmt_df,"FOSL2",type_="tmt_phos")


tmt_ace_site_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_ace_protein\tmt_ace_normalized_log2_statistic_table_remove_error_sample.csv",index_col=0,sep=",")
chemo_radio_tmt_df = tmt_ace_site_df[sub_group_df.index]
chemo_radio_tmt_df = chemo_radio_tmt_df.dropna(thresh=32)
#gene_survival(chemo_radio_tmt_df,"PRKCB:S11",type_="tmt_phos_site")
#gene_survival(chemo_radio_tmt_df,"PRKCB:T642",type_="tmt_phos_site")
gene_survival(chemo_radio_tmt_df,"FOSL2:K222",type_="tmt_ace_site")




tmt_ace_df = pd.read_csv(r"G:\cervical_cancer_multiomics\results\TMT_ace_protein\tmt_acetylprotein_normalized_log2_statistic_table_remove_error_sample.csv",index_col=0,sep=",")
chemo_radio_tmt_df = tmt_ace_df[sub_group_df.index]
chemo_radio_tmt_df = chemo_radio_tmt_df.dropna(thresh=32)
#gene_survival(chemo_radio_tmt_df,"PRKCB:S11",type_="tmt_phos_site")
#gene_survival(chemo_radio_tmt_df,"PRKCB:T642",type_="tmt_phos_site")
gene_survival(chemo_radio_tmt_df,"FOSL2",type_="tmt_ace")
#gene_survival(chemo_radio_tmt_df,"PRKCB",type_="tmt_ace")
gene_survival(chemo_radio_tmt_df,"EP300",type_="tmt_ace")












