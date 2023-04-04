# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 18:15:15 2022

@author: ftead
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import math
import re
#from matplotlib_venn import venn2
import gseapy as gp

#%%
ev0_evp = pd.read_csv('ev0-evp.csv', engine = 'python')
ab0_abp = pd.read_csv('ab0-abp.csv', engine = 'python')
#%%

def group_them_2(row):
    if row['padj'] < 0.05:    
        if row['log2FoldChange'] < -0.58 or row['log2FoldChange'] > 0.58:
            return 'segnificantly differentially expressed'
        else:
            return 'segnificant and not differentially expressed'
    else:
        return 'not segnificant'
        
#%%

ev0_evp['-log(padj)'] = np.log10(ev0_evp['padj'])*-1
ab0_abp['-log(padj)'] = np.log10(ab0_abp['padj'])*-1


#%%

ev0_evp['group'] = ev0_evp.apply(group_them_2,axis = 1)
ab0_abp['group'] = ab0_abp.apply(group_them_2,axis = 1)
#%%
ev0_evp_seg = ev0_evp.loc[ev0_evp['group'] == 'segnificantly differentially expressed']
ab0_abp_seg = ab0_abp.loc[ab0_abp['group'] == 'segnificantly differentially expressed']
#%%
ev0_evp_seg.to_csv('ev0_evp_seg.csv')
ab0_abp_seg.to_csv('ab0_abp_seg.csv')


#%%


plt.figure(figsize=(10,15), dpi = 300)
sns.set(color_codes = True)
sns.set_style("whitegrid")
ax1 = sns.scatterplot(x = 'log2FoldChange', y = '-log(padj)', hue = 'group', palette = "Blues", data = ab0_ev0)
ax1.set(ylim =(-0.1,6), xlim = (-6,6), xlabel='log2(fold change)', ylabel = '-log(p Adjusted)')
ax1.axhline(np.log10(0.05)*-1, ls = '--', color = 'black')

plt.tight_layout()
plt.savefig('ab0_ev0_vol.png',bbox_inches='tight')

#%%


libs = gp.get_library_name(organism="Worm")

#%%

gene_sets = ['GO_Biological_Process_2018','GO_Cellular_Component_2018',
             'GO_Molecular_Function_2018','KEGG_2019',
             'RNAi_Phenotypes_WormBase_2018','Human_Diseases_from_WormBase_2018']

#%%

ab0_ev0_diff_seg = ab0_ev0.loc[ab0_ev0['group'] == 'segnificantly differentially expressed']

#%%

ab0_ev0_diff_seg_enr = gp.enrichr(gene_list = list(ab0_ev0_diff_seg['Gene_name'].unique()), description='test_name',organism='worm', gene_sets= gene_sets, cutoff = 0.5)

#%%

ev0_evp_results = ev0_evp_diff_seg_enr.results
ab0_abp_results = ab0_abp_diff_seg_enr.results

#%%

ev0_evp_results.to_csv('ev0_evp_results_GSEA.csv')
ab0_abp_results.to_csv('ab0_abp_results_GSEA.csv')

#%%
ev0_evp_results_seg = ev0_evp_results.loc[ev0_evp_results['Adjusted P-value']<0.05]
ab0_abp_results_seg = ab0_abp_results.loc[ab0_abp_results['Adjusted P-value']<0.05]
ev0_evp_results_seg['-log(Adjusted P-value)'] =np.log10(ev0_evp_results_seg['Adjusted P-value'])*-1
ab0_abp_results_seg['-log(Adjusted P-value)'] =np.log10(ab0_abp_results_seg['Adjusted P-value'])*-1

#%%
fig, ax = plt.subplots()
plt.figure(figsize=(6,8), dpi = 300)
sns.set_theme()
sns.set(font_scale = 1.4)
sns.set_style("whitegrid")
ax = sns.barplot(data = ev0_evp_results_seg, x = '-log(Adjusted P-value)', y = 'Term', hue = 'Gene_set', palette = 'mako',edgecolor="black")
plt.xlim(0)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('ev0_evp_results_seg.png',dpi = 300,bbox_inches='tight',transparent=True)


fig, ax = plt.subplots()
plt.figure(figsize=(6,8), dpi = 300)
sns.set_theme()
sns.set(font_scale = 1.4)
sns.set_style("whitegrid")
ax = sns.barplot(data = ab0_abp_results_seg, x = '-log(Adjusted P-value)', y = 'Term', hue = 'Gene_set', palette = 'mako',edgecolor="black")
plt.xlim(0)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('ab0_abp_results_seg.png',dpi = 300,bbox_inches='tight',transparent=True)




