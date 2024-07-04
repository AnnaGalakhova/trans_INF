#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 10:48:36 2022

This script is used to analysze and plot subsets of the AIBS RNA-seq database (https://portal.brain-map.org/atlases-and-data/rnaseq).
MTG - SMART-SEQ (2019) seperate matrices (exons and introns) is used to sub-select expression values
for specific genes of interest per transcriptomic type within MTG. 

@author: annagalakhova & stand
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns





#%% load int he region data for mouse
samples_metadata = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/Multiple areas_mouse/metadata.csv')
usecols = ['sample_name', 'Ryr2', 'Rtn1', 'Jph3', 'Jph4']
data = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/Multiple areas_mouse//matrix_log10_CPM.csv', usecols=usecols)
#%% subselecting and adding additional metadata to the output
metadata_cleared = samples_metadata[samples_metadata['class_label'].notnull()]
output = pd.merge(metadata_cleared, data, on='sample_name', how='inner')
regions = ['TEa-PERI-ECT','ACA','MOp','VISp','SSp','AUD']
output=output[output.region_label.isin(regions)]
output.at[output['region_label'] == 'TEa-PERI-ECT', 'region_label'] = 'TEa'
tea=output
#%% plot per region, and then zoom in into tea for cell classes
gene='Jph4'
plt.figure()

ax = sns.violinplot(data=tea, x='class_label', y=tea[gene], legend=None, inner = "quart", 
                    order=['Glutamatergic','GABAergic', 'Non_neuronal',], color='lightgrey', linewidth = 0.5)
plt.ylim(bottom = 0)
plt.ylabel(f'{gene} expression log$_{10}$(CPM+1)')
plt.xlabel('Class')
for l in ax.lines:
    l.set_linestyle('--')
    l.set_linewidth(0)
    l.set_color('red')
    l.set_alpha(0.5)
for l in ax.lines[1::3]:
    l.set_linestyle('-')
    l.set_linewidth(1)
    l.set_color('black')
    l.set_alpha(0.5)

plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/RNAseq_'+str(gene)+'_tea_violin_classes_nn_last.eps')

#%% LOAD AND PROCESS METADATA FILES for Tea region [homologous to MTG in humans]
#load the metadata to select samples 
samples_metadata = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/VISp_mouse/mouse_VISp_2018-06-14_samples-columns.csv')
samples_metadata['class'] = samples_metadata.rename(columns = {'class':'class_label'}, inplace = True)
#select samples based on layers  
layers = ['L2/3']
samples_metadata_L23 = samples_metadata[samples_metadata['brain_subregion'].isin(layers)] 
types = ['L2/3 IT VISp Adamts2', 'L2/3 IT VISp Agmat', 'L2/3 IT VISp Rrad', 'Astro Aqp4']

samples_metadata_L23types = samples_metadata_L23[samples_metadata_L23['cluster'].isin(types)]

#%% LOAD GENES METADATA AND THE GENE SET OF INTEREST 
#load the gene metadata, containing information about all the genes. 
gene_metadata = pd.read_csv("/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/VISp_mouse/mouse_VISp_2018-06-14_genes-rows.csv")

#load in the geneset from text_file NOTE: add quotechar if apostrophe's are present in the original text file (i.e. 'gene' instead of gene)
geneset = gene_metadata[gene_metadata['gene_symbol'].str.contains('Grm')]
gene_selection = list(geneset.gene_entrez_id)


#%% filters to open exons and introns
usecols = list(samples_metadata_L23types.sample_name)
usecols.append('Unnamed: 0')

exons = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/VISp_mouse/mouse_VISp_2018-06-14_exon-matrix.csv', 
                        usecols=usecols, index_col = 0)
introns = pd.read_csv("/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/VISp_mouse/mouse_VISp_2018-06-14_intron-matrix.csv",
                     usecols=usecols, index_col = 0)
#%% further process the table - transpose, normalise, add metadata
data = exons + introns
data_T = data.T
data_T['total_reads'] = data_T.sum(axis=1)
data_T_c = data_T.div(data_T['total_reads'], axis=0)
data_T_cpm = data_T_c*1e6
data_T_cpm1 = data_T_cpm+1
data_log10cpm = np.log10(data_T_cpm1)

output = data_log10cpm[gene_selection]
output=output.reset_index()

for column in output.columns:
    if geneset['gene_entrez_id'].isin([column]).any():
        print(type(column))
        new_name = geneset[geneset.gene_entrez_id == column]['gene_symbol'].values[0]
        print(type(new_name))
        output=output.rename(columns ={column:new_name})
output=output.rename(columns ={'index':'sample_name'})        
#add the corresponding metadata
output=pd.merge(output, samples_metadata_L23types, on='sample_name')
output=output.drop(columns='class')
output['cluster_new']=list(map(lambda x: x.split()[-1], output.cluster))
output.at[output['cluster_new'] == 'Aqp4', 'cluster_new'] = 'Astrocyte'

output.to_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/RNAseq_mouse_tea.csv')

#%%plot the data for t_types
cluster_palette = {'Adamts2': '#BFD579', 'Agmat': '#5DB6AF', 'Rrad': '#94BDD6'}
gene='Grm3'
# plot the grey violins and donors in there 

plt.figure()
ax = sns.violinplot(data=output, x='cluster_new', y=output[gene], legend=None, inner = "quart", 
                order=cluster_palette.keys(), palette=cluster_palette, linewidth = 1)
plt.ylim(bottom = 0)
plt.ylabel(r'Grm3 expression log$_{10}$(CPM+1)')
plt.xlabel('L2/L3 Cell type')
for l in ax.lines:
    l.set_linestyle('--')
    l.set_linewidth(0)
    l.set_color('red')
    l.set_alpha(0.5)
for l in ax.lines[1::3]:
    l.set_linestyle('-')
    l.set_linewidth(1)
    l.set_color('black')
    l.set_alpha(0.5)
    
    
plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/RNAseq_'+str(gene)+'tea_types_violin.eps')













