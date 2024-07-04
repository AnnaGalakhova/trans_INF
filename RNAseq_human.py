#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 12:06:56 2024

@author: annagalakhova
"""

"""
Created on Wed Mar 1 10:16:26 2023

This script is used to analyze and plot subsets of the AIBS RNA-seq database (https://portal.brain-map.org/atlases-and-data/rnaseq).
MULTIPLE CORTICAL AREAS - SMART-SEQ (2019) full matrix (exons+introns) is used to sub-select expression values
for specific genes of interest per region or cell classes. 
Script structure: 
24 - 38   selecting the data based on gene set of interest.
39 - 43   create datasets based on gene-set of interest (can be save and used later)    
44 - 61   Load in selected data and add values to metadata 
62 - 69   Group data by donor, allowing for plots of donordata seperate 
70 - 97  plot the data
98 - end statistical evalutation 

note. script on normaliztion of this dataset can be found in FINAL_RNA_seq_data_nromalization.py 

@author: Stan Driessens & Anna Galakhova 
"""
#load packages 
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np
#read and import data from regions 
#import metadata
metadata = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/Multiple areas SMART_Seq/metadata.csv')
#metadata = metadata.set_index('sample_name')
#clear the metadata
metadata_cleared = metadata[metadata['class_label'].notnull()]
#%% create gene-set data from normalized dataframe
geneset = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/geneset.txt', names=['gene'], header=None   , quotechar="'")
#in case there are duplicates present in the gene set, remove and update the gene set 
geneset = geneset.gene.drop_duplicates()
#convert geneset to a dataframe to improve compatibillity with pandas 
geneset = pd.DataFrame(geneset)
#%% load and filter data based on gene-set (ONLY NEED TO DO THIS ONCE, THEN SAVE PRE MADE DATA AS .CSV)
# selecting expression data of the genes present in the gene set of interest (geneset)
skiprows_set1 = metadata.loc[~metadata.sample_name.isin(metadata_cleared.sample_name)].index
#load in the data, NOTE. here we used pre-normalized log10 CPM data

data_log10CPM = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/Multiple areas SMART_Seq/matrix_normalised_log10.csv',
          usecols=geneset) #skiprows=skiprows_set1, 

data_log10CPM =data_log10CPM[data_log10CPM.sample_name.isin(metadata_cleared.sample_name)]
#save the pre made data 
data_log10CPM.to_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/RNAseq_human_areas.csv')

#%%load in the pre-made normalized and pre-assigned gene-set data 
data_log10CPM = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/RNAseq_human_areas.csv')
#concatenate metadata to avg values from this geneset 
b = pd.merge(metadata_cleared, data_log10CPM, on='sample_name', how='inner')
#create new column for region values where M1 and S1 sub-regions are merged to new M1 and S1
b.loc[b['region_label'].str.contains('M1'), 'region']  = 'M1'
b.loc[b['region_label'].str.contains('MTG'), 'region']  = 'MTG'
b.loc[b['region_label'].str.contains('A1'), 'region']  = 'A1'
b.loc[b['region_label'].str.contains('CgG'), 'region']  = 'CgG'
b.loc[b['region_label'].str.contains('S1'), 'region']  = 'S1'
b.loc[b['region_label'].str.contains('V1'), 'region']  = 'V1'
#select only mtg
mtg=b.loc[b['region'] == 'MTG']
#%% plotting in violin for all donors;all regions
#Figure 1_panel
gene='HTR1A'
plt.figure()
ax = sns.violinplot(data=mtg, x='class_label', y=mtg[gene], legend=None, inner = "quart", 
                    order=[  'GABAergic', 'Non-neuronal','Glutamatergic'], color='lightgrey', linewidth = 0.5)
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
#plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/RNAseq_'+str(gene)+'_mtg_violins_glu_last.eps')

#%% look at the different clusters in layer 23
#read and import data for MTG RNA-seq data  
#import metadata
metadata = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_samples-columns.csv')
metadata=metadata.rename(columns={'class': 'class_label'})
metadata_cleared = metadata[metadata['class_label'].notnull()]
gene_metadata=pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_genes-rows.csv')
#load int the genes of interest
gene_selection = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/geneset.txt', names=['gene'], header=None   , quotechar="'")

usecols = list(metadata_cleared.sample_name)
usecols.append('Unnamed: 0')

#%%load in reads from exon and intron 
exons = pd.read_csv("/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_exon-matrix.csv", 
                        usecols=usecols, index_col = 0)
introns = pd.read_csv("/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/transcriptomics/RNAseq/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_intron-matrix.csv", 
                        usecols=usecols, index_col = 0)

#%%process the data - add the exons and introns and normalise to CPM and log transform
data = exons + introns
#transpose the data to make genes columns 
data_T = data.T
#create a column with total reads per cell
data_T['total_reads'] = data_T.sum(axis=1)
#divide all reads by total reads 
data_T_c = data_T.div(data_T['total_reads'], axis=0)
#multiply by 1e6 to get 'counts per million'
data_T_cpm = data_T_c*1e6
#add 1 to prevent inf values when log transforming
data_T_cpm1 = data_T_cpm+1
#transorm 10 log 
data_log10cpm = np.log10(data_T_cpm1)
#check if the correct genes are in the dataframe 
selection_i = gene_metadata[gene_metadata.gene.isin(gene_selection.gene)]
selection = list(selection_i.entrez_id)
#subselect the data to only contain genes of interest
output = data_log10cpm[selection]
output=output.reset_index()
for column in output.columns:
    if column in gene_metadata.entrez_id:
        new_name = gene_metadata[gene_metadata.entrez_id == column]['gene'].values[0]
        output=output.rename(columns ={column:new_name})
output=output.rename(columns ={'index':'sample_name'})        
#add the corresponding metadata
output=pd.merge(output, metadata_cleared, on='sample_name')
#%% save the pre made data 
output.to_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/RNAseq_human_MTG.csv')
#%% process the data - subselecting the types of interest, giving a new label takig into account layer
final_data = output.loc[(output['brain_subregion'] == 'L2') | (output['brain_subregion'] == 'L3')]
#select the t-types of interest
final_data_types = final_data[(final_data['cluster'].str.contains('FREM3')) | (final_data['cluster'].str.contains('GLP2'))
                  | (final_data['cluster'].str.contains('LTK')) | (final_data['cluster'].str.contains('CARM1P1'))
                  | (final_data['cluster'].str.contains('COL22')) | (final_data['cluster'].str.contains('Astro'))]
#subselect frem into deep and superficial 
final_data_types['cluster_new']=list(map(lambda x: x.split()[3], final_data_types.cluster))
final_data_types.at[final_data_types['cluster_new'] == 'GFAP', 'cluster_new'] = 'Astrocyte'
final_data_types.loc[(final_data_types.cluster_new.values == 'FREM3') & (final_data_types.brain_subregion.values == 'L2'), 'cluster_new'] = 'L2 FREM3'
final_data_types.loc[(final_data_types.cluster_new.values == 'FREM3') & (final_data_types.brain_subregion.values == 'L3'), 'cluster_new'] = 'L3 FREM3'


final_data_types.to_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/processed_RNAseq_human_MTG.csv')
#%%plot the data per donor and per class
cluster_palette = {'LTK': '#BFD579', 'GLP2R': '#5DB6AF', 'L2 FREM3': '#94BDD6', 'L3 FREM3': '#6B95CA', 'CARM1P1': '#A36192', 'COL22A1': '#EE75A2', 'GFAP' : 'lightgrey'}

gene='HTR2A'
plt.figure()
ax = sns.violinplot(data=final_data_types, x='cluster_new', y=final_data_types[gene], legend=None, inner = "quart", 
                order=['LTK', 'GLP2R', 'L2 FREM3', 'L3 FREM3' , 'CARM1P1', 'COL22A1'], palette=cluster_palette, linewidth = 1)
plt.ylim(bottom = 0)
plt.ylabel(r'GRM3 expression log$_{10}$(CPM+1)')
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

plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/RNAseq_'+str(gene)+'mtg_types_violin.eps')

#%%plot histograms
counts=[]
for cluster in cluster_palette.keys():
    non_zero_data = final_data_types[(final_data_types['GRM3'] != 0) & (final_data_types['cluster_new'] == cluster)]
    counts_loc,_,_=plt.hist(non_zero_data['GRM3'], bins=30)
    counts.append(max(counts_loc))
    plt.close()

ylim=max(counts)

#plt.figure() 
for cluster in cluster_palette.keys():
    # Create a new figure
    plt.figure() 
    # Plot histograms for =and 'GRM3'
    type_data = final_data_types[final_data_types['cluster_new'] == cluster]
    sns.histplot(type_data['GRM2'], bins=30, kde=True,  color=cluster_palette[cluster], alpha=0.4, label=cluster)
    sns.histplot(type_data['GRM3'], bins=30, kde=True,  color=cluster_palette[cluster], alpha=1, label=cluster)
    #sns.kdeplot(type_data['GRM2'], color=cluster_palette[cluster], label=cluster, linestyle = '--')
    #sns.kdeplot(type_data['GRM3'], color=cluster_palette[cluster], label=cluster, linestyle = '-')
    # Add title and legend to the figure
    #plt.ylim(0, ylim+10)
    #plt.title(f'{cluster} distribution')
    plt.legend()
    plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/'+'GRM2pale3dark_'+'hist_mtg_types'+f'{cluster}'+'_alldonors.eps')
    # Show the plot for the current class label
    plt.show()
    
del counts, cluster, non_zero_data, counts_loc
#%%statistical evalutation 


