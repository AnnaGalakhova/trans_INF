#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:03:52 2024

@author: annagalakhova
"""
#load packages 
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
#%% define a function to do LMM on correlation
def XY_LMM_correlation(df_file_metadata, Xcorrelate,Ycorrelate):
    df_file_metadata = df_file_metadata.dropna(subset=[Ycorrelate, Xcorrelate])
    df_file_metadata[Xcorrelate] = pd.to_numeric(df_file_metadata[Xcorrelate], errors='coerce')
    # LMM fitting
    md = smf.mixedlm(f"{Ycorrelate} ~ {Xcorrelate}", df_file_metadata, groups=df_file_metadata["NewLabel"])
    mdf = md.fit()
    print(mdf.summary())
    p_value = mdf.pvalues[Xcorrelate]
    # Annotate p-value
    print(p_value)
    return mdf.summary(), p_value
#%%load in the data for mouse
mouse_data = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/old version/PatchSeq/anna_220919_genes_metadata.csv')
#%% plot the violins
mouse_data['cluster_new']=list(map(lambda x: x.split()[-1], mouse_data.cluster_label))

cluster_palette = {'Adamts2': '#BFD579', 'Agmat': '#5DB6AF', 'Rrad': '#94BDD6'}

# plot the grey violins and donors in there 
plt.figure()
ax = sns.violinplot(data=mouse_data, x='cluster_new', y='Grm3', legend=None, inner = 'quart',
                order=cluster_palette.keys(), palette=cluster_palette, linewidth = 1)
                
maxy = plt.ylim()[1]
x_positions = np.arange(len(cluster_palette))
#ax= sns.stripplot(data=human_data, x='NewLabel', y='log10_GRM3',marker='o',order=cluster_palette.keys(), color='black', size=5, facecolors='none')
jitter_amount = 0.05 # Adjust this value for more or less jitter
for x_pos, label in zip(x_positions, cluster_palette.keys()):
    subset = mouse_data[mouse_data['cluster_new'] == label]
    # Add jitter by adding a small amount of random noise to the x positions
    jittered_x = np.random.normal(x_pos, jitter_amount, size=len(subset))
    plt.scatter(jittered_x, subset['Grm3'], marker='o', edgecolor='black', facecolor='none', s=20)
    median = subset['Grm3'].median()
    iqr = subset['Grm3'].quantile(0.75) - subset['Grm3'].quantile(0.25)

    # Annotate the plot with median and IQR
    plt.text(x_pos, maxy, f'Median: {median:.2f}\nIQR: {iqr:.2f}', ha='center', va='bottom')
    
    
plt.ylim(0, maxy)

plt.ylabel(r'Grm3 expression log$_{10}$(CPM+1)')
plt.xlabel('L2/3 Cell type')
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

plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/PS_mouse_Grm3_types_violin.eps')


#%% histograms for mouse
counts=[]
for cluster in cluster_palette.keys():
    type_data = mouse_data[mouse_data['cluster_new'] == cluster]
    counts_loc,_,_=plt.hist(type_data['Grm3'], bins=30)
    counts.append(max(counts_loc))
    plt.close()

ylim=max(counts)

#plt.figure() 
for cluster in cluster_palette.keys():
    # Create a new figure
    plt.figure() 
    # Plot histograms for =and 'GRM3'
    type_data = mouse_data[mouse_data['cluster_new'] == cluster]
    sns.histplot(type_data['Grm3'], bins=30, kde=True,  color=cluster_palette[cluster], alpha=1, label=cluster)
    #sns.kdeplot(type_data['Grm3'], color=cluster_palette[cluster], label=cluster, linestyle = '-')
    # Add title and legend to the figure
   
    plt.legend()
    plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/mouse/patchseq/'+'Grm3_'+'hist_tea_types'+f'{cluster}'+'.eps')
    #plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/mouse/'+'Grm2_'+'hist_tea_types'+'.eps')
    # Show the plot for the current class label
    plt.show()
    
del counts, cluster, non_zero_data, counts_loc


#%% load in human data
human_data=pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/patchseq data_HARpaper/cpm_all_data.csv')
human_data.at[human_data['NewLabel'] == 1, 'NewLabel'] = 'LTK'
human_data.at[human_data['NewLabel'] == 2, 'NewLabel'] = 'L2 FREM3'
human_data.at[human_data['NewLabel'] == 3, 'NewLabel'] = 'L3 FREM3'
human_data.at[human_data['NewLabel'] == 4, 'NewLabel'] = 'GLP2R'
human_data.at[human_data['NewLabel'] == 5, 'NewLabel'] = 'CARM1P1'
human_data.at[human_data['NewLabel'] == 6, 'NewLabel'] = 'COL22A1'

#%% log transform the data for GRM3
human_data['GRM3'] =human_data['GRM3']+1
human_data['log10_GRM3']= np.log10(human_data['GRM3'])
#add other parameters
features = pd.read_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/APtable.csv')

features_table = pd.merge(human_data, features, on='sample_id', how='inner')
features_table.to_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/data_correlations.csv')
#%% plot the expressions per type
cluster_palette = {'LTK': '#BFD579', 'GLP2R': '#5DB6AF', 'L2 FREM3': '#94BDD6', 'L3 FREM3': '#6B95CA', 'CARM1P1': '#A36192', 'COL22A1': '#EE75A2'}

# plot the grey violins and donors in there 
plt.figure()
ax = sns.violinplot(data=human_data, x='NewLabel', y='log10_GRM3', legend=None, inner = 'quart',
                order=cluster_palette.keys(), palette=cluster_palette, linewidth = 1)
                
maxy = plt.ylim()[1]
x_positions = np.arange(len(cluster_palette))
#ax= sns.stripplot(data=human_data, x='NewLabel', y='log10_GRM3',marker='o',order=cluster_palette.keys(), color='black', size=5, facecolors='none')
jitter_amount = 0.05 # Adjust this value for more or less jitter
for x_pos, label in zip(x_positions, cluster_palette.keys()):
    subset = human_data[human_data['NewLabel'] == label]
    # Add jitter by adding a small amount of random noise to the x positions
    jittered_x = np.random.normal(x_pos, jitter_amount, size=len(subset))
    plt.scatter(jittered_x, subset['log10_GRM3'], marker='o', edgecolor='black', facecolor='none', s=20)
    median = subset['log10_GRM3'].median()
    iqr = subset['log10_GRM3'].quantile(0.75) - subset['log10_GRM3'].quantile(0.25)

    # Annotate the plot with median and IQR
    plt.text(x_pos, maxy, f'Median: {median:.2f}\nIQR: {iqr:.2f}', ha='center', va='bottom')
    
    
plt.ylim(0, maxy)

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

plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/PS_human_GRM3_types_violin.eps')


#%% plot the data in histogram
#drop the non correct spike
features_table=features_table.drop(features_table['AmpFAPthresh'].idxmin())
counts=[]
for cluster in cluster_palette.keys():
    type_data = human_data[human_data['NewLabel'] == cluster]
    counts_loc,_,_=plt.hist(type_data['log10_GRM3'], bins=30)
    counts.append(max(counts_loc))
    plt.close()

ylim=max(counts)

plt.figure() 
for cluster in cluster_palette.keys():
    # Create a new figure
    #plt.figure() 
    # Plot histograms for =and 'GRM3'
    type_data = human_data[human_data['NewLabel'] == cluster]
    #sns.histplot(type_data['log10_GRM3'], bins=30, kde=True,  color=cluster_palette[cluster], alpha=1, label=cluster)
    sns.kdeplot(type_data['log10_GRM3'], color=cluster_palette[cluster], bw_adjust=1,label=cluster, linestyle = '-')
    # Add title and legend to the figure
    plt.xlim(left=0)
    plt.legend()
    #plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/human/patchseq/'+'GRM3_'+'hist_mtg_types'+f'{cluster}'+'.eps')
    #plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/mouse/'+'Grm2_'+'hist_tea_types'+'.eps')
    # Show the plot for the current class label
    plt.show()

#%% PLOT AND CALCULATE CORRELATIONS OF EXPRESSION TO PARAMETERS
correlations=pd.DataFrame(columns = ['param1', 'param2', 'test', 'corr', 'p', 'sign?', 'normality', 'N', 'null?', 'mode', 'lmm_p', 'lmm_sign'])
cluster_palette = {'LTK': '#BFD579', 'GLP2R': '#5DB6AF', 'L2 FREM3': '#94BDD6', 'L3 FREM3': '#6B95CA', 'CARM1P1': '#A36192', 'COL22A1': '#EE75A2'}
params = ['TDL_um', 'Rheobase', 'vmbaseM', 'InputR', 'FreqMax', 'FreqTrSwp', 'ThreshFrstAP', 'FAPbasetothresh', 'AmpFAPthresh', 'UpStrkFrstAP', 'DwnStrkFrstAP']
for param in params:
    #drop nans for the correlation
    columns_to_check = ['log10_GRM3', param]
    for column in columns_to_check:
        has_nans = features_table[column].isna().any()
        has_infs = np.isinf(features_table[column]).any()
    data=features_table.replace([np.inf, -np.inf], np.nan).dropna(subset=columns_to_check)
    #check normality
    normal = (stats.normaltest(data['log10_GRM3'])[1]>0.05) & (stats.normaltest(data[param])[1]>0.05)
    # Plot
    plt.figure()
    sns.scatterplot(data=data, x='log10_GRM3', y=param, hue='NewLabel',palette=cluster_palette)
    ax = sns.regplot(data=data, x='log10_GRM3', y=param, color='black', scatter=False)
    plt.xlabel(r'GRM3 expression log$_{10}$(CPM+1)')
    lmm_results, lmm_p =XY_LMM_correlation(features_table, 'log10_GRM3',param)
    if lmm_p <0.05:
        lmm_sign  ='yes'
    else: 
        lmm_sign= 'NS'
    #stats
    alpha = 0.05
    if normal == True:
        pcorrelation, pp_value = stats.pearsonr(data['log10_GRM3'],data[param])
        if pp_value < alpha:
            sign = 'yes'
        else:
            sign =  'NS'
        plt.title(f'Pearsons R = {pcorrelation:.3f}, p = {pp_value:.3g} \n normal')
        correlations=correlations.append({'param1' :'log10_GRM3' , 'param2':param, 'test':'Pearson', 'corr':pcorrelation, 'p':pp_value, 'sign?':sign,'normality':normal,'N':len(data), 'null?':'yes', 'mode':'all_data', 'lmm_p':lmm_p, 'lmm_sign':lmm_sign},ignore_index=True)
    elif normal == False:
        scorrelation, sp_value = stats.spearmanr(data['log10_GRM3'], data[param])
        if sp_value < alpha:
            ssign = 'yes'
        else:
                ssign =  'NS'
        plt.title(f'Spearman R = {scorrelation:.3f}, p = {sp_value:.3g} \n not normal')
        correlations=correlations.append({'param1' :'log10_GRM3' , 'param2':param, 'test':'Spearman', 'corr':scorrelation, 'p':sp_value,'sign?':ssign, 'normality':normal,'N':len(data), 'null?':'yes','mode':'all_data', 'lmm_p':lmm_p, 'lmm_sign':lmm_sign},ignore_index=True)
    plt.savefig('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/'+'GRM3_vs_'f'{param}_all_cells_null_colored'+'.eps')    
correlations['p_adj'] = multipletests(correlations['p'], method='bonferroni', alpha=alpha)[1]
correlations['sign_adj'] = correlations['p_adj'].apply(lambda x: 'yes' if x < alpha else 'NS')
    
    
correlations.to_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/correlations.csv')
    
#%% linear mixed model
    
lmm_results=XY_LMM_correlation(features_table, 'log10_GRM3',param)

#%% subdivide into types
correlationsi=pd.DataFrame(columns = ['param1', 'param2', 'test', 'corr', 'p', 'sign?', 'normality', 'N', 'null?', 'mode'])
params = ['TDL_um', 'Rheobase', 'vmbaseM', 'InputR', 'FreqMax', 'FreqTrSwp', 'ThreshFrstAP', 'AmpFAPthresh', 'UpStrkFrstAP', 'DwnStrkFrstAP']
for cluster in features_table.NewLabel.unique():
    for param in params:
        columns_to_check = ['log10_GRM3', param]
        for column in columns_to_check:
            has_nans = features_table[column].isna().any()
            has_infs = np.isinf(features_table[column]).any()
        data=features_table.replace([np.inf, -np.inf], np.nan).dropna(subset=columns_to_check)
        to_test=data[data.NewLabel == cluster]
        plt.figure()
        sns.scatterplot(data=to_test, x='log10_GRM3', y=param, color=cluster_palette[cluster])
        ax = sns.regplot(data=to_test, x='log10_GRM3', y=param, color='black', scatter=False)
        plt.xlabel(r'GRM3 expression log$_{10}$(CPM+1)')
        try:
            normal = (stats.normaltest(to_test['log10_GRM3'])[1]>0.05) & (stats.normaltest(to_test[param])[1]>0.05)
            if normal == True:
                corr, pval = stats.pearsonr(to_test['log10_GRM3'],to_test[param])
                if pval < alpha:
                    sign = 'yes'
                else: 
                    sign = 'NS'
                plt.title(f'Pearson R = {corr:.3f}, p = {pval:.3g} \n normal')
                correlationsi=correlationsi.append({'param1' :'log10_GRM3' , 'param2':param, 'test':'Pearson', 'corr':corr, 'p':pval, 'sign?':sign,'normality':normal,'N':len(to_test), 'null?':'yes', 'mode':cluster},ignore_index=True)
            else:
                corr, pval = stats.spearmanr(to_test['log10_GRM3'],to_test[param])
                if pval < alpha:
                    sign = 'yes'
                else: 
                    sign = 'NS'
                plt.title(f'Spearman R = {corr:.3f}, p = {pval:.3g} \n not normal')
                correlationsi=correlationsi.append({'param1' :'log10_GRM3' , 'param2':param, 'test':'Spearman', 'corr':corr, 'p':pval, 'sign?':sign,'normality':normal,'N':len(to_test), 'null?':'yes', 'mode':cluster},ignore_index=True)
        except Exception as e:
                print(f"Skipping {param} in cluster {cluster} due to error: {e}")
                plt.title('no tests')
                correlationsi=correlationsi.append({'param1' :'log10_GRM3' , 'param2':param, 'test':'no stats', 'corr':np.nan, 'p':np.nan, 'sign?':np.nan,'normality':np.nan,'N':len(to_test), 'null?':'yes', 'mode':cluster},ignore_index=True)
correlationsi['p_adj'] = multipletests(correlationsi['p'], method='bonferroni', alpha=alpha)[1]
correlationsi['sign_adj'] = correlationsi['p_adj'].apply(lambda x: 'yes' if x < alpha else 'NS')       
        
        
        

correlationsi.to_csv('/Users/annagalakhova/Library/Mobile Documents/com~apple~CloudDocs/PhD INF CNCR VU/DATA/GRM3 project/genes/eps/figure1/correlations_types.csv')
    
    
    
    
    