# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 08:50:04 2025

@author: AChao

This script contains code for processing the WW2DW results generated from the NTA WebApp for the purposes of the intro paper
"""

import time
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle 
import csv

time_list = []
start_time = time.time()
time_list.append(time.time())

'''
# Grab the System suitability chemicals
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 input files')
df1 = pd.read_csv('WW2DW_SSM_Amenable.csv')
'''

'''
##################################################################################################################
# Grab the initial MS1 chemical results from the Excel results file, pare down to SSM chemicals and save/export
##################################################################################################################

# Grab the WebApp MS1 chemical results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (tracer 7ppm MRL 10 parent 7 ppm) w CHM')
df2 = pd.read_excel('Example_nta_NTA_WebApp_results.xlsx', sheet_name='Chemical Results')
# 1/21/2025: Need to get most up to date metadata from Tony/API?


# Grab the WebApp MS1 Hazard Comparison Module results - Don't need to grab, already in the chemical results sheet
#df3 = pd.read_excel('WW2DW_tracer_7ppm_MRL_10_parent_6ppm_20250121_NTA_WebApp_results.xlsx', sheet_name='hcd_search')

# Alternatively grab SSM from the WebApp MS1 results which ran using SSM in the tracer input file
df1_b = pd.read_excel('Example_nta_NTA_WebApp_results.xlsx', sheet_name='Tracer Summary')
# Rename DTXSID column of tracer summary to not collide with DTXSID's of DSSTox/CHEM results
df1_b.rename(columns={'DTXSID':'Tracer DTXSID'}, inplace=True) 


df4 = pd.merge(df1_b, df2, how='left', on='Feature ID')

df4['SSM chemical'] = df4['Tracer DTXSID'] == df4['DTXSID']

# Export the MS1 results for just the SSM chemicals to reduce dataframe overhead
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df4.to_csv('WW2DW MS1 results (chemical and hazard) for SSM chemicals.csv', sep=',', encoding='utf-8', index=False)
'''


'''
##################################################################################################################
# Merge on additional searches/data onto initial exported chemical results from above
##################################################################################################################

# Read in MS1 results for just the SSM chemicals
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df_ms1 = pd.read_csv('WW2DW MS1 results (chemical and hazard) for SSM chemicals.csv')

# Read in AMOS search results for the SSM chemicals
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (tracer 7ppm MRL 10 parent 7 ppm) w CHM/AMOS search results of SSM')
# Read in the AMOS source counts
df_amos = pd.read_excel('AMOS search results of chemical results 20250123.xlsx', sheet_name='AMOS search results of chemical')
df_amos = df_amos.drop_duplicates(subset=['DTXSID'])
# Read in the AMOS methods/fact sheets/spectral counts
df_amos_b = pd.read_excel('AMOS search results of chemical results 20250123.xlsx', sheet_name='Records')
# Split into methods / fact sheets / spectra dataframes for separate counting
df_amos_methods = df_amos_b[df_amos_b['Record Type'] == 'Method']
df_amos_fact_sheets = df_amos_b[df_amos_b['Record Type'] == 'Fact Sheet']
df_amos_spectra = df_amos_b[df_amos_b['Record Type'] == 'Spectrum']
# Tabulate AMOS methods/fact sheets/spectral counts for a given DTXSID
methods_grouped = df_amos_methods.groupby('DTXSID').size()
fact_sheets_grouped = df_amos_fact_sheets.groupby('DTXSID').size()
spectra_grouped = df_amos_spectra.groupby('DTXSID').size()
df_amos_methods_grouped = methods_grouped.reset_index(name='AMOS methods count')
df_amos_fact_sheets_grouped = fact_sheets_grouped.reset_index(name='AMOS fact sheets count')
df_amos_spectra_grouped = spectra_grouped.reset_index(name='AMOS spectra count')

# Read in Comptox search results for the SSM chemicals on water lists
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (tracer 7ppm MRL 10 parent 7 ppm) w CHM/Comptox search results of SSM')
df_water_lists = pd.read_csv('CCD-Batch-Search_2025-01-23_04_22_45.csv')
# De-duplicate on DTXSID
df_water_lists = df_water_lists.drop_duplicates(subset=['DTXSID'])
df_water_lists = df_water_lists.drop(['INPUT', 'FOUND_BY', 'PREFERRED_NAME', 'DTXCID'], axis=1)
# Tabulate list counts for each chemical
df_water_lists['Presence in water lists count'] = df_water_lists.iloc[:, 1:].apply(lambda row: row.eq('Y').sum(), axis=1)

# Merge search results onto MS1 export
df_merge = pd.merge(df_ms1, df_amos[['DTXSID', 'Sources', 'Patents', 'Articles', 'PubMed Record Count']], how='left', on='DTXSID')
df_merge = pd.merge(df_merge, df_amos_methods_grouped, how='left', on='DTXSID')
df_merge = pd.merge(df_merge, df_amos_fact_sheets_grouped, how='left', on='DTXSID')
df_merge = pd.merge(df_merge, df_amos_spectra_grouped, how='left', on='DTXSID')
df_merge = pd.merge(df_merge, df_water_lists[['DTXSID', 'Presence in water lists count']], how='left', on='DTXSID')


# Export the merged results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df_merge.to_csv('WW2DW MS1 results (chemical and hazard) for SSM chemicals with search results included.csv', sep=',', encoding='utf-8', index=False)
'''


'''

# MS2 plots hereherehere
##################################################################################################################
# Merge MS2 results onto merged MS1 chemical results from above, debug/check for matches
##################################################################################################################
# Read in MS1 results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')
#df1 = df1.drop_duplicates(subset=['Feature ID'])

df1.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 

df1_pos = df1[df1['Ionization Mode'] == 'ESI+']
df1_neg = df1[df1['Ionization Mode'] == 'ESI-']

# Read in the MS2 results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/neg')
df_ms2_neg = pd.read_csv('20210131_WW2DW_MS2_neg_CFMID_results_neg.csv')
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/pos')
df_ms2_pos = pd.read_csv('20250131_WW2DW_MS2_pos_CFMID_results_pos.csv')


# # Do some manual checking of merge process
# df_ms2_neg = df_ms2_neg.drop_duplicates(subset=['ID'])
# df_ms2_pos = df_ms2_pos.drop_duplicates(subset=['ID'])

# ms1_mass_pos = df1_pos['Mass'].tolist()
# ms1_RT_pos = df1_pos['Retention Time'].tolist()
# ms1_mass_neg = df1_neg['Mass'].tolist()
# ms1_RT_neg = df1_neg['Retention Time'].tolist()

# ms2_mass_neg = df_ms2_neg['MASS_NEUTRAL'].tolist()
# ms2_RT_neg = df_ms2_neg['RT'].tolist()
# ms2_mass_pos = df_ms2_pos['MASS_NEUTRAL']
# ms2_RT_pos = df_ms2_pos['RT'].tolist()

# for mass_1, RT_1 in zip(ms1_mass_pos, ms1_RT_pos):
#     for mass_2, RT_2 in zip(ms2_mass_pos, ms2_RT_pos):
#         if abs(mass_1 - mass_2) < 0.05:
#             if abs(RT_1 - RT_2) < 0.5:
#                 print('Pos Match found between MS1:', mass_1, RT_1, ' and MS2:', mass_2, RT_2)

# mass_error = []
# RT_diff = []

# # Create function for calculating mass error
# def calc_mass_error_ppm(mass_1, mass_2):
#     mass_error = abs(mass_1 - mass_2)/mass_1*1000000
#     return mass_error

# for mass_1, RT_1 in zip(ms1_mass_neg, ms1_RT_neg):
#     for mass_2, RT_2 in zip(ms2_mass_neg, ms2_RT_neg):
#         #if abs(mass_1 - mass_2) < 0.05:
#         if calc_mass_error_ppm(mass_1, mass_2) < 5:
#             if abs(RT_1 - RT_2) < 0.5:
#                 #temp_mass_error = (mass_1 - mass_2) / mass_1 * 1000000
#                 mass_error.append(calc_mass_error_ppm(mass_1, mass_2))
#                 RT_diff.append(RT_1 - RT_2)
#                 print('neg Match found between MS1:', mass_1, RT_1, ' and MS2:', mass_2, RT_2)


matched_df_pos = df1_pos.merge(df_ms2_pos[
                [
                    "DTXCID",
                    "MASS_MGF",
                    "MASS_NEUTRAL",
                    "RT",
                    "SUM_SCORE",
                    "Q-SCORE",
                    "PERCENTILE",
                ]
            ],
            how="left",
            on="DTXCID",
        )

matched_df_neg = df1_neg.merge(df_ms2_neg[
                [
                    "DTXCID",
                    "MASS_MGF",
                    "MASS_NEUTRAL",
                    "RT",
                    "SUM_SCORE",
                    "Q-SCORE",
                    "PERCENTILE",
                ]
            ],
            how="left",
            on="DTXCID",
        )

# Select an SSM chemical
truncate = 0
#select_chem = 'atenolol'
select_chem = 'theophylline'




matched_pos_plot = matched_df_pos[['Chemical Name', 'DTXCID', 'SUM_SCORE']]

matched_pos_plot = matched_pos_plot[matched_pos_plot['Chemical Name'] == select_chem]
#df5_plot['PREFERRED_NAME'] = df5_plot['PREFERRED_NAME'].str[0:30]

matched_pos_plot = matched_pos_plot.drop_duplicates(subset='DTXCID', keep='first')

matched_pos_plot['SUM_SCORE'] = matched_pos_plot['SUM_SCORE'].fillna(0)

matched_pos_plot = matched_pos_plot.sort_values(by='SUM_SCORE', ascending=True)

# Identify the DTXCID associated with the SSM chemical
df1_pos_true = df1_pos[df1_pos['Chemical Name'] == select_chem]
SSM_true = df1_pos_true[df1_pos_true['SSM chemical'] == True]
SSM_DTXCID = SSM_true['DTXCID'].tolist()[0]


if truncate > 0:
    matched_pos_plot = matched_pos_plot.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
    SSM_font_size = 12
    label_sizes = 12
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    SSM_font_size = 22
    label_sizes = 22
    ax.tick_params(axis='y', labelsize=label_sizes - 4)
    ax.tick_params(axis='x', labelsize=label_sizes + 12)
    
#colors = ['blue'] * len()

#ax.barh(matched_pos_plot['DTXCID'], matched_pos_plot['SUM_SCORE'])
#bars = ax.barh(matched_pos_plot['DTXCID'], matched_pos_plot['SUM_SCORE'], color='skyblue')
bars = ax.barh(matched_pos_plot['DTXCID'], matched_pos_plot['SUM_SCORE'], alpha=0.3)


#plt.legend()
plt.ylabel('Candidate', fontsize=label_sizes+12)
plt.xlabel('MS2 Score', fontsize=label_sizes+12)

# # Set size and color of DTXCID corresponding to SSM chemical
# for label in ax.get_yticklabels():
#     if label.get_text() == SSM_DTXCID:
#         label.set_fontsize(SSM_font_size)  # Set font size to 16
#         label.set_color('green')

# Set size and color of DTXCID corresponding to SSM chemical
for bar, label in zip(bars, ax.get_yticklabels()):
    if label.get_text() == SSM_DTXCID:
        label.set_fontsize(SSM_font_size)  # Set font size to 16
        label.set_color('blue')
        bar.set_alpha(1)
        #bar.set_color('blue')


#ax.legend(fontsize=label_sizes + 12, loc='lower right')

title_string = 'MS2 Results for SSM Chemical: ' + select_chem
plt.title(title_string, fontsize=label_sizes + 12)
plt.show()

'''

'''
##################################################################################################################
# Merge MS2 results onto merged MS1 chemical results from above
# Combine Metadata and MS2 results into a stacked bar chart
# Right now only handles one mode of data
# I don't really like these
##################################################################################################################
# Read in MS1 results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')
#df1 = df1.drop_duplicates(subset=['Feature ID'])

df1.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 

df1_pos = df1[df1['Ionization Mode'] == 'ESI+']
df1_neg = df1[df1['Ionization Mode'] == 'ESI-']

# Read in the MS2 results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/neg')
df_ms2_neg = pd.read_csv('20210131_WW2DW_MS2_neg_CFMID_results_neg.csv')
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/pos')
df_ms2_pos = pd.read_csv('20250131_WW2DW_MS2_pos_CFMID_results_pos.csv')

matched_df_pos = df1_pos.merge(df_ms2_pos[["DTXCID", "MASS_MGF", "MASS_NEUTRAL", "RT", "SUM_SCORE", "Q-SCORE", "PERCENTILE",]], how="left", on="DTXCID")
matched_df_neg = df1_neg.merge(df_ms2_neg[["DTXCID", "MASS_MGF", "MASS_NEUTRAL", "RT", "SUM_SCORE", "Q-SCORE", "PERCENTILE",]], how="left", on="DTXCID")

matched_df_pos['Q-SCORE'] = matched_df_pos['Q-SCORE'].fillna(0)
matched_df_neg['Q-SCORE'] = matched_df_neg['Q-SCORE'].fillna(0)


# Select an SSM chemical
truncate = 0
#select_chem = 'warfarin'
#select_chem = 'trimethoprim'
#select_chem = 'ranitidine'
#select_chem = 'norfluoxetine'
select_chem = '10-hydroxy-amitriptyline'


# Start off with just grabbing a single feature
df1 = matched_df_pos[matched_df_pos['Chemical Name'] == select_chem]
# Remove feature duplicates if SSM chemical in more than one mode
df1 = df1.sort_values(by='Feature ID')
df1 = df1.drop_duplicates(subset='DTXSID', keep='first')

# Get leg 1 data
metadata_columns = ['Sources', 'Patents', 'Articles', 'PubMed Record Count', 'AMOS methods count', 'AMOS fact sheets count', 'AMOS spectra count', 'Presence in water lists count']
#df_pos = matched_df_pos.groupby('DTXCID').agg({col: 'sum' for col in metadata_columns})
df_pos = df1.groupby(['Feature ID', 'DTXCID']).agg({col: 'sum' for col in metadata_columns})
df_pos = pd.merge(df1[['Feature ID', 'DTXCID', 'Q-SCORE']], df_pos, on=['Feature ID', 'DTXCID'], how='right')
df_pos = df_pos.drop_duplicates(subset=['Feature ID', 'DTXCID'])
df_pos = df_pos.reset_index()


df_pos['Sources_norm'] = df_pos['Sources'] / df_pos['Sources'].max()
df_pos['Patents_norm'] = df_pos['Patents'] / df_pos['Patents'].max()
df_pos['Articles_norm'] = df_pos['Articles'] / df_pos['Articles'].max()
df_pos['PubMed Record Count_norm'] = df_pos['PubMed Record Count'] / df_pos['PubMed Record Count'].max()
df_pos['AMOS methods count_norm'] = df_pos['AMOS methods count'] / df_pos['AMOS methods count'].max()
df_pos['AMOS fact sheets count_norm'] = df_pos['AMOS fact sheets count'] / df_pos['AMOS fact sheets count'].max()
df_pos['AMOS spectra count_norm'] = df_pos['AMOS spectra count'] / df_pos['AMOS spectra count'].max()
df_pos['Presence in water lists count_norm'] = df_pos['Presence in water lists count'] / df_pos['Presence in water lists count'].max()

df_pos.fillna(0, inplace=True)

df_pos['Total_norm'] = df_pos['Sources_norm'] + df_pos['Patents_norm'] + df_pos['Articles_norm'] + df_pos['PubMed Record Count_norm'] + df_pos['AMOS methods count_norm'] + df_pos['AMOS fact sheets count_norm'] + df_pos['AMOS spectra count_norm'] + df_pos['Presence in water lists count_norm']

df_pos['Total_norm_norm'] = df_pos['Total_norm'] / df_pos['Total_norm'].max()

df_pos['metadata_plus_ms2'] = df_pos['Total_norm_norm'] + df_pos['Q-SCORE']
df_pos = df_pos.sort_values(by='metadata_plus_ms2', ascending=True)
#df_pos = df_pos.sort_values(by='Total_norm_norm', ascending=True)

df_plot = df_pos[['DTXCID', 'Q-SCORE', 'Total_norm_norm']]

if truncate > 0:
    df_plot = df_plot.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    ax.tick_params(axis='y', labelsize=22)  # Set y-axis tick label size to 14
# Initialize bottom for stacking
bottom = [0] * len(df_plot)

# Plot each column
for col in df_plot.columns[1:]:
    ax.barh(df_plot['DTXCID'], df_plot[col], left=bottom, label=col)
    bottom += df_plot[col]

plt.legend()
plt.ylabel('Candidate')
plt.xlabel('Total Metadata Score + MS2 Score')



title_string = 'Metadata/MS2 Results for SSM Chemical: ' + select_chem
plt.title(title_string)
plt.show()
'''


'''
##################################################################################################################
# Do some real bar chart plotting - Multi-bar chart
##################################################################################################################
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df_merge = pd.read_csv('WW2DW MS1 results (chemical and hazard) for SSM chemicals with search results included.csv')
# Start off with just grabbing a single feature
df5 = df_merge[df_merge['Chemical Name'] == '10-hydroxy-amitriptyline']
# Normalize all columns to graph
df5['Sources_norm'] = df5['Sources'] / df5['Sources'].max()
df5['Patents_norm'] = df5['Patents'] / df5['Patents'].max()
df5['Articles_norm'] = df5['Articles'] / df5['Articles'].max()
df5['PubMed Record Count_norm'] = df5['PubMed Record Count'] / df5['PubMed Record Count'].max()
df5['AMOS methods count_norm'] = df5['AMOS methods count'] / df5['AMOS methods count'].max()
df5['AMOS fact sheets count_norm'] = df5['AMOS fact sheets count'] / df5['AMOS fact sheets count'].max()
df5['AMOS spectra count_norm'] = df5['AMOS spectra count'] / df5['AMOS spectra count'].max()
df5['Presence in water lists count_norm'] = df5['Presence in water lists count'] / df5['Presence in water lists count'].max()

# Set up the figure and axes
fig, ax = plt.subplots()

# Width of each bar
bar_width = 0.08

chemicals = df5['DTXSID'].tolist()

truncate = 132

#chemicals_reverse = chemicals[::-1]
chemicals_reverse = chemicals[:truncate]
chemicals_reverse = chemicals_reverse[::-1]

metadata_fields = ['Sources_norm', 'Patents_norm', 'Articles_norm', 'PubMed Record Count_norm', 'AMOS methods count_norm', 'AMOS fact sheets count_norm', 
                   'AMOS spectra count_norm', 'Presence in water lists count_norm']

# # Create the bars
# for i, group in enumerate(metadata_fields):
#     ax.bar(x_pos + i * bar_width, values[:, i], width=bar_width, label=group)


# Create horizontal bars
for i, group in enumerate(metadata_fields):
    y_pos = np.arange(len(chemicals_reverse)) + i * bar_width
    temp_data = df5[group][:truncate]
    temp_data = temp_data[::-1]
    ax.barh(y_pos, temp_data, height=bar_width, label=group)

test = np.arange(len(chemicals_reverse))

# Set y-axis labels and title
ax.set_yticks(np.arange(len(chemicals_reverse)) + bar_width * (len(metadata_fields) - 1) / 2)
ax.set_yticklabels(chemicals_reverse)
#ax.set_xlabel('Values')
#ax.set_title('Multi-Bar Chart with Horizontal Bars')

ax.tick_params(axis='x', labelsize=3)

# Find index of SSM chemical
temp_boolean = df5[['SSM chemical']][:truncate]

temp_boolean = temp_boolean[::-1]
temp_boolean.reset_index(inplace=True, drop=True)
#SSM_chemical_index = df5[df5['SSM chemical']].index[0]
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]

# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)

# Increase the figure size
fig.set_size_inches(12, 40) 
#fig.set_size_inches(6, 10) 
# Add a legend, show plot
plt.legend()
plt.show()

'''




'''

# Chemical metadata plots hereherehere

truncate = 0
select_chem = 'atenolol'
#select_chem = 'theophylline'
select_SSM_structure_only = False

##################################################################################################################
# Leg 1 stacked bar chart - USE THIS ONE
##################################################################################################################
# Read in MS1 results for just the SSM chemicals
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df_merge = pd.read_csv('WW2DW MS1 results (chemical and hazard) for SSM chemicals with search results included.csv')

# Start off with just grabbing a single feature
df5 = df_merge[df_merge['Chemical Name'] == select_chem]
# Remove feature duplicates if SSM chemical in more than one mode
df5 = df5.sort_values(by='Feature ID')
df5 = df5.drop_duplicates(subset='DTXSID', keep='first')


if select_SSM_structure_only:
    temp_df = df5[df5['SSM chemical'] == True]
    SSM_DTXCID = temp_df['DTXCID_INDIVIDUAL_COMPONENT'].tolist()[0]
    #print(SSM_DTXCID)
    df5 = df5[df5['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID]


# Normalize all columns to graph
df5['Sources_norm'] = df5['Sources'] / df5['Sources'].max()
df5['Patents_norm'] = df5['Patents'] / df5['Patents'].max()
df5['Articles_norm'] = df5['Articles'] / df5['Articles'].max()
df5['PubMed Record Count_norm'] = df5['PubMed Record Count'] / df5['PubMed Record Count'].max()
df5['AMOS methods count_norm'] = df5['AMOS methods count'] / df5['AMOS methods count'].max()
df5['AMOS fact sheets count_norm'] = df5['AMOS fact sheets count'] / df5['AMOS fact sheets count'].max()
df5['AMOS spectra count_norm'] = df5['AMOS spectra count'] / df5['AMOS spectra count'].max()
df5['Presence in water lists count_norm'] = df5['Presence in water lists count'] / df5['Presence in water lists count'].max()

df5.fillna(0, inplace=True)

df5['Total_norm'] = df5['Sources_norm'] + df5['Patents_norm'] + df5['Articles_norm'] + df5['PubMed Record Count_norm'] + df5['AMOS methods count_norm'] + df5['AMOS fact sheets count_norm']  +  df5['AMOS spectra count_norm'] + df5['Presence in water lists count_norm']

df5 = df5.sort_values(by='Total_norm', ascending=True)

df5_plot = df5[['PREFERRED_NAME', 'Sources_norm', 'Patents_norm', 'Articles_norm', 'PubMed Record Count_norm', 'AMOS methods count_norm', 'AMOS fact sheets count_norm', 
                   'AMOS spectra count_norm', 'Presence in water lists count_norm', 'SSM chemical']]

#df5_plot['PREFERRED_NAME'] = df5_plot['PREFERRED_NAME'].str[0:30]

if truncate > 0:
    df5_plot = df5_plot.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
else:
    fig, ax = plt.subplots(figsize=(12, 40))

# Find index of SSM chemical
temp_boolean = df5_plot[['SSM chemical']]
temp_boolean.reset_index(inplace=True, drop=True)
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]
df5_plot = df5_plot.drop('SSM chemical', axis=1)

# Initialize bottom for stacking
bottom = [0] * len(df5_plot)

# Plot each column
for col in df5_plot.columns[1:]:
    ax.barh(df5_plot['PREFERRED_NAME'], df5_plot[col], left=bottom, label=col)
    bottom += df5_plot[col]

# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)

plt.legend()
plt.ylabel('Candidate')
plt.xlabel('Total Metadata Score')

title_string = 'Metadata Results for SSM Chemical: ' + select_chem
plt.title(title_string)
plt.show()

##################################################################################################################
# Leg 1 stacked bar chart - Collapse on DTXCID
# Sum together metadata for all substances associated with a structure/DTXCID
##################################################################################################################
metadata_columns = ['Sources', 'Patents', 'Articles', 'PubMed Record Count', 'AMOS methods count', 'AMOS fact sheets count', 'AMOS spectra count', 'Presence in water lists count']

use_transparency = False

# Identify the DTXCID associated with the SSM chemical
df5_true = df5[df5['SSM chemical'] == True]
SSM_DTXCID = df5_true['DTXCID_INDIVIDUAL_COMPONENT'].tolist()[0]

#df6 = df5.groupby('DTXCID_INDIVIDUAL_COMPONENT').agg({col: 'sum' for col in df5.columns if col.endswith('_norm')})
df6 = df5.groupby('DTXCID_INDIVIDUAL_COMPONENT').agg({col: 'sum' for col in metadata_columns})

df6 = df6.reset_index()

#df5_select = ['Feature ID', 'Chemical Name', 'DTXCID_INDIVIDUAL_COMPONENT'] + [col for col in metadata_columns]
#debug_5 = df5[df5_select]


df6['Sources_norm'] = df6['Sources'] / df6['Sources'].max()
df6['Patents_norm'] = df6['Patents'] / df6['Patents'].max()
df6['Articles_norm'] = df6['Articles'] / df6['Articles'].max()
df6['PubMed Record Count_norm'] = df6['PubMed Record Count'] / df6['PubMed Record Count'].max()
df6['AMOS methods count_norm'] = df6['AMOS methods count'] / df6['AMOS methods count'].max()
df6['AMOS fact sheets count_norm'] = df6['AMOS fact sheets count'] / df6['AMOS fact sheets count'].max()
df6['AMOS spectra count_norm'] = df6['AMOS spectra count'] / df6['AMOS spectra count'].max()
df6['Presence in water lists count_norm'] = df6['Presence in water lists count'] / df6['Presence in water lists count'].max()

df6['Total_norm'] = df6['Sources_norm'] + df6['Patents_norm'] + df6['Articles_norm'] + df6['PubMed Record Count_norm'] + df6['AMOS methods count_norm'] + df6['AMOS fact sheets count_norm'] + df6['AMOS spectra count_norm'] + df6['Presence in water lists count_norm']

# Do some filtering on candidates with only 1 source
df6 = df6[df6['Total_norm'] != 0]
df6 = df6[df6['Sources'] > 1]


df6 = df6.sort_values(by='Total_norm', ascending=True)

df6_plot = df6[['DTXCID_INDIVIDUAL_COMPONENT', 'Sources_norm', 'Patents_norm', 'Articles_norm', 'PubMed Record Count_norm', 'AMOS methods count_norm', 'AMOS fact sheets count_norm', 
                   'AMOS spectra count_norm', 'Presence in water lists count_norm']]

#df5_plot['PREFERRED_NAME'] = df5_plot['PREFERRED_NAME'].str[0:30]

if truncate > 0:
    df6 = df6.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
    SSM_font_size = 12
    label_sizes = 12
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    SSM_font_size = 22
    label_sizes = 22
    ax.tick_params(axis='y', labelsize=label_sizes - 4)
    ax.tick_params(axis='x', labelsize=label_sizes + 12)

# Initialize bottom for stacking
bottom = [0] * len(df6_plot)

# Loop through each row (each group of stacked bars)
for i, label in enumerate(df6_plot['DTXCID_INDIVIDUAL_COMPONENT']):
    # Set different alpha if the label matches
    #alpha_value = 1.0 if label == SSM_DTXCID else 0.4
    if use_transparency == True:
        alpha_value = 1.0 if label == SSM_DTXCID else 0.4
    else:
        alpha_value = 1
    # Reset bottom for each row
    bottom_value = 0  

    # Loop through each column for stacking
    for col in df6_plot.columns[1:]:
        ax.barh(label, df6_plot[col].iloc[i], left=bottom_value, label=col if i == 0 else "", alpha=alpha_value)
        bottom_value += df6_plot[col].iloc[i]  # Update bottom for stacking


plt.legend()
plt.ylabel('Candidate', fontsize=label_sizes+12)
plt.xlabel('Total Metadata Score', fontsize=label_sizes+12)

# Set size and color of DTXCID corresponding to SSM chemical
for label in ax.get_yticklabels():
    if label.get_text() == SSM_DTXCID:
        label.set_fontsize(SSM_font_size)  # Set font size to 16
        label.set_color('blue')

legend = ax.legend(fontsize=label_sizes + 4, loc='lower right')

# Iterate through the legend handles (patches)
for handle in legend.legendHandles:
    # Set the alpha value for each handle
    handle.set_alpha(1)  # Adjust the value between 0 (fully transparent) and 1 (fully opaque)


plt.ylabel('Candidate', fontsize=label_sizes+12)
plt.xlabel('Metadata Score', fontsize=label_sizes+12)

title_string = 'Metadata Results for SSM Chemical: ' + select_chem
plt.title(title_string, fontsize=label_sizes + 12)
plt.show()

'''
# # Plot each column
# for col in df6_plot.columns[1:]:
#     ax.barh(df6_plot['DTXCID_INDIVIDUAL_COMPONENT'], df6_plot[col], left=bottom, label=col)
#     bottom += df6_plot[col]


# plt.legend()
# plt.ylabel('Candidate', fontsize=label_sizes+12)
# plt.xlabel('Total Metadata Score', fontsize=label_sizes+12)

# # Set size and color of DTXCID corresponding to SSM chemical
# for label in ax.get_yticklabels():
#     if label.get_text() == SSM_DTXCID:
#         label.set_fontsize(SSM_font_size)  # Set font size to 16
#         label.set_color('green')

# ax.legend(fontsize=label_sizes + 4, loc='lower right')

# title_string = 'Metadata Results for SSM Chemical: ' + select_chem
# plt.title(title_string, fontsize=label_sizes + 12)
# plt.show()



'''
# Read in MS1 results for just the SSM chemicals
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df_merge = pd.read_csv('WW2DW MS1 results (chemical and hazard) for SSM chemicals with search results included.csv')

# Start off with just grabbing a single feature
df5 = df_merge[df_merge['Chemical Name'] == 'acetaminophen']
# Normalize all columns to graph
df5['Sources_norm'] = df5['Sources'] / df5['Sources'].max()
df5['Patents_norm'] = df5['Patents'] / df5['Patents'].max()
df5['Articles_norm'] = df5['Articles'] / df5['Articles'].max()
df5['PubMed Record Count_norm'] = df5['PubMed Record Count'] / df5['PubMed Record Count'].max()
df5['AMOS methods count_norm'] = df5['AMOS methods count'] / df5['AMOS methods count'].max()
df5['AMOS fact sheets count_norm'] = df5['AMOS fact sheets count'] / df5['AMOS fact sheets count'].max()
df5['AMOS spectra count_norm'] = df5['AMOS spectra count'] / df5['AMOS spectra count'].max()
df5['Presence in water lists count_norm'] = df5['Presence in water lists count'] / df5['Presence in water lists count'].max()

df5['Total_norm'] =df5['Sources_norm'] + df5['Patents_norm'] + df5['Articles_norm'] + df5['PubMed Record Count_norm'] + df5['AMOS methods count_norm'] + df5['AMOS fact sheets count_norm']  +  df5['AMOS spectra count_norm'] + df5['Presence in water lists count_norm']

#chemicals = df5['DTXSID'].tolist()

df5_plot = df5[['DTXSID', 'Sources_norm', 'Patents_norm', 'Articles_norm', 'PubMed Record Count_norm', 'AMOS methods count_norm', 'AMOS fact sheets count_norm', 
                   'AMOS spectra count_norm', 'Presence in water lists count_norm']]

truncate = 25
df5_plot = df5_plot[:truncate]
df5_plot = df5_plot.iloc[::-1]

df5_plot.fillna(0, inplace=True)

# Plotting
fig, ax = plt.subplots()

# Initialize bottom for stacking
bottom = [0] * len(df5_plot)

# Plot each column
for col in df5_plot.columns[1:]:
    print(col)
    print(bottom)
    ax.barh(df5_plot['DTXSID'], df5_plot[col], left=bottom, label=col)
    bottom += df5_plot[col]

# Find index of SSM chemical
temp_boolean = df5[['SSM chemical']][:truncate]
temp_boolean = temp_boolean[::-1]
temp_boolean.reset_index(inplace=True, drop=True)
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]
# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)

# Increase the figure size
#fig.set_size_inches(12, 40) 
fig.set_size_inches(6, 10) 
plt.legend()
plt.ylabel('Candidate')
plt.show()
'''


'''
##################################################################################################################
# Convert hazard data into numeric values
##################################################################################################################
# Read in MS1 results with hazard data
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df_merge = pd.read_csv('WW2DW MS1 results (chemical and hazard) for SSM chemicals with search results included.csv')
keep_strings = ['DTXSID', '_score', '_authority']
# Filter columns based on the list of strings
df_hazard = df_merge.filter(regex='|'.join(keep_strings))
df_hazard = df_hazard.drop('Tracer DTXSID', axis=1)
df_hazard_b = df_hazard.replace('ND', np.nan)
df_hazard_b = df_hazard_b.replace('Authoritative', 3)
df_hazard_b = df_hazard_b.replace('Authoritative', 3)
df_hazard_b = df_hazard_b.replace('Screening', 2)
df_hazard_b = df_hazard_b.replace('QSAR Model', 1)
df_hazard_b = df_hazard_b.replace('VH', 4)
df_hazard_b = df_hazard_b.replace('H', 3)
df_hazard_b = df_hazard_b.replace('M', 2)
df_hazard_b = df_hazard_b.replace('L', 1)

for col in df_hazard_b.columns:
    if col.endswith("_authority"):
        score_col = col.replace("_authority", "_score")
        if score_col in df_hazard_b.columns:
            df_hazard_b.loc[df_hazard_b[score_col] == "I", col] = np.nan
            
df_hazard_b = df_hazard_b.replace('I', np.nan)

score_cols = [col for col in df_hazard_b.columns if col.endswith("_score")]
authority_cols = [col for col in df_hazard_b.columns if col.endswith("_authority")]

df_hazard_b['Average Hazard Score'] = df_hazard_b[score_cols].mean(axis=1)
df_hazard_b['Average Quality Score'] = df_hazard_b[authority_cols].mean(axis=1)
df_hazard_b['Quality-Adjusted Hazard Score'] = df_hazard_b['Average Hazard Score'] * df_hazard_b['Average Quality Score']
df_hazard_b['Number of end points with data available'] = df_hazard_b[[col for col in df_hazard_b.columns if col.endswith("_score")]].notna().sum(axis=1)
df_hazard_b['Completeness Score'] = df_hazard_b['Number of end points with data available'] / len(score_cols)

df_hazard_b = df_hazard_b.drop_duplicates()

# Merge numeric values back onto original results
df_merge_b = df_merge.drop(columns=score_cols)
df_merge_b = df_merge_b.drop(columns=authority_cols)

df_merge_b = pd.merge(df_merge_b, df_hazard_b, on='DTXSID', how='left')

# Export the converted hazard results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df_merge_b.to_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv', sep=',', encoding='utf-8', index=False)
'''

'''
##################################################################################################################
# Hazard results plotting - Multi-bar chart
##################################################################################################################
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')



# Start off with just grabbing a single feature
df1 = df1[df1['Chemical Name'] == 'atenolol']
df1 = df1.sort_values(by='Quality-Adjusted Hazard Score', ascending=False)

# Set up the figure and axes
fig, ax = plt.subplots()

# Width of each bar
bar_width = 0.4

chemicals = df1['DTXSID'].tolist()

truncate = 25

#chemicals_reverse = chemicals[::-1]
chemicals_reverse = chemicals[:truncate]
chemicals_reverse = chemicals_reverse[::-1]

metadata_fields = ['Quality-Adjusted Hazard Score', 'Completeness Score']


# Try normalizing completeness score to highest QA-hazard score
max_QA = df1['Quality-Adjusted Hazard Score'].max()
max_completeness = df1['Completeness Score'].max()
#df1['Completeness Score'] = df1['Completeness Score'] * max_QA / max_completeness
df1['Completeness Score'] = df1['Completeness Score'] * max_QA


# # Create the bars
# for i, group in enumerate(metadata_fields):
#     ax.bar(x_pos + i * bar_width, values[:, i], width=bar_width, label=group)


# Create horizontal bars
for i, group in enumerate(metadata_fields):
    y_pos = np.arange(len(chemicals_reverse)) + i * bar_width
    temp_data = df1[group][:truncate]
    print(df1[group])
    print(temp_data)
    temp_data = temp_data[::-1]
    ax.barh(y_pos, temp_data, height=bar_width, label=group)

test = np.arange(len(chemicals_reverse))

# Set y-axis labels and title
ax.set_yticks(np.arange(len(chemicals_reverse)) + bar_width * (len(metadata_fields) - 1) / 2)
ax.set_yticklabels(chemicals_reverse)

ax.tick_params(axis='x', labelsize=3)

# Find index of SSM chemical
temp_boolean = df1[['SSM chemical']][:truncate]

temp_boolean = temp_boolean[::-1]
temp_boolean.reset_index(inplace=True, drop=True)
#SSM_chemical_index = df5[df5['SSM chemical']].index[0]
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]

# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)

# Increase the figure size
#fig.set_size_inches(12, 40) 
fig.set_size_inches(6, 10) 
# Add a legend, show plot
plt.legend()
plt.show()
'''

'''
##################################################################################################################
# Hazard results plotting - Bar chart slashes based on completeness
##################################################################################################################
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

# Start off with just grabbing a single feature
df2 = df1[df1['Chemical Name'] == 'atenolol']
df2['Quality-Adjusted Hazard Score'] = df2['Quality-Adjusted Hazard Score'].fillna(0)
df2 = df2.sort_values(by='Quality-Adjusted Hazard Score', ascending=True)

#truncate = 25
#df2 = df2.tail(truncate)

#plt.figure(figsize=(12, 40))
#plt.yticks(fontsize=20)

# Generate list of linewidths and step values based on completeness score
completeness_list = df2['Completeness Score'].tolist()
line_widths = []
step_list = []
for score in completeness_list:
    if score < 0.2:
        line_widths.append(2)
        step_list.append(0.2)
    elif score < 0.4:
        line_widths.append(1)
        step_list.append(0.3)
    elif score < 0.6:
        line_widths.append(0.5)
        step_list.append(0.35)
    elif score < 0.8:
        line_widths.append(0.1)
        step_list.append(0.4)
    else:
        line_widths.append(0)
        step_list.append(0.1)        
        
        
#fig, ax = plt.subplots(figsize=(12, 40))
fig, ax = plt.subplots(figsize=(6, 10))
bars = ax.barh(df2['DTXSID'], df2['Quality-Adjusted Hazard Score'], color='skyblue', edgecolor='white')

#plt.barh(df2['PREFERRED_NAME'], df2['Quality-Adjusted Hazard Score'], color='skyblue')
#plt.barh(df2['DTXSID'], df2['Quality-Adjusted Hazard Score'], color='skyblue', edgecolor='black')

for bar, line_width, step_width in zip(bars, line_widths, step_list):
    x, y, width, height = bar.get_x(), bar.get_y(), bar.get_width(), bar.get_height()
    step=step_width
    rounded_range = round(width / step) + 1
    #for i in range(int(width / step) + 1):
    for i in range(int(rounded_range)):
        ax.plot([x + i * step, x + i * step + height], [y, y + height], color="white", linewidth=line_width)


# Find index of SSM chemical
temp_boolean = df2[['SSM chemical']]

#temp_boolean = temp_boolean[::-1]
temp_boolean.reset_index(inplace=True, drop=True)
#SSM_chemical_index = df5[df5['SSM chemical']].index[0]
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]
# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)

plt.show()
'''
'''
##################################################################################################################
# Create blue gradient rectangle for hazard plot legend
##################################################################################################################
fig, ax = plt.subplots(figsize=(6,3)) 
gradient = np.linspace(0.1, 1, 1000).reshape(1, -1)

blue_color = np.zeros((1, 1000, 4))
blue_color[:, :, 2] = 1
blue_color[:, :, 3] = gradient

ax.imshow(blue_color, extent=[0, 5, 0, 2], aspect = 'auto')

ax.set_xticks([])
ax.set_yticks([])
ax.set_frame_on(False)

plt.show()
'''


'''

# Hazard data plotting hereherehere
truncate = 0
#select_chem = 'atenolol'
select_chem = 'theophylline'
select_SSM_structure_only = False

##################################################################################################################
# Hazard results plotting - Bar chart transparency based on completeness - GO WITH THIS
##################################################################################################################

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

# Start off with just grabbing a single feature
df2 = df1[df1['Chemical Name'] == select_chem]

# Remove feature duplicates if SSM chemical in more than one mode
df2 = df2.sort_values(by='Feature ID')
df2 = df2.drop_duplicates(subset='DTXSID', keep='first')

if select_SSM_structure_only:
    temp_df = df2[df2['SSM chemical'] == True]
    SSM_DTXCID = temp_df['DTXCID_INDIVIDUAL_COMPONENT'].tolist()[0]
    #print(SSM_DTXCID)
    df2 = df2[df2['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID]




debug = df2.copy()
debug = debug[['PREFERRED_NAME', 'Quality-Adjusted Hazard Score']]

df2['Quality-Adjusted Hazard Score'] = df2['Quality-Adjusted Hazard Score'].fillna(0)
df2 = df2.sort_values(by='Quality-Adjusted Hazard Score', ascending=True)

if select_SSM_structure_only:
    fig, ax = plt.subplots(figsize=(3, 10))
elif truncate > 0:
    df2 = df2.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    

new_min, new_max = 0.01, 1
df2['Transparency'] = new_min + df2['Completeness Score'] * (new_max - new_min)

for index, row in df2.iterrows():
    ax.barh(row['PREFERRED_NAME'], row['Quality-Adjusted Hazard Score'], color="blue", alpha=row['Transparency'])
    

# Find index of SSM chemical
temp_boolean = df2[['SSM chemical']]
temp_boolean.reset_index(inplace=True, drop=True)
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]

# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)

# labels = [f"{cat} (Alpha: {alpha:.1f})" for cat, alpha in zip(df2['DTXSID'], df2['Transparency'])]
plt.xlabel('Quality-Adjusted Hazard Score')
plt.ylabel('Candidate')
title_string = 'Hazard Results for SSM Chemical: ' + df2['Chemical Name'].tolist()[0]
plt.title(title_string)

plt.show()

debug_a = df2[['DTXCID_INDIVIDUAL_COMPONENT', 'Quality-Adjusted Hazard Score', 'Completeness Score']]


##################################################################################################################
# Hazard results plotting - Bar chart transparency based on completeness
# Collapse on Structure/DTXCID try 3 - Take max QAH value amongst all substances for a structure
##################################################################################################################
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df6 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

# Start off with just grabbing a single feature
df7 = df6[df6['Chemical Name'] == select_chem]

# Identify the DTXCID associated with the SSM chemical
df7_true = df7[df7['SSM chemical'] == True]
SSM_DTXCID = df7_true['DTXCID_INDIVIDUAL_COMPONENT'].tolist()[0]

chemical_name = df7['Chemical Name'].tolist()[0]

# Remove feature duplicates if SSM chemical in more than one mode
df7 = df7.sort_values(by='Feature ID')
df7 = df7.drop_duplicates(subset='DTXSID', keep='first')

df7['Quality-Adjusted Hazard Score'] = df7['Quality-Adjusted Hazard Score'].fillna(0)
#df7 = df7.sort_values(by='Quality-Adjusted Hazard Score', ascending=True)
df7 = df7.sort_values(by=['Quality-Adjusted Hazard Score', 'Completeness Score'], ascending=True)

# Remove structures with no score
df7 = df7[df7['Quality-Adjusted Hazard Score'] > 0]


# De-duplicate on structure and keep highest QAH score
df8 = df7.drop_duplicates(subset='DTXCID_INDIVIDUAL_COMPONENT', keep='last')
debug_c = df7.drop_duplicates(subset='DTXCID_INDIVIDUAL_COMPONENT', keep='last')


if truncate > 0:
    df8 = df8.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
    SSM_font_size = 12
    label_sizes = 12
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    SSM_font_size = 22
    label_sizes = 22
    ax.tick_params(axis='y', labelsize=label_sizes - 4)
    ax.tick_params(axis='x', labelsize=label_sizes + 12)
    
new_min, new_max = 0.01, 1
df8['Transparency'] = new_min + df8['Completeness Score'] * (new_max - new_min)

use_select_color = False


for index, row in df8.iterrows():
    #bar_color = "magenta" if row['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID else "blue"
    if use_select_color == True:
        bar_color = "blue" if row['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID else "red"
        ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color=bar_color, alpha=row['Transparency'])
        bar_outline = 10 if row['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID else 0
    else:
        bar_color = "blue"
        ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color=bar_color, alpha=row['Transparency'])
        bar_outline = 0
    #ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color='blue', edgecolor='red', linewidth=bar_outline, alpha=row['Transparency'])  
    #ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color='red', edgecolor='blue', linewidth=bar_outline, alpha=row['Transparency'])  

# Set size and color of DTXCID corresponding to SSM chemical
for label in ax.get_yticklabels():
    if label.get_text() == SSM_DTXCID:
        label.set_fontsize(SSM_font_size)
        label.set_color('blue')



plt.xlabel('Quality-Adjusted Hazard Score', fontsize=label_sizes+12)
plt.ylabel('Candidate', fontsize=label_sizes+12)
title_string = 'Hazard Results for SSM Chemical (approach 1): ' + chemical_name
plt.title(title_string, fontsize=label_sizes+12)

plt.show()
'''

'''
##################################################################################################################
# SSM feature results plotting hereherehere
# Collapse on structure:
#   Calculate Summed metadata
#   Calculate max QAH
# Count structures
# Identify rank of true SSM
# Plot line charts which show number of candidates and the true SSM rank for each category 
##################################################################################################################

# Read in MS1 results for just the SSM chemicals
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

# Remove feature duplicates if SSM chemical in more than one mode
df2 = df1.sort_values(by='Feature ID')
df2 = df2.drop_duplicates(subset='DTXSID', keep='first')

col_list = df2.columns.tolist()

df2 = df2[['Feature ID', 'Chemical Name', 'Tracer DTXSID', 'DTXCID_INDIVIDUAL_COMPONENT', 'DTXSID', 'PREFERRED_NAME', 'SSM chemical', 
           'Sources', 'Patents', 'Articles', 'PubMed Record Count', 'AMOS methods count', 'AMOS fact sheets count', 'AMOS spectra count', 
           'Presence in water lists count', 'Quality-Adjusted Hazard Score', 'Completeness Score']]

df2 = df2.fillna(0)

# Create a column of the SSM chemical structure True/False so we can eventually collapse into structure
df2['SSM structure'] = df2.groupby(['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])['SSM chemical'].transform('any')

# Now collapse hazard/metadata on structure
# HAZARD: Create a column that has the maximum QAH score for a given structure within a feature
df2['Max_structure_QAH'] = df2.groupby(['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])['Quality-Adjusted Hazard Score'].transform('max')

# METADATA: sum across substances for structures
metadata_columns = ['Sources', 'Patents', 'Articles', 'PubMed Record Count', 'AMOS methods count', 'AMOS fact sheets count', 'AMOS spectra count', 'Presence in water lists count']

# Within each structure within each feature, sum the substances
df3 = df2.groupby(['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT']).agg({col: 'sum' for col in metadata_columns})
df3 = df3.reset_index()

# Within each feature, normalize each metadata value to the maximum value
for col in metadata_columns:
    df3[f'{col}_norm'] = df3.groupby('Feature ID')[col].transform(lambda x: x / x.max())

df3 = df3.fillna(0)
df3['Total_norm'] = df3['Sources_norm'] + df3['Patents_norm'] + df3['Articles_norm'] + df3['PubMed Record Count_norm'] + df3['AMOS methods count_norm'] + df3['AMOS fact sheets count_norm'] + df3['AMOS spectra count_norm'] + df3['Presence in water lists count_norm']


# Get count of structures for each group
df3['structure_count'] = df3.groupby('Feature ID')['DTXCID_INDIVIDUAL_COMPONENT'].transform('count')

# Combine and grab only the needed columns
#df4 = pd.merge(df2[['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT', 'SSM structure', 'Max_structure_QAH']], df3[['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT', 'structure_count', 'Total_norm']], how='left', on=['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])
df4 = pd.merge(df3[['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT', 'structure_count', 'Total_norm']], df2[['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT', 'SSM structure', 'Max_structure_QAH']], how='left', on=['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])
df4 = df4.drop_duplicates()

df4.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 

# Bring in MS2
# Read in the MS2 results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/neg')
df_ms2_neg = pd.read_csv('20210131_WW2DW_MS2_neg_CFMID_results_neg.csv')
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/pos')
df_ms2_pos = pd.read_csv('20250131_WW2DW_MS2_pos_CFMID_results_pos.csv')

# Need to add a RT tolerance for merging... or fix merge workflow ac 2/18/2025
df_ms2_neg.rename(columns={'SUM_SCORE':'MS2 neg'}, inplace=True) 
df_ms2_pos.rename(columns={'SUM_SCORE':'MS2 pos'}, inplace=True) 

df5 = pd.merge(df4, df_ms2_neg[['DTXCID', 'MS2 neg']], how='left', on='DTXCID')
df5 = pd.merge(df5, df_ms2_pos[['DTXCID', 'MS2 pos']], how='left', on='DTXCID')

#debug = df5[df5['DTXCID'] == 'DTXCID20208869']
#debug = df5[df5['Feature ID'] == 2901]

# De-duplicate merges where both modes data is present
df5 = df5.sort_values('MS2 neg', ascending=False).drop_duplicates(['Feature ID', 'DTXCID', 'MS2 neg']).sort_index()
df5 = df5.sort_values('MS2 pos', ascending=False).drop_duplicates(['Feature ID', 'DTXCID', 'MS2 pos']).sort_index()

rank_columns = ['Total_norm', 'Max_structure_QAH', 'MS2 pos', 'MS2 neg']

# Generate rank columns
for col in rank_columns:
    df5[f'{col}_rank'] = df5.groupby('Feature ID')[col].rank(method='dense', ascending=False)


# Grab just the SSM structures
df6 = df5[df5['SSM structure'] == True]

# Sort by structure counts
df6 = df6.sort_values(by='structure_count', ascending=True).reset_index(drop=True)
df6.reset_index()


# Plot all together
plt.figure(figsize=(20, 14))
# Draw horizontal lines from x=0 to x=LineValue for each Category
plt.hlines(y=df6['DTXCID'], xmin=0, xmax=df6['structure_count'], color='skyblue', linewidth=5)
# List of marker columns to add to the lines
marker_columns = ['Total_norm_rank', 'Max_structure_QAH_rank', 'MS2 pos_rank', 'MS2 neg_rank']
marker_colors  = ['green', 'red', 'orange', 'purple']  # Distinct colors for each marker column
marker_shapes = ['|', '|', '2', '1']
# Plot markers from each of the four additional columns
for marker, color, shape in zip(marker_columns, marker_colors, marker_shapes):
    plt.scatter(df6[marker], df6['DTXCID'], color=color, label=marker, marker=shape, zorder=3, s=100)
# Add labels, title, and legend
plt.xlabel("Number of Candidate structures")
plt.ylabel("SSM Structure")
plt.title("SSM Chemical Structure Rank Results: All")
plt.legend(loc='lower right')
plt.show()

# Plot just metadata
plt.figure(figsize=(20, 14))
# Draw horizontal lines from x=0 to x=LineValue for each Category
plt.hlines(y=df6['DTXCID'], xmin=0, xmax=df6['structure_count'], color='skyblue', linewidth=5)
# List of marker columns to add to the lines
marker_columns = ['Total_norm_rank']
marker_colors  = ['green']  # Distinct colors for each marker column
marker_shapes = ['|']
# Plot markers from each of the four additional columns
for marker, color, shape in zip(marker_columns, marker_colors, marker_shapes):
    plt.scatter(df6[marker], df6['DTXCID'], color=color, label=marker, marker=shape, zorder=3, s=100)
# Add labels, title, and legend
plt.xlabel("Number of Candidate structures")
plt.ylabel("SSM Structure")
plt.title("SSM Chemical Structure Rank Results: Metadata")
plt.legend(loc='lower right')
plt.show()

# Plot just hazard
plt.figure(figsize=(20, 14))
# Draw horizontal lines from x=0 to x=LineValue for each Category
plt.hlines(y=df6['DTXCID'], xmin=0, xmax=df6['structure_count'], color='skyblue', linewidth=5)
# List of marker columns to add to the lines
marker_columns = ['Max_structure_QAH_rank']
marker_colors  = ['red']  # Distinct colors for each marker column
marker_shapes = ['|']
# Plot markers from each of the four additional columns
for marker, color, shape in zip(marker_columns, marker_colors, marker_shapes):
    plt.scatter(df6[marker], df6['DTXCID'], color=color, label=marker, marker=shape, zorder=3, s=100)
# Add labels, title, and legend
plt.xlabel("Number of Candidate structures")
plt.ylabel("SSM Structure")
plt.title("SSM Chemical Structure Rank Results: Hazard")
plt.legend(loc='lower right')
plt.show()


# Plot all together
plt.figure(figsize=(20, 14))
# Draw horizontal lines from x=0 to x=LineValue for each Category
plt.hlines(y=df6['DTXCID'], xmin=0, xmax=df6['structure_count'], color='skyblue', linewidth=5)
# List of marker columns to add to the lines
marker_columns = ['MS2 pos_rank', 'MS2 neg_rank']
marker_colors  = ['orange', 'purple']  # Distinct colors for each marker column
marker_shapes = ['2', '1']
# Plot markers from each of the four additional columns
for marker, color, shape in zip(marker_columns, marker_colors, marker_shapes):
    plt.scatter(df6[marker], df6['DTXCID'], color=color, label=marker, marker=shape, zorder=3, s=100)
# Add labels, title, and legend
plt.xlabel("Number of Candidate structures")
plt.ylabel("SSM Structure")
plt.title("SSM Chemical Structure Rank Results: MS2")
plt.legend(loc='lower right')
plt.show()
'''

'''
# Normalize all columns to graph
df2['Sources_norm'] = df2['Sources'] / df2['Sources'].max()
df2['Patents_norm'] = df2['Patents'] / df2['Patents'].max()
df2['Articles_norm'] = df2['Articles'] / df2['Articles'].max()
df2['PubMed Record Count_norm'] = df2['PubMed Record Count'] / df2['PubMed Record Count'].max()
df2['AMOS methods count_norm'] = df2['AMOS methods count'] / df2['AMOS methods count'].max()
df2['AMOS fact sheets count_norm'] = df2['AMOS fact sheets count'] / df2['AMOS fact sheets count'].max()
df2['AMOS spectra count_norm'] = df2['AMOS spectra count'] / df2['AMOS spectra count'].max()
df2['Presence in water lists count_norm'] = df2['Presence in water lists count'] / df2['Presence in water lists count'].max()

df2['Total_norm'] = df2['Sources_norm'] + df2['Patents_norm'] + df2['Articles_norm'] + df2['PubMed Record Count_norm'] + df2['AMOS methods count_norm'] + df2['AMOS fact sheets count_norm']  +  df2['AMOS spectra count_norm'] + df2['Presence in water lists count_norm']
'''

'''
df2 = df2.sort_values(by=['Feature ID', 'SSM chemical'], ascending=[True, False])

df2_structure_only = df2.drop_duplicates(subset='DTXCID_INDIVIDUAL_COMPONENT', keep='first')





# Group by 'Category' and get counts using size()
counts_size = df2_structure_only.groupby('Chemical Name').size()
print("Counts using size():\n", counts_size)

# Group by 'Category' and get counts using count()
counts_count = df2_structure_only.groupby('Feature ID').count()
print("\nCounts using count():\n", counts_count)

counts_size = counts_size.sort_values()

counts_size.plot(kind='barh', color='blue', lw=1, alpha=0.5, figsize=(10, 16))
plt.title("Number of Candidate Structures for SSM Chemical Feature")
#plt.figure(figsize=(10, 16))

plt.show()
'''
# Collapse on structure for metadata sum data

# df2['Sources_norm'] = df2['Sources'] / df2['Sources'].max()
# df2['Patents_norm'] = df2['Patents'] / df2['Patents'].max()
# df2['Articles_norm'] = df2['Articles'] / df2['Articles'].max()
# df2['PubMed Record Count_norm'] = df2['PubMed Record Count'] / df2['PubMed Record Count'].max()
# df2['AMOS methods count_norm'] = df2['AMOS methods count'] / df2['AMOS methods count'].max()
# df2['AMOS fact sheets count_norm'] = df2['AMOS fact sheets count'] / df2['AMOS fact sheets count'].max()
# df2['AMOS spectra count_norm'] = df2['AMOS spectra count'] / df2['AMOS spectra count'].max()
# df2['Presence in water lists count_norm'] = df2['Presence in water lists count'] / df2['Presence in water lists count'].max()

# df2['Total_norm'] = df2['Sources_norm'] + df2['Patents_norm'] + df2['Articles_norm'] + df2['PubMed Record Count_norm'] + df2['AMOS methods count_norm'] + df2['AMOS fact sheets count_norm'] + df2['AMOS spectra count_norm'] + df2['Presence in water lists count_norm']



'''
##################################################################################################################
# Hazard results plotting - Bar chart transparency based on completeness
# Collapse on Structure/DTXCID try 2 - Take max value within each endpoint, recalculate QAH
# Group didn't like this
##################################################################################################################
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1_b = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

df2_b = df1_b[df1_b['Chemical Name'] == select_chem]

chemical_name = df2_b['Chemical Name'].tolist()[0]

# Remove feature duplicates if SSM chemical in more than one mode
df2_b = df2_b.sort_values(by='Feature ID')
df3 = df2_b.drop_duplicates(subset='DTXSID', keep='first')

# Grab just score/authority columns
keep_strings = ['DTXCID_INDIVIDUAL_COMPONENT', '_score', '_authority']
# Filter columns based on the list of strings
df4 = df3.filter(regex='|'.join(keep_strings))

df5 = df4.groupby('DTXCID_INDIVIDUAL_COMPONENT').max()


score_cols = [col for col in df5.columns if col.endswith("_score")]
authority_cols = [col for col in df5.columns if col.endswith("_authority")]

# Calculate hazard averages
df5['Average Hazard Score'] = df5[score_cols].mean(axis=1)
df5['Average Quality Score'] = df5[authority_cols].mean(axis=1)
df5['Quality-Adjusted Hazard Score'] = df5['Average Hazard Score'] * df5['Average Quality Score']
df5['Number of end points with data available'] = df5[[col for col in df5.columns if col.endswith("_score")]].notna().sum(axis=1)
df5['Completeness Score'] = df5['Number of end points with data available'] / len(score_cols)

df5 = df5.reset_index()
debug_b = df5[['DTXCID_INDIVIDUAL_COMPONENT', 'Quality-Adjusted Hazard Score', 'Completeness Score']]

df5['Quality-Adjusted Hazard Score'] = df5['Quality-Adjusted Hazard Score'].fillna(0)
df5 = df5.sort_values(by='Quality-Adjusted Hazard Score', ascending=True)

if truncate > 0:
    df5 = df5.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
else:
    fig, ax = plt.subplots(figsize=(12, 40))

new_min, new_max = 0.01, 1
df5['Transparency'] = new_min + df5['Completeness Score'] * (new_max - new_min)

for index, row in df5.iterrows():
    ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color="blue", alpha=row['Transparency'])
    

plt.xlabel('Quality-Adjusted Hazard Score')
plt.ylabel('Candidate')
title_string = 'Hazard Results for SSM Chemical (approach 2): ' + chemical_name
plt.title(title_string)

plt.show()
'''


'''
##################################################################################################################
# Hazard results plotting - Bar chart transparency based on completeness - Sort based on Chemical metadata total score
##################################################################################################################

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

# Start off with just grabbing a single feature
df2 = df1[df1['Chemical Name'] == 'atenolol']
df2['Quality-Adjusted Hazard Score'] = df2['Quality-Adjusted Hazard Score'].fillna(0)


# Normalize all columns to graph
df2['Sources_norm'] = df2['Sources'] / df2['Sources'].max()
df2['Patents_norm'] = df2['Patents'] / df2['Patents'].max()
df2['Articles_norm'] = df2['Articles'] / df2['Articles'].max()
df2['PubMed Record Count_norm'] = df2['PubMed Record Count'] / df2['PubMed Record Count'].max()
df2['AMOS methods count_norm'] = df2['AMOS methods count'] / df2['AMOS methods count'].max()
df2['AMOS fact sheets count_norm'] = df2['AMOS fact sheets count'] / df2['AMOS fact sheets count'].max()
df2['AMOS spectra count_norm'] = df2['AMOS spectra count'] / df2['AMOS spectra count'].max()
df2['Presence in water lists count_norm'] = df2['Presence in water lists count'] / df2['Presence in water lists count'].max()

df2.fillna(0, inplace=True)
df2['Total_norm'] = df2['Sources_norm'] + df2['Patents_norm'] + df2['Articles_norm'] + df2['PubMed Record Count_norm'] + df2['AMOS methods count_norm'] + df2['AMOS fact sheets count_norm'] + df2['AMOS spectra count_norm'] + df2['Presence in water lists count_norm']
#df2 = df2.sort_values(by='Total_norm', ascending=True)

df2 = df2.sort_values(by='Quality-Adjusted Hazard Score', ascending=True)

truncate = 132
df2 = df2.tail(truncate)

new_min, new_max = 0.01, 1
#df2['Transparency'] = df2['Completeness Score'] + 0.3
df2['Transparency'] = new_min + df2['Completeness Score'] * (new_max - new_min)

fig, ax = plt.subplots(figsize=(12, 40))
#fig, ax = plt.subplots(figsize=(6, 10))

for index, row in df2.iterrows():
    ax.barh(row['PREFERRED_NAME'], row['Quality-Adjusted Hazard Score'], color="blue", alpha=row['Transparency'])
    

# Find index of SSM chemical
temp_boolean = df2[['SSM chemical']]

#temp_boolean = temp_boolean[::-1]
temp_boolean.reset_index(inplace=True, drop=True)
#SSM_chemical_index = df5[df5['SSM chemical']].index[0]
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]

# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)


# labels = [f"{cat} (Alpha: {alpha:.1f})" for cat, alpha in zip(df2['DTXSID'], df2['Transparency'])]
# ax.legend(bars, labels, title="Category & Transparency")
plt.xlabel('Quality-Adjusted Hazard Score')
plt.ylabel('Candidate')
title_string = 'Hazard Results for SSM Chemical: ' + df2['Chemical Name'].tolist()[0]
plt.title(title_string)

plt.show()
'''

'''

def create_gradient_legend(color, alpha_range=(0, 1), num_steps=10):
    """
    Creates a list of color patches with a gradient of alpha values for a custom legend. 
    
    Args:
        color (str or tuple): The base color to use for the gradient.
        alpha_range (tuple): Tuple defining the minimum and maximum alpha values (0-1).
        num_steps (int): Number of alpha steps to create in the gradient. 
    
    Returns:
        list: A list of Rectangle patches with varying alpha levels.
    """
    alpha_steps = [i / (num_steps - 1) for i in range(num_steps)]
    gradient_patches = [Rectangle((0, 0), 1, 1, color=color, alpha=alpha) for alpha in alpha_steps] 
    return gradient_patches

# Example usage
color = 'blue'  # Base color for the gradient
data = [10, 20, 30]

fig, ax = plt.subplots()
bars = ax.bar(range(len(data)), data, color=color)

# Generate gradient legend patches
gradient_legend = create_gradient_legend(color) 

# Add custom legend to the plot
ax.legend(gradient_legend, ['Low Alpha' ,'Medium Alpha', 'High Alpha', 'test'], title='Transparency Level')

plt.show()


import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Define the base color and alpha range
base_color = 'blue'  
alpha_range = [i/10 for i in range(11)]  

# Generate color gradient
gradient_colors = [(base_color, alpha) for alpha in alpha_range ]

# Create custom legend handle
custom_legend = Line2D([0], [0], color=base_color, lw=10, alpha=1, cmap='viridis', norm=plt.Normalize(0, 1), gradient=True)

# Sample data for bar plot
data = [5, 10, 15, 20]

fig, ax = plt.subplots()
bars = ax.bar(range(len(data)), data, color=base_color)

# Add legend with custom gradient entry
ax.legend([custom_legend], ['Intensity Gradient'], loc='upper right')

plt.show()

fig, ax = plt.subplots()
bars = ax.bar(['A', 'B', 'C'], [3, 5, 2], color='blue')

gradient = np.zeros((1, 256, 4))
gradient[:, :, 2] = 1
gradient[:, :, 3] = np.linspace(1, 0, 256)

class GradientLegendHandler(HandlerBase):
    def __init__(self, image_array, **kwargs):
        super().__init__(**kwargs)
        self.image = image_array
        
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        image = OffsetImage(self.image, zoom=0.5, transform=trans)
        ab = AnnotationBbox(image, (width / 2, height / 2), frameon=False, box_alignment=(0.5, 0.5), pad=0)
        return [ab]
    
legend_patch = mpatches.Patch(color='blue', label="Transparency (0 -> 1)")
legend = ax.legend(handles=[legend_patch], handler_map={legend_patch: GradientLegendHandler(gradient)}, loc='upper right')

plt.show()
'''
'''
# Set up the figure and axes
fig, ax = plt.subplots()

# Width of each bar
bar_width = 0.4

chemicals = df1['DTXSID'].tolist()

truncate = 100

#chemicals_reverse = chemicals[::-1]
chemicals_reverse = chemicals[:truncate]
chemicals_reverse = chemicals_reverse[::-1]

metadata_fields = ['Quality-Adjusted Hazard Score', 'Completeness Score']


# Try normalizing completeness score to highest QA-hazard score
max_QA = df1['Quality-Adjusted Hazard Score'].max()
max_completeness = df1['Completeness Score'].max()
#df1['Completeness Score'] = df1['Completeness Score'] * max_QA / max_completeness
df1['Completeness Score'] = df1['Completeness Score'] * max_QA


# # Create the bars
# for i, group in enumerate(metadata_fields):
#     ax.bar(x_pos + i * bar_width, values[:, i], width=bar_width, label=group)


# Create horizontal bars
for i, group in enumerate(metadata_fields):
    y_pos = np.arange(len(chemicals_reverse)) + i * bar_width
    temp_data = df1[group][:truncate]
    print(df1[group])
    print(temp_data)
    temp_data = temp_data[::-1]
    ax.barh(y_pos, temp_data, height=bar_width, label=group)

test = np.arange(len(chemicals_reverse))

# Set y-axis labels and title
ax.set_yticks(np.arange(len(chemicals_reverse)) + bar_width * (len(metadata_fields) - 1) / 2)
ax.set_yticklabels(chemicals_reverse)

ax.tick_params(axis='x', labelsize=3)

# Find index of SSM chemical
temp_boolean = df1[['SSM chemical']][:truncate]

temp_boolean = temp_boolean[::-1]
temp_boolean.reset_index(inplace=True, drop=True)
#SSM_chemical_index = df5[df5['SSM chemical']].index[0]
SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]

# Update ticklabel for the SSM chemical
yticklabels = ax.get_yticklabels()
yticklabels[SSM_chemical_index].set_color('green')
yticklabels[SSM_chemical_index].set_fontsize(16)

# Increase the figure size
#fig.set_size_inches(12, 40) 
fig.set_size_inches(6, 10) 
# Add a legend, show plot
plt.legend()
plt.show()
'''



'''
##################################################################################################################
# Hazard results plotting - Bar chart transparency based on completeness - Collapse on Structure/DTXCID
##################################################################################################################

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

#df1_DTXCID_counts = df1['DTXCID_INDIVIDUAL_COMPONENT'].value_counts()

# Start off with just grabbing a single feature
#df2 = df1[df1['Chemical Name'] == 'atenolol']
#df2 = df1[df1['Chemical Name'] == 'acetaminophen']
df2 = df1[df1['Chemical Name'] == 'theophylline']
#df2 = df1[df1['Chemical Name'] == 'caffeine']
#df2 = df1[df1['Chemical Name'] == 'atorvastatin']

#debug = df1[df1['Chemical Name'] == 'atenolol']
#debug = df1[df1['DTXCID_INDIVIDUAL_COMPONENT'] == 'DTXCID80197003'] 


df2['Quality-Adjusted Hazard Score'] = df2['Quality-Adjusted Hazard Score'].fillna(0)
df2 = df2.sort_values(by='Quality-Adjusted Hazard Score', ascending=True)

#df2_DTXCID_counts = df2['DTXCID_INDIVIDUAL_COMPONENT'].value_counts()


# Grab the highest structure count for atenolol
#df2 = df2[df2['DTXCID_INDIVIDUAL_COMPONENT'] == 'DTXCID7024530']
#df2 = df2[df2['DTXCID_INDIVIDUAL_COMPONENT'] == 'DTXCID202628']
#df2 = df2[df2['DTXCID_INDIVIDUAL_COMPONENT'] == 'DTXCID606']
#df2 = df2[df2['DTXCID_INDIVIDUAL_COMPONENT'] == 'DTXCID001336'] # Theophylline
#df2 = df2[df2['DTXCID_INDIVIDUAL_COMPONENT'] == 'DTXCID40232'] # caffeine
#df2 = df2[df2['DTXCID_INDIVIDUAL_COMPONENT'] == 'DTXCID80197003'] # atorvastatin

#truncate = 127
#df2 = df2.tail(truncate)

new_min, new_max = 0.01, 1
#df2['Transparency'] = df2['Completeness Score'] + 0.3
df2['Transparency'] = new_min + df2['Completeness Score'] * (new_max - new_min)

#fig, ax = plt.subplots(figsize=(12, 40))
fig, ax = plt.subplots(figsize=(6, 10))

for index, row in df2.iterrows():
    ax.barh(row['PREFERRED_NAME'], row['Quality-Adjusted Hazard Score'], color="blue", alpha=row['Transparency'])
    

# labels = [f"{cat} (Alpha: {alpha:.1f})" for cat, alpha in zip(df2['DTXSID'], df2['Transparency'])]
# ax.legend(bars, labels, title="Category & Transparency")
plt.xlabel('Quality-Adjusted Hazard Score')
plt.ylabel('Candidate')
title_string = 'Hazard Results for SSM Chemical: ' + df2['Chemical Name'].tolist()[0]
plt.title(title_string)

plt.show()

#df2_export = df2[df2['Feature ID'] == 793]
#df2_export.to_csv('WW2DW - Theophylline DTXCID001336 hazard results.csv', sep=',', encoding='utf-8', index=False)
'''


'''
##################################################################################################################
# Hazard results plotting - Multi-bar chart: Iterate through all SSM chemicals
##################################################################################################################
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

df2 = df1.drop_duplicates(subset=['Chemical Name', 'Ionization Mode'])
chemical_list = df2['Chemical Name'].tolist()
ionization_mode = df2['Ionization Mode'].tolist()

chemical_list_b = list(zip(chemical_list, ionization_mode))

bar_width = 0.4

for chemical in chemical_list_b:
    print('Working on chemical', chemical[0])
    df_temp = df1[(df1['Chemical Name'] == chemical[0]) * (df1['Ionization Mode'] == chemical[1])]
    df_temp = df_temp.sort_values(by='Quality-Adjusted Hazard Score', ascending=False)
    
    # Plot the non-truncated figures
    fig, ax = plt.subplots()
    # Width of each bar
    chemicals = df_temp['DTXSID'].tolist()
    truncate = len(chemicals)

    chemicals_reverse = chemicals[:truncate]
    chemicals_reverse = chemicals_reverse[::-1]
    metadata_fields = ['Quality-Adjusted Hazard Score', 'Completeness Score']

    # Try normalizing completeness score to highest QA-hazard score
    max_QA = df_temp['Quality-Adjusted Hazard Score'].max()
    max_completeness = df_temp['Completeness Score'].max()
    df_temp['Completeness Score'] = df_temp['Completeness Score'] * max_QA
    
    # Create horizontal bars
    for i, group in enumerate(metadata_fields):
        y_pos = np.arange(len(chemicals_reverse)) + i * bar_width
        temp_data = df_temp[group][:truncate]
        temp_data = temp_data[::-1]
        ax.barh(y_pos, temp_data, height=bar_width, label=group)
    
    test = np.arange(len(chemicals_reverse))
    
    # Set y-axis labels and title
    ax.set_yticks(np.arange(len(chemicals_reverse)) + bar_width * (len(metadata_fields) - 1) / 2)
    ax.set_yticklabels(chemicals_reverse)
    #ax.set_xlabel('Values')
    title_string = 'Hazard Results for SSM Chemical: ' + chemical[0]
    ax.set_title(title_string)
    
    ax.tick_params(axis='x', labelsize=3)
    
    # Find index of SSM chemical
    temp_boolean = df_temp[['SSM chemical']][:truncate]
    temp_boolean = temp_boolean[::-1]
    temp_boolean.reset_index(inplace=True, drop=True)
    
    yticklabels = ax.get_yticklabels()
    # Check if a SSM chemical was found and if so label on the graph
    if temp_boolean['SSM chemical'].sum() > 0:
    #if sum(temp_boolean) > 0:
        SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]
        # Update ticklabel for the SSM chemical
        yticklabels[SSM_chemical_index].set_color('green')
        yticklabels[SSM_chemical_index].set_fontsize(16)
    
    # Increase the figure size
    fig.set_size_inches(12, 40) 
    # Add a legend, show plot
    plt.legend()
    plt.show()    

    df_temp = df1[(df1['Chemical Name'] == chemical[0]) * (df1['Ionization Mode'] == chemical[1])]
    df_temp = df_temp.sort_values(by='Quality-Adjusted Hazard Score', ascending=False)

    # Plot the truncated figures
    # Set up the figure and axes
    fig, ax = plt.subplots()
    # Width of each bar
    chemicals = df_temp['DTXSID'].tolist()
    truncate = 25

    chemicals_reverse = chemicals[:truncate]
    chemicals_reverse = chemicals_reverse[::-1]
    metadata_fields = ['Quality-Adjusted Hazard Score', 'Completeness Score']

    # Try normalizing completeness score to highest QA-hazard score
    max_QA = df_temp['Quality-Adjusted Hazard Score'].max()
    max_completeness = df_temp['Completeness Score'].max()
    df_temp['Completeness Score'] = df_temp['Completeness Score'] * max_QA
    
    # Create horizontal bars
    for i, group in enumerate(metadata_fields):
        y_pos = np.arange(len(chemicals_reverse)) + i * bar_width
        temp_data = df_temp[group][:truncate]
        temp_data = temp_data[::-1]
        ax.barh(y_pos, temp_data, height=bar_width, label=group)
    
    test = np.arange(len(chemicals_reverse))
    
    # Set y-axis labels and title
    ax.set_yticks(np.arange(len(chemicals_reverse)) + bar_width * (len(metadata_fields) - 1) / 2)
    ax.set_yticklabels(chemicals_reverse)
    #ax.set_xlabel('Values')
    title_string = 'Hazard Results for SSM Chemical: ' + chemical[0]
    ax.set_title(title_string)
    
    ax.tick_params(axis='x', labelsize=3)
    
    # Find index of SSM chemical
    temp_boolean = df_temp[['SSM chemical']][:truncate]
    temp_boolean = temp_boolean[::-1]
    temp_boolean.reset_index(inplace=True, drop=True)
    
    yticklabels = ax.get_yticklabels()
    # Check if a SSM chemical was found and if so label on the graph
    if temp_boolean['SSM chemical'].sum() > 0:
    #if sum(temp_boolean) > 0:
        SSM_chemical_index = temp_boolean[temp_boolean['SSM chemical']].index[0]
        # Update ticklabel for the SSM chemical
        yticklabels[SSM_chemical_index].set_color('green')
        yticklabels[SSM_chemical_index].set_fontsize(16)
    
    # Increase the figure size
    #fig.set_size_inches(12, 40) 
    fig.set_size_inches(6, 10) 
    # Add a legend, show plot
    plt.legend()
    plt.show()
'''



##################################################################################################################
# Generate cleaned file - substance level for folks
##################################################################################################################
# Read in MS1 results for just the SSM chemicals
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df1 = pd.read_csv('WW2DW MS1 results for SSM chemicals (converted hazard values).csv')

# Remove feature duplicates if SSM chemical in more than one mode
df2 = df1.sort_values(by='Feature ID')
#df2 = df2.drop_duplicates(subset='DTXSID', keep='first')

col_list = df2.columns.tolist()

df2 = df2[['Feature ID', 'Chemical Name', 'Tracer DTXSID', 'Ionization Mode', 'DTXCID_INDIVIDUAL_COMPONENT', 'DTXSID', 'PREFERRED_NAME', 'SSM chemical', 
           'Sources', 'Patents', 'Articles', 'PubMed Record Count', 'AMOS methods count', 'AMOS fact sheets count', 'AMOS spectra count', 
           'Presence in water lists count', 'Quality-Adjusted Hazard Score', 'Completeness Score']]

df2 = df2.fillna(0)

# Create a column of the SSM chemical structure True/False so we can eventually collapse into structure
df2['SSM structure'] = df2.groupby(['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])['SSM chemical'].transform('any')

# Now collapse hazard/metadata on structure
# HAZARD: Create a column that has the maximum QAH score for a given structure within a feature
df2['Max_structure_QAH'] = df2.groupby(['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])['Quality-Adjusted Hazard Score'].transform('max')

df2_SSM_count = df2.drop_duplicates(subset=['Feature ID'])

# METADATA: sum across substances for structures
metadata_columns = ['Sources', 'Patents', 'Articles', 'PubMed Record Count', 'AMOS methods count', 'AMOS fact sheets count', 'AMOS spectra count', 'Presence in water lists count']

# Within each structure within each feature, sum the substances
df3 = df2.groupby(['Feature ID', 'Chemical Name', 'Ionization Mode', 'DTXCID_INDIVIDUAL_COMPONENT']).agg({col: 'sum' for col in metadata_columns})
df3 = df3.reset_index()

# Within each feature, normalize each metadata value to the maximum value
for col in metadata_columns:
    df3[f'{col}_norm'] = df3.groupby('Feature ID')[col].transform(lambda x: x / x.max())

df3 = df3.fillna(0)
df3['Total_norm'] = df3['Sources_norm'] + df3['Patents_norm'] + df3['Articles_norm'] + df3['PubMed Record Count_norm'] + df3['AMOS methods count_norm'] + df3['AMOS fact sheets count_norm'] + df3['AMOS spectra count_norm'] + df3['Presence in water lists count_norm']


# Get count of structures for each group
df3['structure_count'] = df3.groupby('Feature ID')['DTXCID_INDIVIDUAL_COMPONENT'].transform('count')

df3_SSM_count = df3.drop_duplicates(subset=['Feature ID'])

# Combine and grab only the needed columns
#df4 = pd.merge(df2[['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT', 'SSM structure', 'Max_structure_QAH']], df3[['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT', 'structure_count', 'Total_norm']], how='left', on=['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])
df4 = pd.merge(df3[['Feature ID', 'Chemical Name', 'Ionization Mode', 'DTXCID_INDIVIDUAL_COMPONENT', 'structure_count', 'Total_norm']], df2[['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT', 'SSM structure', 'Max_structure_QAH']], how='left', on=['Feature ID', 'DTXCID_INDIVIDUAL_COMPONENT'])
df4 = df4.drop_duplicates()

df4.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 

# Bring in MS2
# Read in the MS2 results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/neg')
df_ms2_neg = pd.read_csv('20210131_WW2DW_MS2_neg_CFMID_results_neg.csv')
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/pos')
df_ms2_pos = pd.read_csv('20250131_WW2DW_MS2_pos_CFMID_results_pos.csv')

# Need to add a RT tolerance for merging... or fix merge workflow ac 2/18/2025
df_ms2_neg.rename(columns={'SUM_SCORE':'MS2 neg'}, inplace=True) 
df_ms2_pos.rename(columns={'SUM_SCORE':'MS2 pos'}, inplace=True) 

df5 = pd.merge(df4, df_ms2_neg[['DTXCID', 'MS2 neg']], how='left', on='DTXCID')
df5 = pd.merge(df5, df_ms2_pos[['DTXCID', 'MS2 pos']], how='left', on='DTXCID')

#debug = df5[df5['DTXCID'] == 'DTXCID20208869']
#debug = df5[df5['Feature ID'] == 2901]

# De-duplicate merges where both modes data is present
df5 = df5.sort_values('MS2 neg', ascending=False).drop_duplicates(['Feature ID', 'DTXCID', 'MS2 neg']).sort_index()
df5 = df5.sort_values('MS2 pos', ascending=False).drop_duplicates(['Feature ID', 'DTXCID', 'MS2 pos']).sort_index()

rank_columns = ['Total_norm', 'Max_structure_QAH', 'MS2 pos', 'MS2 neg']

# Generate rank columns
for col in rank_columns:
    df5[f'{col}_rank'] = df5.groupby('Feature ID')[col].rank(method='dense', ascending=False)



df5_SSM_count = df5.drop_duplicates(subset=['Feature ID'])


df3_subset = df3.drop(columns=['Total_norm', 'structure_count'])
df3_subset.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 
df6 = pd.merge(df5, df3_subset, how='left', on=['Feature ID', 'Chemical Name', 'Ionization Mode', 'DTXCID'])

df6 = df6.rename(columns={'structure_count': 'Number of candidate structures for feature', 'Total_norm': 'Chemical metadata score',
                          'SSM structure': 'Correct SSM structure', 'Max_structure_QAH': 'Hazard score', 
                          'MS2 neg': 'MS2 score (neg)', 'MS2 pos': 'MS2 score (pos)', 'Total_norm_rank': 'Chemical metadata rank',
                          'Max_structure_QAH_rank': 'Hazard rank', 'MS2 neg_rank': 'MS2 rank (neg)', 'MS2 pos_rank': 'MS2 rank (pos)'
                          })

df6 = df6.apply(lambda x: pd.to_numeric(x, errors='ignore') if x.dtype == 'object' else x)

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis')
df6.to_csv('20250220 WW2DW tripod leg results.csv', sep=',', float_format="%.15g", quoting=csv.QUOTE_NONE, encoding='utf-8', index=False)
df6.to_excel('20250220 WW2DW tripod leg results.xlsx', index=False)

print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))  