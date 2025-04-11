# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 07:18:54 2025

@author: AChao
"""

import time
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap

time_list = []
start_time = time.time()
time_list.append(time.time())


# Export dataframe with all the raw and collapsed values to CSV
#os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
#df1 = pd.read_csv('WW2DW Data Analysis file 4 (Collapsed to structure-level)(Metadata summed and normalized, MS2 normalized, Hazard max value).csv')

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v3')
#df1 = pd.read_csv('WW2DW Data Analysis file 4 (Collapsed to structure-level)(Metadata summed and normalized, MS2 normalized, Hazard max value) v3.csv')
df1 = pd.read_csv('WW2DW Data Analysis file 4 (Structure level only, all results) v3 (20250227).csv')

df_substance = pd.read_csv('WW2DW Data Analysis file 3 (Substance and structure level, all results) v3 (20250226).csv')
df_substance_features = df_substance['Feature ID'].tolist()
df_substance_features = list(set(df_substance_features))
debug = df_substance[df_substance['Feature ID'] == 200]
df_substance_SSM = df_substance[df_substance['SSM chemical feature'] == 'Y']

df_structure_SSM = df1[df1['SSM chemical feature'] == 'Y']


df1 = pd.read_csv('WW2DW Data Analysis file 4 (Structure level only, all results) v3 (20250226).csv')
df1_cols = df1.columns.tolist()

df2 = df1.rename(columns={'Structure_Quality-Adjusted Hazard Score': 'Hazard score', 'MS2 quotient score': 'MS2 score', 'Structure_total_norm': 'Metadata score', 
                          'Final Occurrence Count':'Occurrence Count', 'Median blanksub mean feature abundance':'Median abundance'})
df2_cols = df2.columns.tolist()

df2['log Median abundance'] = df2['Median abundance'].apply(lambda x: np.log10(x) if x > 0 else np.nan)

debug = df2[df2['Feature ID'] == 200]

df4 = df2.drop_duplicates(subset=['Feature ID']) 

# Grab the level 1 identifications
df_IDs = pd.read_csv('WW2DW Level 1s with DSSTox identifiers.csv')
df_IDs['Level 1 ID'] = True



# Select SSM chemical features only
df_SSM = df1[df1['SSM chemical feature'] == 'Y']

# Select non-SSM chemical features and only include those that have MS2 data
# df2 = df1.groupby('Feature ID').filter(lambda x: (x['MS2 score'] > 0).any())
# df2 = df2[df2['SSM chemical feature'] != 'Y']

#df3 = df2[df2['Feature ID'] == 247]

df3 = df2.groupby('Feature ID').filter(lambda x: (x['MS2 score'] > 0).any())
#df3 = df3[df3['Structure_Sources'] > 5]

#####df3 = df3[df3['Metadata score'] > 0.5]
#####df3 = df3[df3['Structure_Presence in water lists count'] > 0]


df3 = pd.merge(df3, df_IDs[['DTXCID', 'Level 1 ID']], how='left', on='DTXCID')

# Check how many features match level 1 ID's in overall dataset
test_ID = pd.merge(df2, df_IDs[['DTXCID', 'Level 1 ID']], how='left', on='DTXCID')
test_ID_b = test_ID[test_ID.groupby('Feature ID')['Level 1 ID'].transform('any')]
test_SSM = test_ID_b[test_ID_b['SSM chemical feature'] == 'Y']
test_non_SSM = test_ID_b[test_ID_b['SSM chemical feature'] != 'Y']
test_non_SSM_features = test_non_SSM['Feature ID'].tolist()
test_non_SSM_features = list(set(test_non_SSM_features))
test_SSM_features = test_SSM['Feature ID'].tolist()
test_SSM_features = list(set(test_SSM_features))

df3_SSM = df3[df3['SSM chemical feature'] == 'Y']

df3_ID = df3[df3.groupby('Feature ID')['Level 1 ID'].transform('any')]
df3_ID_b = df3_ID[df3_ID['SSM chemical feature'] != 'Y']

SSM_features = df3_SSM['Feature ID'].tolist()
SSM_features = list(set(SSM_features))

df3['color'] = np.where(df3['SSM structure'], 'red', 'blue')
df4['color'] = np.where(df4['SSM structure'], 'red', 'blue')

# colors = ['#00008B', "#4169E1", "#FF4500", "#FF0000", "#FF6347"]
# colors = ['#00008B', "#4169E1", "#FFFF00", "#850101", "#FF0000"]
# colors = ['#00008B', "#4169E1", "#FFFF00", "#FF8C00", "#FF0000"]
# colors = ['#00008B', "#4169E1", "#DDAA00", "#FFEE00", "#FF0000"]
# cmp = ListedColormap(colors)

#cmp='viridis'
# cmp='inferno'

cmp='YlOrRd'
font_size = 20

# Metadata vs MS2, size = occurrence, color = hazard
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.8, s=df3['Occurrence Count'] * 50, c=df3['Hazard score'], cmap=cmp)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)
plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
#plt.title('Size = Occurrence Count, Color = Hazard', fontsize=font_size)


# Print all SSM features via loop
for feature in SSM_features:
    if feature == 2136:
        temp_df = df3[df3['Feature ID'] == feature]
        #temp_df = df2[df2['Feature ID'] == feature]
        
        dummy_x = np.nan
        dummy_y = np.nan
        dummy_color = df3['Hazard score'].max()
        dummy_df = pd.DataFrame({'Metadata score':[8], 'MS2 score':[1], 'Occurrence Count':0.0000001, 'Hazard score':[12]})
        
        temp_df = pd.concat([temp_df, dummy_df], ignore_index=True)
        
        plt.figure(figsize=(10, 8))
        plt.xticks(fontsize=font_size-4)
        plt.yticks(fontsize=font_size-4)
        scatter = plt.scatter(temp_df['Metadata score'], temp_df['MS2 score'], alpha = 0.8, s=temp_df['Occurrence Count'] * 50, c=temp_df['Hazard score'], cmap=cmp)
        cbar = plt.colorbar(scatter)
        cbar.set_label('Hazard Score', size=font_size) 
        cbar.ax.tick_params(labelsize=font_size-4)
        plt.xlim(-0.5, 8.5)
        plt.ylim(-0.05, 1.05)
        plt.xlabel("Metadata score", fontsize=font_size)
        plt.ylabel("MS2 score", fontsize=font_size)
        #plt.title('Feature '+str(feature)+': Size = Occurrence Count', fontsize=font_size)
    


# Plot A - True chemicals, true structures
df5 = df3[df3['SSM structure'] == True]

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df3['Hazard score'].max()
dummy_df = pd.DataFrame({'Metadata score':[8], 'MS2 score':[1], 'Occurrence Count':0.0000001, 'Hazard score':[12]})
df5 = pd.concat([df5, dummy_df], ignore_index=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df5['Metadata score'], df5['MS2 score'], alpha = 0.8, s=df5['Occurrence Count'] * 50, c=df5['Hazard score'], cmap=cmp)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)
plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('Plot A: SSM features, true candidates', fontsize=font_size)

# Plot 8 - True chemicals, other structures
df5 = df3[df3['SSM chemical feature'] == 'Y']
df5 = df5[df5['SSM structure'] == False]

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df3['Hazard score'].max()
dummy_df = pd.DataFrame({'Metadata score':[8], 'MS2 score':[1], 'Occurrence Count':0.0000001, 'Hazard score':[12]})
df5 = pd.concat([df5, dummy_df], ignore_index=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df5['Metadata score'], df5['MS2 score'], alpha = 0.8, s=df5['Occurrence Count'] * 50, c=df5['Hazard score'], cmap=cmp)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)
plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('Plot B: SSM features, incorrect candidates', fontsize=font_size)

# Plot C - Other chemical features
df5 = df3[df3['SSM chemical feature'] != 'Y']

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df3['Hazard score'].max()
dummy_df = pd.DataFrame({'Metadata score':[8], 'MS2 score':[1], 'Occurrence Count':0.0000001, 'Hazard score':[12]})
df5 = pd.concat([df5, dummy_df], ignore_index=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df5['Metadata score'], df5['MS2 score'], alpha = 0.8, s=df5['Occurrence Count'] * 50, c=df5['Hazard score'], cmap=cmp)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)
plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('Plot C: Non-SSM features', fontsize=font_size)


# Plot D - Level 1's, true structures
df5 = df3_ID_b[df3_ID_b['Level 1 ID'] == True]

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df3['Hazard score'].max()
dummy_df = pd.DataFrame({'Metadata score':[8], 'MS2 score':[1], 'Occurrence Count':0.0000001, 'Hazard score':[12]})
dummy_df_b = pd.DataFrame({'Metadata score':[0], 'MS2 score':[0], 'Occurrence Count':0.0000001, 'Hazard score':[0]})
df5 = pd.concat([df5, dummy_df], ignore_index=True)
df5 = pd.concat([df5, dummy_df_b], ignore_index=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df5['Metadata score'], df5['MS2 score'], alpha = 0.8, s=df5['Occurrence Count'] * 50, c=df5['Hazard score'], cmap=cmp)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)
plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('Plot D: Level 1 IDs, true candidates', fontsize=font_size)

# Plot E - Level 1's, other structures
df5 = df3_ID_b[df3_ID_b['Level 1 ID'] != True]

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df3['Hazard score'].max()
dummy_df = pd.DataFrame({'Metadata score':[8], 'MS2 score':[1], 'Occurrence Count':0.0000001, 'Hazard score':[12]})
dummy_df_b = pd.DataFrame({'Metadata score':[0], 'MS2 score':[0], 'Occurrence Count':0.0000001, 'Hazard score':[0]})
df5 = pd.concat([df5, dummy_df], ignore_index=True)
df5 = pd.concat([df5, dummy_df_b], ignore_index=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df5['Metadata score'], df5['MS2 score'], alpha = 0.8, s=df5['Occurrence Count'] * 50, c=df5['Hazard score'], cmap=cmp)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)
plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('Plot E: Level 1 IDs, incorrect candidates', fontsize=font_size)



# 2D plot - SSM true chemical vs others
plt.figure(figsize=(8, 8))
# Change tick size
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.3, s=df3['Hazard score'] ** 3 * 5 + 10, c=df3['color'])
plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('Size = Hazard, Red = True SSM chemical', fontsize=font_size)



# # plt.figure(figsize=(10, 8))
# # scatter = plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.5, s=df3['log Median abundance'] ** 3 * 10 + 10, c=df3['Hazard score'], cmap=cmp)
# # plt.colorbar(scatter, label='Hazard Score')
# # plt.xlabel("Metadata score")
# # plt.ylabel("MS2 score")
# # plt.title('Size = log Median abundance, Color = Hazard')

# plt.figure(figsize=(10, 8))
# scatter = plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.5, s=df3['Hazard score'] ** 3 * 5 + 10, c=df3['Occurrence Count'], cmap=cmp)
# plt.colorbar(scatter, label='Occurrence Count')
# plt.xlabel("Metadata score")
# plt.ylabel("MS2 score")
# plt.title('Size = Hazard, Color = Occurrence Count')

# plt.figure(figsize=(10, 8))
# scatter = plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.5, s=df3['Hazard score'] ** 3 * 5 + 10, c=df3['Median abundance'], cmap=cmp)
# plt.colorbar(scatter, label='Median abundance')
# plt.xlabel("Metadata score")
# plt.ylabel("MS2 score")
# plt.title('Size = Hazard, Color = Median abundance')

# plt.figure(figsize=(10, 8))
# scatter = plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.5, s=df3['Hazard score'] ** 3 * 5 + 10, c=df3['log Median abundance'], cmap=cmp)
# plt.colorbar(scatter, label='log Median abundance')
# plt.xlabel("Metadata score")
# plt.ylabel("MS2 score")
# plt.title('Size = Hazard, Color = log Median abundance')

# plt.figure(figsize=(10, 8))
# scatter = plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.5, s=df3['Occurrence Count'] ** 3 * 5 + 10, c=df3['Median abundance'], cmap=cmp)
# plt.colorbar(scatter, label='Median abundance')
# plt.xlabel("Metadata score")
# plt.ylabel("MS2 score")
# plt.title('Size = Occurrence Count, Color = Median abundance')

# plt.figure(figsize=(10, 8))
# scatter = plt.scatter(df3['Metadata score'], df3['MS2 score'], alpha = 0.5, s=df3['Occurrence Count'] * 200 + 10, c=df3['Median abundance'], cmap=cmp)
# plt.colorbar(scatter, label='Median abundance')
# plt.xlabel("Metadata score")
# plt.ylabel("MS2 score")
# plt.title('Size = Occurrence Count, Color = Median abundance')



# plt.figure(figsize=(8, 8))
# scatter = plt.scatter(df3['Occurrence Count'], df3['log Median abundance'], alpha = 0.3, s=df3['Hazard score'] ** 3 * 5 + 10, c=df3['color'])
# #plt.colorbar(scatter, label='Median abundance')
# plt.xlabel("Occurrence Count")
# plt.ylabel("log Median abundance")
# plt.title('Size = Hazard Score')




# plt.figure(figsize=(8, 8))
# scatter = plt.scatter(df4['Occurrence Count'], df4['log Median abundance'], alpha = 0.2, s=200, c=df4['color'])
# plt.xlabel("Occurrence Count")
# plt.ylabel("log Median abundance")



# 3D plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Create the scatter plot
scatter = ax.scatter(df3['Metadata score'], df3['MS2 score'], df3['Hazard score'], cmap='viridis', s=15, alpha = 0.4, c=df3['color'])

# Add labels and title
ax.set_xlabel('Metadata score')
ax.set_ylabel('MS2 score')
ax.set_zlabel('Hazard score')
#ax.set_title('3D Scatter Plot')




##################################################################################################################
# Tripod leg 1: Metadata bar chart
##################################################################################################################
truncate = 25
select_SSM_structure_only = False
use_transparency = False

df6 = df1[df1['Feature ID'] == 2136]
df6_cols = df6.columns.tolist()

# Identify the DTXCID associated with the SSM chemical
df6_true = df6[df6['SSM structure'] == True]
SSM_DTXCID = df6_true['DTXCID'].tolist()[0]


# Do some filtering on candidates with only 1 source
df6 = df6[df6['Structure_total_norm'] != 0]
df6 = df6[df6['Structure_Sources'] > 1]


df6 = df6.sort_values(by='Structure_total_norm', ascending=True)

df6_plot = df6[['DTXCID', 'Structure_Sources_norm', 'Structure_Patents_norm', 'Structure_Articles_norm', 'Structure_PubMed Record Count_norm', 
                'Structure_AMOS methods count_norm', 'Structure_AMOS fact sheets count_norm', 
                'Structure_AMOS spectra count_norm', 'Structure_Presence in water lists count_norm']]

df6_plot = df6_plot.rename(columns={'Structure_Sources_norm': 'Sources', 'Structure_Patents_norm': 'Patents', 'Structure_Articles_norm': 'Articles', 
                                    'Structure_PubMed Record Count_norm': 'PubMed Records', 'Structure_AMOS methods count_norm': 'Methods',
                                    'Structure_AMOS fact sheets count_norm': 'Fact sheets', 'Structure_AMOS spectra count_norm': 'Spectra',
                                    'Structure_Presence in water lists count_norm':'Water Lists'})

#df5_plot['PREFERRED_NAME'] = df5_plot['PREFERRED_NAME'].str[0:30]

if truncate > 0:
    df6_plot = df6_plot.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
    SSM_font_size = 12
    label_sizes = 12
    ax.tick_params(axis='y', labelsize=label_sizes)
    ax.tick_params(axis='x', labelsize=label_sizes + 4)
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    SSM_font_size = 22
    label_sizes = 22
    ax.tick_params(axis='y', labelsize=label_sizes - 4)
    ax.tick_params(axis='x', labelsize=label_sizes + 12)

# Initialize bottom for stacking
bottom = [0] * len(df6_plot)

# Loop through each row (each group of stacked bars)
for i, label in enumerate(df6_plot['DTXCID']):
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


#plt.legend(fontsize=18)
# plt.ylabel('Candidate', fontsize=label_sizes+6)
# plt.xlabel('Total Metadata Score', fontsize=label_sizes+6)

# Set size and color of DTXCID corresponding to SSM chemical
for label in ax.get_yticklabels():
    if label.get_text() == SSM_DTXCID:
        label.set_fontsize(SSM_font_size + 6)  # Set font size to 16
        label.set_color('red')

legend = ax.legend(fontsize=label_sizes, loc='lower right')

# Iterate through the legend handles (patches)
for handle in legend.legendHandles:
    # Set the alpha value for each handle
    handle.set_alpha(1)  # Adjust the value between 0 (fully transparent) and 1 (fully opaque)


# plt.ylabel('Candidate', fontsize=label_sizes+12)
# plt.xlabel('Metadata Score', fontsize=label_sizes+12)

# title_string = 'Metadata Results for Feature 2136'
#plt.title(title_string, fontsize=label_sizes + 12)
plt.show()

##################################################################################################################
# Tripod leg 2: MS2 bar chart
##################################################################################################################
truncate = 25
df8= df1[df1['Feature ID'] == 2136]

# Identify the DTXCID associated with the SSM chemical
df8_true = df8[df8['SSM structure'] == True]
SSM_DTXCID = df8_true['DTXCID'].tolist()[0]

df8['MS2 quotient score'] = df8['MS2 quotient score'].fillna(0)
df8 = df8.sort_values(by='MS2 quotient score', ascending=True)

if truncate > 0:
    df8 = df8.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
    SSM_font_size = 12
    label_sizes = 12
    ax.tick_params(axis='y', labelsize=label_sizes)
    ax.tick_params(axis='x', labelsize=label_sizes + 4)
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    SSM_font_size = 22
    label_sizes = 22
    ax.tick_params(axis='y', labelsize=label_sizes - 4)
    ax.tick_params(axis='x', labelsize=label_sizes + 12)

bars = ax.barh(df8['DTXCID'], df8['MS2 quotient score'], alpha=0.6)

# plt.ylabel('Candidate', fontsize=label_sizes+12)
# plt.xlabel('MS2 Score', fontsize=label_sizes+12)

# Set size and color of DTXCID corresponding to SSM chemical
for bar, label in zip(bars, ax.get_yticklabels()):
    if label.get_text() == SSM_DTXCID:
        label.set_fontsize(SSM_font_size + 6)  # Set font size to 16
        label.set_color('red')
        bar.set_alpha(1)

# title_string = 'MS2 Results for Feature 2136'
# plt.title(title_string, fontsize=label_sizes + 12)
plt.show()




##################################################################################################################
# Tripod leg 3: Hazard bar chart
##################################################################################################################
# Start off with just grabbing a single feature
truncate = 25
select_SSM_structure_only = False
use_transparency = False

df7 = df1[df1['Feature ID'] == 2136]

# Identify the DTXCID associated with the SSM chemical
df7_true = df7[df7['SSM structure'] == True]
SSM_DTXCID = df7_true['DTXCID'].tolist()[0]

#chemical_name = df7['Chemical Name'].tolist()[0]

# # Remove feature duplicates if SSM chemical in more than one mode
# df7 = df7.sort_values(by='Feature ID')
# df7 = df7.drop_duplicates(subset='DTXSID', keep='first')

# Need to grab completeness scores from other file
#df_substances = pd.read_csv('WW2DW Data Analysis file 3 (Substance and structure level, all results) v3 (20250226).csv')
df_substances_col = df_substances.columns.tolist()
df7 = pd.merge(df7, df_substances[['Feature ID', 'DTXCID', 'Hazard Completeness Score']], how='left', on=['Feature ID', 'DTXCID'])
df7 = df7.drop_duplicates(subset=['Feature ID', 'DTXCID'], keep='first')

df7['Structure_Quality-Adjusted Hazard Score'] = df7['Structure_Quality-Adjusted Hazard Score'].fillna(0)
df7 = df7.sort_values(by=['Structure_Quality-Adjusted Hazard Score', 'Hazard Completeness Score'], ascending=True)

# Remove structures with no score
#df8 = df7[df7['Structure_Quality-Adjusted Hazard Score'] > 0]
df8 = df7.copy()
debug = df8[['DTXCID', 'Structure_Quality-Adjusted Hazard Score', 'Hazard Completeness Score']]

if truncate > 0:
    df8 = df8.tail(truncate)
    fig, ax = plt.subplots(figsize=(4, 10))
    SSM_font_size = 12
    label_sizes = 12
    ax.tick_params(axis='y', labelsize=label_sizes)
    ax.tick_params(axis='x', labelsize=label_sizes + 4)
else:
    fig, ax = plt.subplots(figsize=(12, 40))
    SSM_font_size = 22
    label_sizes = 22
    ax.tick_params(axis='y', labelsize=label_sizes - 4)
    ax.tick_params(axis='x', labelsize=label_sizes + 12)
    
new_min, new_max = 0.01, 1
df8['Transparency'] = new_min + df8['Hazard Completeness Score'] * (new_max - new_min)

use_select_color = False


# for index, row in df8.iterrows():
#     #bar_color = "magenta" if row['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID else "blue"
#     if use_select_color == True:
#         bar_color = "blue" if row['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID else "red"
#         ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color=bar_color, alpha=row['Transparency'])
#         bar_outline = 10 if row['DTXCID_INDIVIDUAL_COMPONENT'] == SSM_DTXCID else 0
#     else:
#         bar_color = "red"
#         ax.barh(row['DTXCID'], row['Structure_Quality-Adjusted Hazard Score'], color=bar_color, alpha=row['Transparency'])
#         bar_outline = 0
#     #ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color='blue', edgecolor='red', linewidth=bar_outline, alpha=row['Transparency'])  
#     #ax.barh(row['DTXCID_INDIVIDUAL_COMPONENT'], row['Quality-Adjusted Hazard Score'], color='red', edgecolor='blue', linewidth=bar_outline, alpha=row['Transparency'])  


for index, row in df8.iterrows():
    ax.barh(row['DTXCID'], row['Structure_Quality-Adjusted Hazard Score'], color='red', alpha=row['Transparency'])

# Set size and color of DTXCID corresponding to SSM chemical
for label in ax.get_yticklabels():
    if label.get_text() == SSM_DTXCID:
        label.set_fontsize(SSM_font_size + 6)  # Set font size to 16
        label.set_color('red')


# plt.xlabel('Hazard Score', fontsize=label_sizes+12)
# plt.ylabel('Candidate', fontsize=label_sizes+12)
# title_string = 'Hazard Results for Feature 2136'
# plt.title(title_string, fontsize=label_sizes+12)

plt.show()

##################################################################################################################
# Create red gradient rectangle for hazard plot legend
##################################################################################################################
fig, ax = plt.subplots(figsize=(6,2))
gradient = np.linspace(0.1, 1, 1000).reshape(1, -1)

red_color = np.zeros((1, 1000, 4))
red_color[:, :, 0] = 1  # Set the red channel to 1
red_color[:, :, 3] = gradient  # Apply the gradient to the alpha channel

ax.imshow(red_color, extent=[0, 5, 0, 2], aspect='auto')

ax.set_xticks([])
ax.set_yticks([])
ax.set_frame_on(False)

plt.show()


# Plot MS2 availability vs. abundance
df_ms2 = df1[['Feature ID', 'DTXCID', 'MS2 quotient score', 'Median blanksub mean feature abundance', 'Final Occurrence Count', 'SSM chemical feature']]
df_ms2 = df_ms2.drop_duplicates()

df_ms2 = df_ms2[df_ms2['SSM chemical feature'] != 'Y']

df_ms2['Feature has MS2'] = df_ms2.groupby('Feature ID')['MS2 quotient score'].transform(lambda x: x.notna().any())

df_ms2['log Median abundance'] = df_ms2['Median blanksub mean feature abundance'].apply(lambda x: np.log10(x) if x > 0 else np.nan)

debug = df_ms2[df_ms2['Feature has MS2'] == True]
debug = df1[df1['Feature ID'] == 5838]
debug = df_ms2[df_ms2['Feature has MS2'] == False]

df_ms2_b = df_ms2.drop_duplicates(subset=['Feature ID'], keep='first')

df_ms2_b['color'] = np.where(df_ms2_b['Feature has MS2'], 'red', 'blue')

df_ms2_true = df_ms2_b[df_ms2_b['Feature has MS2'] == True]
df_ms2_false = df_ms2_b[df_ms2_b['Feature has MS2'] == False]

# 2D plot - Has MS2 vs. not
plt.figure(figsize=(8, 8))
# Change tick size
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
#plt.scatter(df_ms2_b['Median blanksub mean feature abundance'], df_ms2_b['Final Occurrence Count'], alpha = 0.3, s=50, c=df_ms2_b['color'])
plt.scatter(df_ms2_b['log Median abundance'], df_ms2_b['Final Occurrence Count'], alpha = 0.3, s=50, c=df_ms2_b['color'])
plt.xlabel("Median feature abundance", fontsize=font_size)
plt.ylabel("Feature occurrence count", fontsize=font_size)
plt.title('Size = Hazard, Red = Has MS2l', fontsize=font_size)

# violin plots
df_ms2_b['Feature has MS2'] = df_ms2_b['Feature has MS2'].astype(bool)
plt.figure(figsize=(10,6))
sns.violinplot(x=df_ms2_b['Feature has MS2'], y=df_ms2_b['Median blanksub mean feature abundance'])
plt.xlabel("Feature has MS2 acquired", fontsize=font_size)
plt.ylabel("Feature median abundance", fontsize=font_size)
#plt.title('Size = Hazard, Red = Has MS2l', fontsize=font_size)

df_ms2_b['Feature has MS2'] = df_ms2_b['Feature has MS2'].astype(bool)
plt.figure(figsize=(10,6))
sns.violinplot(x=df_ms2_b['Feature has MS2'], y=df_ms2_b['log Median abundance'])
plt.xlabel("Feature has MS2 acquired", fontsize=12)
plt.ylabel("Feature median abundance (log-scale)", fontsize=12)
#plt.title('Size = Hazard, Red = Has MS2l', fontsize=font_size)




print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))  