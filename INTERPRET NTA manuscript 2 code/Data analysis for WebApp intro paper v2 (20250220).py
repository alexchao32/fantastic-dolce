# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:48:19 2025

@author: AChao

2/20/2025
Create output files from newly processed Ww2DW that has the newest database snapshot

Don't down-select SSM chemical features, should end up with:

-	Metadata scores: Collapsed on structure by summing across structure substances, normalized across all candidates for a given feature
-	MS2 scores: Normalized across all candidates for a given feature
-	Hazard scores: Collapsed on structure by taking the highest hazard score from the structure substances, NOT normalized but kept as raw values
-	Detection frequencies of features
"""

import time
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle 
import csv
import glob

time_list = []
start_time = time.time()
time_list.append(time.time())

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (20250220)')
'''
# Grab SSM chemicals from tracer summary sheet
df1 = pd.read_excel('20250220_WW2DW_MS1_analysis_NTA_WebApp_results.xlsx', sheet_name='Tracer Summary')
# Rename DTXSID column of tracer summary to not collide with DTXSID's of DSSTox/CHEM results
df1.rename(columns={'DTXSID':'Tracer DTXSID'}, inplace=True) 

# Grab the WebApp MS1 chemical results
df2 = pd.read_excel('20250220_WW2DW_MS1_analysis_NTA_WebApp_results.xlsx', sheet_name='Chemical Results')
'''
df2_cols = df2.columns.tolist()
# List of columns to drop
columns_to_drop = ['INPUT', 'MONOISOTOPIC_MASS_INDIVIDUAL_COMPONENT', 'SMILES_INDIVIDUAL_COMPONENT', 'CASRN', 'INCHIKEY', 'IUPAC_NAME', 'MOLECULAR_FORMULA', 'MONOISOTOPIC_MASS', 
                   'EXPOCAST_MEDIAN_EXPOSURE_PREDICTION_MG/KG-BW/DAY', 'EXPOCAST', 'NHANES', 'TOXCAST_PERCENT_ACTIVE', 'TOXCAST_NUMBER_OF_ASSAYS/TOTAL', 'MASS_DIFFERENCE', 
                   'FOUND_BY']

# Add on column indicating which features correspond to SSM chemicals
df3 = df2.drop(columns=columns_to_drop)
df3_cols = df3.columns.tolist()
df3_SSM = pd.merge(df1, df2, how='left', on='Feature ID')

df3_SSM['SSM chemical feature'] = 'Y'
df3_SSM = df3_SSM[['Feature ID', 'SSM chemical feature']]
df3_SSM = df3_SSM.drop_duplicates()

df4 = pd.merge(df3, df3_SSM, how='left', on='Feature ID')


# Add on column indicating which substance matches the true SSM chemical substance
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
df_SSM_names = pd.read_csv('SSM chemical DTXSID and names.csv')
df_SSM_names['SSM chemical substance'] = 'Y'
df4_SSM = df4[df4['SSM chemical feature'] == 'Y']
df4_SSM = df4_SSM.drop_duplicates()
df5 = pd.merge(df4_SSM, df_SSM_names[['PREFERRED_NAME', 'SSM chemical substance']], how='left', on='PREFERRED_NAME') 

df5 = df5.drop_duplicates()

#debug = df5[['Feature ID', 'PREFERRED_NAME', 'SSM chemical feature', 'SSM chemical substance']]
df6 = pd.merge(df4, df5[['Feature ID', 'PREFERRED_NAME', 'SSM chemical substance']], how='left', on=['Feature ID', 'PREFERRED_NAME'])

#debug = df6[df6['Feature ID'] == 200]

# This first file has the chemical substances for all features. It also has the raw hazard data, and columns indicating which feature is a SSM feature, and which substance is the SSM chemical
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
#df6.to_csv('WW2DW Data Analysis file 1 (substances, raw hazard values).csv', sep=',', encoding='utf-8', index=False)


# Grab all the Comptox water lists searches together
def combine_csvs():
    """Reads all CSV files in the current directory and concatenates them into a single DataFrame."""
    all_files = glob.glob("*.csv")
    
    if not all_files:
        raise FileNotFoundError("No CSV files found in the current directory.")
    
    all_df = []
    for f in all_files:
        df = pd.read_csv(f)
        all_df.append(df)
    
    merged_df = pd.concat(all_df, ignore_index=True)
    return merged_df


os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2/Comptox search for water lists')
df_water_lists = combine_csvs()
df_water_lists_cols = df_water_lists.columns.tolist()

# Now calculate/count presence in water lists
df_water_lists['Presence in water lists count'] = df_water_lists.iloc[:, 1:].apply(lambda row: row.eq('Y').sum(), axis=1)
#debug = df_water_lists.head(20)


df6 = df6.drop('DTXSID', axis=1)
debug = df6.head(20)
debug = df6[df6['Feature ID'] == 200]

# Merge water lists count and DTXSID's onto chemicals & raw hazard scores dataframe
df7 = pd.merge(df6, df_water_lists[['PREFERRED_NAME', 'DTXSID', 'Presence in water lists count']], how='left', on='PREFERRED_NAME')
debug = df7.head(20)

dtxsid_only = df7[['DTXSID']]
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
#dtxsid_only.to_csv('WW2DW Data Analysis - DTXSID list all features.csv', sep=',', encoding='utf-8', index=False)


# Read in AMOS files and concatenate
def combine_excel_sheets(directory, sheet_name):
    """
    Reads a specific sheet from all XLSX files in a directory and combines them into one DataFrame.

    Args:
        directory (str): The path to the directory containing the XLSX files.
        sheet_name (str): The name of the sheet to read from each file.

    Returns:
        pd.DataFrame: A DataFrame containing the combined data from all files, 
                      or an empty DataFrame if no files are found.
    """
    all_files = glob.glob(os.path.join(directory, "*.xlsx"))
    
    if not all_files:
        print(f"No XLSX files found in directory: {directory}")
        return pd.DataFrame()

    all_df = []
    for file in all_files:
        print('Working on file: ', file, sheet_name)
        try:
            df = pd.read_excel(file, sheet_name=sheet_name)
            all_df.append(df)
        except Exception as e:
            print(f"Error reading sheet '{sheet_name}' from file '{file}': {e}")
    
    if not all_df:
         return pd.DataFrame()
    
    combined_df = pd.concat(all_df, ignore_index=True)
    return combined_df

# Example usage:
directory_path = 'L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2/AMOS search for sources'
target_sheet = "Substances"
target_sheet_b = "Records"

df_amos = combine_excel_sheets(directory_path, target_sheet)
df_amos_b = combine_excel_sheets(directory_path, target_sheet_b)

debug = df_amos.head(20)
debug = df_amos_b.head(20)

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

# Merge AMOS results onto dataframe with chemicals, DTXSID's, raw hazard scores, and water list counts
df8 = pd.merge(df7, df_amos[['DTXSID', 'Sources', 'Patents', 'Articles', 'PubMed Record Count']], how='left', on='DTXSID')
df8 = pd.merge(df8, df_amos_methods_grouped, how='left', on='DTXSID')
df8 = pd.merge(df8, df_amos_fact_sheets_grouped, how='left', on='DTXSID')
df8 = pd.merge(df8, df_amos_spectra_grouped, how='left', on='DTXSID')

debug = df8.head(20)
df8_cols = df8.columns.tolist()

# Read in MS2 files and merge onto SSM chemical features
# First grab ionization mode and merge onto feature ID dataframe
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (20250220)')
df_feature_ion_mode = pd.read_excel('20250220_WW2DW_MS1_analysis_NTA_WebApp_results.xlsx', sheet_name='Final Occurrence Matrix')
df9 = pd.merge(df8, df_feature_ion_mode[['Feature ID', 'Ionization Mode']], how='left', on='Feature ID')

# Grab just SSM chemical features
df9_SSM = df9[df9['SSM chemical feature'] == 'Y']
df9_SSM.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 
df9_SSM = df9_SSM[['Feature ID', 'DTXCID', 'Ionization Mode']]
df9_SSM = df9_SSM.drop_duplicates()

df9_SSM_pos = df9_SSM[df9_SSM['Ionization Mode'] == 'ESI+']
df9_SSM_neg = df9_SSM[df9_SSM['Ionization Mode'] == 'ESI-']

# Read in the MS2 results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/neg')
df_ms2_neg = pd.read_csv('20210131_WW2DW_MS2_neg_CFMID_results_neg.csv')
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/pos')
df_ms2_pos = pd.read_csv('20250131_WW2DW_MS2_pos_CFMID_results_pos.csv')

df9_SSM_pos = pd.merge(df9_SSM_pos, df_ms2_pos[['DTXCID', 'Q-SCORE']], how='left', on='DTXCID')
df9_SSM_neg = pd.merge(df9_SSM_neg, df_ms2_neg[['DTXCID', 'Q-SCORE']], how='left', on='DTXCID')

df9_SSM_both = pd.concat([df9_SSM_pos, df9_SSM_neg], ignore_index=True)

# Merge MS2 results onto dataframe with chemcials, DTXSID's, raw hazard scores, AMOS results, and water list counts
df8.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 
df10 = pd.merge(df8, df9_SSM_both[['Feature ID', 'DTXCID', 'Q-SCORE']], how='left', on=['Feature ID', 'DTXCID'])


debug = df10.head(2000)
debug = df10[df10['Feature ID'] == 200]
debug = df10.drop_duplicates(subset=['DTXSID'])
debug = df10[df10.duplicated(subset=['DTXSID'], keep=False)]
# Duplicate dtxsid: DTXSID001000390


# Export dataframe with all the raw values to CSV
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
#df10.to_csv('WW2DW Data Analysis file 2 (substances, raw data all three legs).csv', sep=',', encoding='utf-8', index=False)


# Grab hazard columns, convert hazard values into numeric values
df10_cols = df10.columns.tolist()
keep_strings = ['DTXSID', '_score', '_authority']
# Filter columns based on the list of strings
df_hazard = df10.filter(regex='|'.join(keep_strings))
df_hazard_b = df_hazard.replace('ND', np.nan)
df_hazard_b = df_hazard_b.replace('Authoritative', 3)
df_hazard_b = df_hazard_b.replace('Screening', 2)
df_hazard_b = df_hazard_b.replace('QSAR Model', 1)
df_hazard_b = df_hazard_b.replace('VH', 4)
df_hazard_b = df_hazard_b.replace('H', 3)
df_hazard_b = df_hazard_b.replace('M', 2)
df_hazard_b = df_hazard_b.replace('L', 1)

debug = df_hazard.head(200)
debug = df_hazard_b.head(2000)

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

debug = df_hazard_b[df_hazard_b['DTXSID'].isna()]
df_hazard_b = df_hazard_b.dropna(subset=['DTXSID'])

# Merge numeric values back onto original results
df11 = df10.drop(columns=score_cols)
df11 = df11.drop(columns=authority_cols)


df12 = pd.merge(df11, df_hazard_b, on='DTXSID', how='left')

debug = df11.head(2000)
debug_b = df12.head(2000)
debug = df11[df11['Feature ID'] == 64]
debug_b = df12[df12['Feature ID'] == 64]

debug = df11[df11['DTXSID'] == 'DTXSID001000390']
debug_b = df12[df12['DTXSID'] == 'DTXSID001000390']

debug = df12[df12['_merge'] != 'both']

df12_cols = df12.columns.tolist()

debug = df12[df12['Q-SCORE'] > 0]

# Export dataframe with all the raw values to CSV
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
df12 = df12.apply(lambda x: pd.to_numeric(x, errors='ignore') if x.dtype == 'object' else x)
#df12.to_csv('WW2DW Data Analysis file 2 (substances, raw data all three legs).csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)


# Collapse onto structure, using sum approach for metadata and max approach for hazard
df13 = df12[df12.columns.drop(list(df12.filter(regex='_score')))]
df13 = df13[df13.columns.drop(list(df13.filter(regex='_authority')))]
df13 = df13[df13.columns.drop(list(df13.filter(regex='Average ')))]
df13 = df13.drop(columns=['Mass', 'Retention Time'])

debug = df13.head(2000)
debug = df13[df13['Feature ID'] == 200]


# Identify columns to fill with 0
columns_to_fill = df13.columns.difference(['SSM chemical feature', 'SSM chemical substance', 'DTXSID', 'Q-SCORE'])

# Fill NaN values in selected columns with 0
df13[columns_to_fill] = df13[columns_to_fill].fillna(0)

df_substance_norm = df13.copy()

# Create a column of the SSM chemical structure True/False so we can eventually collapse into structure
df13['SSM structure'] = df13.groupby(['Feature ID', 'DTXCID'])['SSM chemical substance'].transform('any')


# Now collapse hazard/metadata on structure
# HAZARD: Create a column that has the maximum QAH score for a given structure within a feature
df13['Max_structure_QAH'] = df13.groupby(['Feature ID', 'DTXCID'])['Quality-Adjusted Hazard Score'].transform('max')

debug= df13.drop_duplicates(subset=['DTXSID'], keep='first')
debug_b= df13.drop_duplicates(subset=['DTXCID'], keep='first')

# METADATA: sum across substances for structures
metadata_columns = ['Sources', 'Patents', 'Articles', 'PubMed Record Count', 'AMOS methods count', 'AMOS fact sheets count', 'AMOS spectra count', 'Presence in water lists count']

# Within each structure within each feature, sum the substances
df14 = df13.groupby(['Feature ID', 'DTXCID']).agg({col: 'sum' for col in metadata_columns})
df14 = df14.reset_index()

debug = df13[df13['Feature ID'] == 200]
debug = df14.head(2000)
debug = df14[df14['Feature ID'] == 200]

df14_cols = df14.columns.tolist()

# Within each feature, normalize each metadata value to the maximum value
for col in metadata_columns:
    df14[f'{col}_norm'] = df14.groupby('Feature ID')[col].transform(lambda x: x / x.max())

df14 = df14.fillna(0)
df14['Total_norm'] = df14['Sources_norm'] + df14['Patents_norm'] + df14['Articles_norm'] + df14['PubMed Record Count_norm'] + df14['AMOS methods count_norm'] + df14['AMOS fact sheets count_norm'] + df14['AMOS spectra count_norm'] + df14['Presence in water lists count_norm']

# Get count of structures for each group
df14['structure_count'] = df14.groupby('Feature ID')['DTXCID'].transform('count')

# Merge normalized scores back onto main dataframe
df15 = pd.merge(df13, df14[['Feature ID', 'DTXCID', 'Sources_norm', 'Patents_norm', 'Articles_norm', 'PubMed Record Count_norm', 'AMOS methods count_norm', 
                            'AMOS fact sheets count_norm', 'AMOS spectra count_norm', 'Presence in water lists count_norm', 'Total_norm', 'structure_count']], 
                how='left', on=['Feature ID', 'DTXCID'])


# Export dataframe with all the raw and collapsed values to CSV
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
df15 = df15.apply(lambda x: pd.to_numeric(x, errors='ignore') if x.dtype == 'object' else x)
df15_test = df15.head(2000)
#df15.to_csv('WW2DW Data Analysis file 3 (substances and structures, raw data all three legs).csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)


# Generate rank columns
rank_columns = ['Total_norm', 'Max_structure_QAH', 'Q-SCORE']

for col in rank_columns:
    df15[f'{col}_rank'] = df15.groupby('Feature ID')[col].rank(method='dense', ascending=False)

df15_cols = df15.columns.tolist()

debug = df15[df15['Feature ID'] == 200]
debug = df15.head(2000)

# Pull in occurrence counts/percentages
df16 = pd.merge(df15, df_feature_ion_mode[['Feature ID', 'Final Occurrence Count', 'Final Occurrence Percentage']], how='left', on='Feature ID')

df17 = df16[['Feature ID', 'DTXCID', 'SSM chemical feature', 'SSM structure', 'structure_count', 'Final Occurrence Count', 'Final Occurrence Percentage', 
             'Total_norm', 'Q-SCORE', 'Max_structure_QAH', 'Total_norm_rank', 'Q-SCORE_rank', 'Max_structure_QAH_rank']]

# Export dataframe with all the raw and collapsed values to CSV
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v2')
#df17.to_csv('WW2DW Data Analysis file 4 (Collapsed to structure-level)(Metadata summed and normalized, MS2 normalized, Hazard max value).csv', sep=',', encoding='utf-8', index=False)

df17_test = df17.head(200)

df17 = df17.rename(columns={'Total_norm': 'Metadata score', 'Q-SCORE': 'MS2 score', 'Max_structure_QAH': 'Hazard Score', 'Total_norm_rank':'Metadata rank', 
                            'Q-SCORE_rank':'MS2 rank', 'Max_structure_QAH_rank':'Hazard rank'})


df17['MS2 score'] = df17['MS2 score'].fillna(0)

df18 = df17.drop_duplicates()

# df18.to_csv('WW2DW Data Analysis file 4 (Collapsed to structure-level)(Metadata summed and normalized, MS2 normalized, Hazard max value).csv', 
#             sep=',', float_format="%.15g", quoting=csv.QUOTE_NONE, encoding='utf-8', index=False)

print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))  
