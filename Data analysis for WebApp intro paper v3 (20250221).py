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
import numpy as np
from scipy.stats import spearmanr
import matplotlib.cm as cm
from matplotlib.ticker import LogLocator

time_list = []
start_time = time.time()
time_list.append(time.time())

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (20250221) v3')
# Grab SSM chemicals from tracer summary sheet
df1 = pd.read_excel('WW2DW_MS1_no_cals_20250221_NTA_WebApp_results.xlsx', sheet_name='Tracer Summary')
# Rename DTXSID column of tracer summary to not collide with DTXSID's of DSSTox/CHEM results
df1.rename(columns={'DTXSID':'Tracer DTXSID'}, inplace=True) 

SSM_features = df1['Feature ID'].tolist()

# Grab the WebApp MS1 chemical results
df2 = pd.read_excel('WW2DW_MS1_no_cals_20250221_NTA_WebApp_results.xlsx', sheet_name='Chemical Results')

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
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v4')
#df6.to_csv('WW2DW Data Analysis file 1 (substances, raw hazard values) v4.csv', sep=',', encoding='utf-8', index=False)


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
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v3')
#dtxsid_only.to_csv('WW2DW Data Analysis - DTXSID list all features v3.csv', sep=',', encoding='utf-8', index=False)


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

debug = df8.head(2000)
df8_cols = df8.columns.tolist()

df8_features = df8['Feature ID'].tolist()
df8_features = list(set(df8_features))

debug = list(set(SSM_features) & set(df8_features))

# Read in MS2 files and merge onto SSM chemical features
# First grab ionization mode and merge onto feature ID dataframe
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (20250221) v3')
df_feature_ion_mode = pd.read_excel('WW2DW_MS1_no_cals_20250221_NTA_WebApp_results.xlsx', sheet_name='Final Occurrence Matrix (flags)')
df9 = pd.merge(df8, df_feature_ion_mode[['Feature ID', 'Ionization Mode']], how='left', on='Feature ID')

df8.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 
df9.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 

debug = df9.head(2000)

# # Grab just SSM chemical features
# df9_SSM = df9[df9['SSM chemical feature'] == 'Y']
# df9_SSM.rename(columns={'DTXCID_INDIVIDUAL_COMPONENT':'DTXCID'}, inplace=True) 
# #df9_SSM = df9_SSM[['Feature ID', 'DTXCID', 'Ionization Mode']]
# df9_SSM = df9_SSM[['Feature ID', 'Mass', 'Retention Time', 'DTXCID', 'Ionization Mode']]
# df9_SSM = df9_SSM.drop_duplicates()
# df9_SSM_pos = df9_SSM[df9_SSM['Ionization Mode'] == 'ESI+']
# df9_SSM_neg = df9_SSM[df9_SSM['Ionization Mode'] == 'ESI-']

# Grab all chemical features
df9_all = df9[['Feature ID', 'Mass', 'Retention Time', 'DTXCID', 'Ionization Mode']]
df9_all = df9_all.drop_duplicates()

# Code for MS2 all frag results after MS2 fix 3/20/2025
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250320 All fragmentation, post MS2 fix/neg')
df_ms2_neg = pd.read_csv('WW2DW_all_frag_neg_CFMID_results_neg.csv')
df_ms2_neg_masses = df_ms2_neg['MASS_MGF'].tolist()
df_ms2_neg_masses = list(set(df_ms2_neg_masses))

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250320 All fragmentation, post MS2 fix/pos')
df_ms2_pos = pd.read_csv('all_frag_CFMID_results_pos.csv')
df_ms2_pos_masses = df_ms2_pos['MASS_MGF'].tolist()
df_ms2_pos_masses = list(set(df_ms2_pos_masses))

# # Code for MS2 v3 results (1000 split)
# # Read in the MS2 results - all split into separate files from separate runs
# # These include the SSM results so no need to bring those in
# # First the neg
# os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250318 All fragmentation, split files (1000) v3/neg')
# neg_csv_files = [file for file in os.listdir() if file.endswith('.csv')]
# neg_df_list = [pd.read_csv(file) for file in neg_csv_files]
# df_ms2_neg_old = pd.concat(neg_df_list, ignore_index=True)
# df_ms2_neg_old_masses = df_ms2_neg_old['MASS_MGF'].tolist()
# df_ms2_neg_old_masses = list(set(df_ms2_neg_old_masses))
# # Now do the positive directory
# os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250318 All fragmentation, split files (1000) v3/pos')
# pos_csv_files = [file for file in os.listdir() if file.endswith('.csv')]
# pos_df_list = [pd.read_csv(file) for file in pos_csv_files]
# df_ms2_pos_old = pd.concat(pos_df_list, ignore_index=True)
# df_ms2_pos_old_masses = df_ms2_pos_old['MASS_MGF'].tolist()
# df_ms2_pos_old_masses = list(set(df_ms2_pos_old_masses))
# df_ms2_pos_old_RT = df_ms2_pos_old['RT'].tolist()

# df_ms2_pos_masses_diff = list(set(df_ms2_pos_masses) ^ set(df_ms2_pos_old_masses))
# df_ms2_neg_masses_diff = list(set(df_ms2_neg_masses) ^ set(df_ms2_neg_old_masses))
#### There are some differences between the split files and the single file results, because when the files are split, some of the MS2 features that were de-duplicated in the single file
#### Are kept in the split files as they are in a separate file

# Old code for MS2 v2 results (3000 split)
# # Read in the MS2 results - all split into separate files from separate runs
# # These include the SSM results so no need to bring those in
# # First the neg
# os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250310 All fragmentation, split files/neg')
# neg_csv_files = [file for file in os.listdir() if file.endswith('.csv')]
# neg_df_list = [pd.read_csv(file) for file in neg_csv_files]
# df_ms2_neg = pd.concat(neg_df_list, ignore_index=True)
# # Now do the positive directory
# os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250310 All fragmentation, split files/pos')
# pos_csv_files = [file for file in os.listdir() if file.endswith('.csv')]
# pos_df_list = [pd.read_csv(file) for file in pos_csv_files]
# df_ms2_pos = pd.concat(pos_df_list, ignore_index=True)

# Old code for MS2 single run files
# # Read in the MS2 results - All fragmentation
# os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250221 All fragmentation')
# df_ms2_neg = pd.read_csv('WW2DW_all_fragmentation_neg_CFMID_results_neg.csv')
# df_ms2_pos = pd.read_csv('WW2DW_all_fragmentation_pos_CFMID_results_pos.csv')


# # Read in the MS2 results - Just SSM
# os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/neg')
# df_ms2_neg_SSM = pd.read_csv('20210131_WW2DW_MS2_neg_CFMID_results_neg.csv')
# os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250131 data from Laura results/pos')
# df_ms2_pos_SSM = pd.read_csv('20250131_WW2DW_MS2_pos_CFMID_results_pos.csv')

# Combine MS2 files into one dataframe
df_ms2_neg['Ionization Mode'] = 'ESI-'
df_ms2_pos['Ionization Mode'] = 'ESI+'
# df_ms2_neg_SSM['Ionization Mode'] = 'ESI-'
# df_ms2_pos_SSM['Ionization Mode'] = 'ESI+'
# df_ms2 = pd.concat([df_ms2_pos, df_ms2_neg, df_ms2_pos_SSM, df_ms2_neg_SSM], ignore_index=True)
df_ms2 = pd.concat([df_ms2_pos, df_ms2_neg], ignore_index=True)


# df9_merged['RT diff'] = df9_merged['Retention Time'] - df9_merged['RT']
# df9_merged['mass diff'] = df9_merged['Mass'] - df9_merged['MASS_NEUTRAL']

# df9_SSM_pos = pd.merge(df9_SSM_pos, df_ms2_pos[['DTXCID', 'Q-SCORE', 'MASS_NEUTRAL', 'RT']], how='left', on='DTXCID')
# df9_SSM_neg = pd.merge(df9_SSM_neg, df_ms2_neg[['DTXCID', 'Q-SCORE', 'MASS_NEUTRAL', 'RT']], how='left', on='DTXCID')

# df9_SSM_both = pd.concat([df9_SSM_pos, df9_SSM_neg], ignore_index=True)

# Merge based on DTXCID ond ionization mode.... AFTER determining MS1 and MS2 are within tolerance windows for mass and RT
# Define tolerance thresholds
RT_window = 0.5  # 0.5 min RT window for matching MS1/MS2
Mass_window = 0.05   # 0.05 min RT window for matching MS1/MS2 feature masses
df9_merged = pd.merge(df9_all, df_ms2[['DTXCID', 'Ionization Mode', 'Q-SCORE', 'SUM_SCORE', 'MASS_NEUTRAL', 'RT']], how='inner', on=['DTXCID', 'Ionization Mode'])

debug_features = df9_merged['Feature ID'].tolist()
debug_features = list(set(debug_features))

# Step 2: Filter rows based on the tolerance windows
df9_merged_b = df9_merged[
    (abs(df9_merged['Retention Time'] - df9_merged['RT']) <= RT_window) &
    (abs(df9_merged['Mass'] - df9_merged['MASS_NEUTRAL']) <= Mass_window)
]

df9_merged_b = df9_merged[
    (abs(df9_merged['Retention Time'] - df9_merged['RT']) <= RT_window)
]

# Merge MS2 results onto dataframe with chemicals, DTXSID's, raw hazard scores, AMOS results, and water list counts
df10 = pd.merge(df8, df9_merged_b[['Feature ID', 'Ionization Mode', 'DTXCID', 'Q-SCORE', 'SUM_SCORE']], how='left', on=['Feature ID', 'DTXCID'])

debug = df10.head(2000)
debug = df10[df10['Feature ID'] == 200]
debug = df10[df10['SSM chemical feature'] == 'Y']

debug = df10.groupby('Feature ID').filter(lambda x: x['Q-SCORE'].notna().any())
debug = debug.groupby('Feature ID').filter(lambda x: x['Q-SCORE'].notna().any())
debug_features = debug['Feature ID'].tolist()
debug_features = list(set(debug_features))

# df_test = pd.merge(df8, df9_SSM_both[['Feature ID', 'DTXCID', 'Q-SCORE', 'MASS_NEUTRAL', 'RT']], how='left', on=['Feature ID', 'DTXCID'])
# df_test = df_test[df_test['SSM chemical substance'] == 'Y']
# df_test['RT error'] = df_test['Retention Time'] - df_test['RT']
# df_test['Mass error'] = df_test['Mass'] - df_test['MASS_NEUTRAL']
# debug = df_test[['Feature ID', 'RT error', 'Mass error']]


# #df10 = pd.merge(df8, df9_SSM_both[['Feature ID', 'DTXCID', 'Q-SCORE']], how='left', on=['Feature ID', 'DTXCID'])
# df10 = pd.merge(df8, df9_merged[['Feature ID', 'DTXCID', 'Q-SCORE']], how='left', on=['Feature ID', 'DTXCID'])

debug = df10.head(2000)
debug = df10[df10['Feature ID'] == 200]
debug = df10.drop_duplicates(subset=['DTXSID'])
debug = df10[df10.duplicated(subset=['DTXSID'], keep=False)]

debug = df10[df10['SSM chemical feature'] == 'Y']
debug_b = debug.groupby('Feature ID').filter(lambda x: x['SUM_SCORE'].isna().all())
debug_c = debug_b.drop_duplicates(subset=['Feature ID'], keep='first')
debug_c = debug_c[['Feature ID', 'Ionization Mode', 'Mass', 'Retention Time', 'DTXCID', 'PREFERRED_NAME', 'Q-SCORE', 'SUM_SCORE']]
debug = df10[df10['SSM chemical substance'] == 'Y']

debug_ms2 = df_ms2.drop_duplicates(subset=['ID', 'MASS_MGF', 'RT'], keep='first')

debug = df10[df10['Q-SCORE'] > 0]
debug_feature = debug['Feature ID'].tolist()
debug_feature = list(set(debug_feature))


debug_b = debug[debug['Q-SCORE'] > 0]
debug_b_feature = debug_b['Feature ID'].tolist()
debug_b_feature = list(set(debug_b_feature))
# Duplicate dtxsid: DTXSID001000390


# Export dataframe with all the raw values to CSV
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v4')



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

# debug = df11.head(2000)
# debug_b = df12.head(2000)
# debug = df11[df11['Feature ID'] == 64]
# debug_b = df12[df12['Feature ID'] == 64]

# debug = df11[df11['DTXSID'] == 'DTXSID001000390']
# debug_b = df12[df12['DTXSID'] == 'DTXSID001000390']

# debug = df12[df12['_merge'] != 'both']

df12_cols = df12.columns.tolist()

#debug = df12[df12['Q-SCORE'] > 0]

# Export dataframe with all the raw values to CSV
#os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v3')
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v5')
df12 = df12.apply(lambda x: pd.to_numeric(x, errors='ignore') if x.dtype == 'object' else x)
df12.to_csv('WW2DW Data Analysis file 2 (substances, raw data all three legs) v5.csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)


# Collapse onto structure, using sum approach for metadata and max approach for hazard
df13 = df12[df12.columns.drop(list(df12.filter(regex='_score')))]
df13 = df13[df13.columns.drop(list(df13.filter(regex='_authority')))]
df13 = df13[df13.columns.drop(list(df13.filter(regex='Average ')))]
df13 = df13.drop(columns=['Mass', 'Retention Time'])

debug = df13.head(2000)
debug = df13[df13['Feature ID'] == 200]
#debug.to_csv('Temp debug.csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)

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


# Within each feature, normalize each substance's metadata values to the maximum value
for col in metadata_columns:
    df13[f'Substance_{col}_norm'] = df13.groupby('Feature ID')[col].transform(lambda x: x / x.max())

# Within each structure within each feature, sum the substances
df14 = df13.groupby(['Feature ID', 'DTXCID']).agg({col: 'sum' for col in metadata_columns})
df14 = df14.reset_index()


df14_cols = df14.columns.tolist()

# Within each feature, normalize each structure's metadata values to the maximum value
for col in metadata_columns:
    df14[f'Structure_{col}_norm'] = df14.groupby('Feature ID')[col].transform(lambda x: x / x.max())

debug = df14.head(2000)
debug = df14[df14['Feature ID'] == 68]

df14 = df14.fillna(0)
#df14['Total_norm'] = df14['Sources_norm'] + df14['Patents_norm'] + df14['Articles_norm'] + df14['PubMed Record Count_norm'] + df14['AMOS methods count_norm'] + df14['AMOS fact sheets count_norm'] + df14['AMOS spectra count_norm'] + df14['Presence in water lists count_norm']
metadata_columns_structures = ["Structure_" + col + "_norm" for col in metadata_columns]
df14['Structure_total_norm'] = df14[metadata_columns_structures].sum(axis=1)

metadata_columns_substances = ["Substance_" + col + "_norm" for col in metadata_columns]
df13['Substance_total_norm'] = df13[metadata_columns_substances].sum(axis=1)

# Get count of structures for each group
df14['structure_count'] = df14.groupby('Feature ID')['DTXCID'].transform('count')

# Rename all columns to structure level columns
#df14.columns = ["Structure_count_" + col for col in metadata_columns]
df14.rename(columns={col: "Structure_" + col for col in df14.columns if col in metadata_columns}, inplace=True)

# Rename all columns to substance level columns
df13.rename(columns={col: "Substance_" + col for col in df13.columns if col in metadata_columns}, inplace=True)

debug = df14.head(2000)

debug = df13[df13['Feature ID'] == 200]
debug_b = df14[df14['Feature ID'] == 200]
#debug.to_csv('Temp debug.csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)



# Merge normalized scores back onto main dataframe
# df15 = pd.merge(df13, df14[['Feature ID', 'DTXCID', 'Sources_norm', 'Patents_norm', 'Articles_norm', 'PubMed Record Count_norm', 'AMOS methods count_norm', 
#                             'AMOS fact sheets count_norm', 'AMOS spectra count_norm', 'Presence in water lists count_norm', 'Total_norm', 'structure_count']], 
#                 how='left', on=['Feature ID', 'DTXCID'])
df15 = pd.merge(df13, df14, how='left', on=['Feature ID', 'DTXCID'])


df15 = df15.rename(columns={'Q-SCORE': 'MS2 quotient score', 'SUM_SCORE':'MS2 raw score'})


# Export dataframe with all the raw and collapsed values to CSV
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v3')
df15 = df15.apply(lambda x: pd.to_numeric(x, errors='ignore') if x.dtype == 'object' else x)
df15_test = df15.head(2000)
####df15.to_csv('WW2DW Data Analysis file 3 (Substance and structure level) v3.csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)


# Generate rank columns
rank_columns = ['Structure_total_norm', 'Max_structure_QAH', 'MS2 quotient score']

for col in rank_columns:
    df15[f'{col}_rank'] = df15.groupby('Feature ID')[col].rank(method='dense', ascending=False)

df15_cols = df15.columns.tolist()

debug = df15[df15['Feature ID'] == 200]
debug = df15.head(2000)

blanksub_columns = [col for col in df_feature_ion_mode.columns if col.startswith('BlankSub')]
df_feature_ion_mode['Median blanksub mean feature abundance'] = df_feature_ion_mode[blanksub_columns].median(axis=1)

df_feature_ion_mode = df_feature_ion_mode.rename(columns={'Final Occurrence Count (with flags)': 'Final Occurrence Count', 'Final Occurrence Percentage (with flags)': 'Final Occurrence Percentage'})

# Pull in occurrence counts/percentages
df16 = pd.merge(df15, df_feature_ion_mode[['Feature ID', 'm/z', 'Mass', 'Retention Time', 'Final Occurrence Count', 'Final Occurrence Percentage', 'Median blanksub mean feature abundance']], how='left', on='Feature ID')
debug = df16[df16['Feature ID'] == 200]
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v5')
#df16.to_csv('WW2DW Data Analysis file 3 (Substance and structure level, all results) v5.csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)



df16_SSM = df16[df16['SSM chemical substance'] == 'Y']
df16_SSM = df16_SSM.rename(columns={'PREFERRED_NAME': 'SSM chemical'})

df17 = df16.drop(columns=[col for col in df16.columns if col.startswith('Substance_')])

df17 = df17.drop(['DTXSID', 'SSM chemical substance', 'PREFERRED_NAME', 'Quality-Adjusted Hazard Score', 'Number of end points with data available', 'Completeness Score'], axis=1)    
df17 = pd.merge(df17, df16_SSM[['Feature ID', 'DTXCID', 'SSM chemical']], how='left', on=['Feature ID', 'DTXCID'])
    
debug = df17[df17['Feature ID'] == 200]
df17_cols = df17.columns.tolist()

df17 = df17.rename(columns={'Total_norm': 'Metadata score', 'Q-SCORE': 'MS2 score', 'Max_structure_QAH': 'Hazard Score', 'Total_norm_rank':'Metadata rank'})

df18 = df17.drop_duplicates(subset=['Feature ID', 'DTXCID'])



# Need to merge completeness score from file 3 onto file 4
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v4')
#df_file_3 = pd.read_csv('WW2DW Data Analysis file 3 (Substance and structure level, all results) v3 (20250226).csv')
#df_file_4 = pd.read_csv('WW2DW Data Analysis file 4 (Structure level only, all results) v3 (20250226).csv')
df_file_3 = pd.read_csv('WW2DW Data Analysis file 3 (Substance and structure level, all results) v4.csv')
df_file_4 = df18.copy()
df_file_3_cols = df_file_3.columns.tolist()
df_file_4_cols = df_file_4.columns.tolist()

debug = df_file_3[df_file_3['Feature ID'] == 200]

# 2/27/2025 - File 4 is missing "Hazard Completeness Score" associated with each hazard score
# Sort file 3 by Feature ID, DTXCID, hazard score, completeness score
# Then de-duplicate by Feature ID, DTXCID - Should keep the DTXCID row with the highest hazard score
# Then merge completeness score onto file 4
df_file_3 = df_file_3.sort_values(by=['Feature ID', 'DTXCID', 'Quality-Adjusted Hazard Score', 'Completeness Score'], ascending=[True, True, False, False])
debug = df_file_3[df_file_3['Feature ID'] == 200]
df_file_3 = df_file_3.rename(columns={'Completeness Score': 'Hazard Completeness Score'})

df_file_3_b = df_file_3.drop_duplicates(subset=['Feature ID', 'DTXCID'])
debug = df_file_3_b[df_file_3_b['Feature ID'] == 200]

df_file_export = pd.merge(df_file_4, df_file_3_b[['Feature ID', 'DTXCID', 'Hazard Completeness Score']], how='left', on=['Feature ID', 'DTXCID'])


# Read in level 1 ID's from Heather's file in order to annotate features as level 1's
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/')
df_IDs = pd.read_csv('List of SSM and level 1 IDs from WW2DW SI and Heather.csv', encoding='cp1252')
df_IDs_DTXSID = df_IDs['Final DTXSID'].tolist()
df_IDs_DTXSID = list(set(df_IDs_DTXSID))

df_IDs = df_IDs.rename(columns={'Final DTXSID': 'DTXSID', 'RT':'ID Retention Time', 'm/z': 'ID m/z'})
df_IDs = df_IDs.drop(columns=['Adduct', 'MS1 hit name', 'XIC review', 'MS2 review'])

# Read in just the level 1's from the WW2DW SI:
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/')
df_level1 = pd.read_csv('WW2DW Level 1s with DSSTox identifiers.csv')

# Merge on RT's/masses from Heather's file
df_level1 = pd.merge(df_level1, df_IDs[['DTXSID', 'ID m/z', 'ID Retention Time', 'PREFERRED_NAME']], how='left', on=['DTXSID', 'PREFERRED_NAME'])
df_level1['match'] = df_level1['DTXCID'].isin(df_file_export['DTXCID'])

# Merge on level 1/SSMs
RT_window = 0.2  # 0.2 min RT window for matching level 1's/SSMs
Mass_window = 0.05   # 0.05 Da mass window for for matching level 1's/SSMs
df19 = pd.merge(df_file_export[['Feature ID', 'DTXCID', 'm/z', 'Retention Time']], df_level1[['DTXCID', 'PREFERRED_NAME', 'ID m/z', 'ID Retention Time']], how='inner', on=['DTXCID'])

debug = df_file_export[df_file_export['SSM structure'] == True]

# Step 2: Filter merged rows based on the tolerance windows
df20 = df19[
    (abs(df19['Retention Time'] - df19['ID Retention Time']) <= RT_window) &
    (abs(df19['m/z'] - df19['ID m/z']) <= Mass_window)
]

df20['Level 1 structure'] = 'Y'
df20 = df20.rename(columns={'PREFERRED_NAME':'Level 1 chemical'})

# Merge MS2 results onto dataframe with chemicals, DTXSID's, raw hazard scores, AMOS results, and water list counts
df21 = pd.merge(df_file_export, df20[['Feature ID', 'DTXCID', 'Level 1 structure', 'Level 1 chemical']], how='left', on=['Feature ID', 'DTXCID'])

debug = df21[df21['Level 1 structure'] == 'Y']
debug = df21[df21['Feature ID'] == 407]

# Identify features associated with level 1's
df21['Level 1 feature'] = df21.groupby('Feature ID')['Level 1 structure'].transform(lambda x: 'Y' in x.values)


df21['SSM or level 1 feature'] = df21['Level 1 feature'] | df21['SSM chemical feature']
df21['SSM or level 1 structure'] = df21['SSM structure'] | (df21['Level 1 structure'] == 'Y')
df21['SSM or level 1 name'] = df21['SSM chemical'].combine_first(df21['Level 1 chemical'])


debug = df21[df21['SSM or level 1 feature'] == True]
debug = df21[df21['SSM or level 1 structure'] == True]

# Grab the ionization mode for the features from the MS1 results and merge onto df21 since it's missing some values
# Grab the WebApp MS1 chemical results
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS1 WebApp results/MS1 WebApp results (20250221) v3')
df_feature_pos = pd.read_excel('WW2DW_MS1_no_cals_20250221_NTA_WebApp_results.xlsx', sheet_name='All Detection Statistics (Pos)')
df_feature_neg = pd.read_excel('WW2DW_MS1_no_cals_20250221_NTA_WebApp_results.xlsx', sheet_name='All Detection Statistics (Neg)')
df_feature_pos = df_feature_pos[['Feature ID', 'Ionization Mode']]
df_feature_neg = df_feature_neg[['Feature ID', 'Ionization Mode']]
df_feature_all = pd.concat([df_feature_pos, df_feature_neg], ignore_index=True)

df21.rename(columns={'Ionization Mode':'Ionization Mode (old)'}, inplace=True) 

df22 = pd.merge(df21, df_feature_all, how='left', on='Feature ID')

debug = df22[['Feature ID', 'Ionization Mode', 'Ionization Mode (old)']]
debug_b = debug[debug['Ionization Mode'] != debug['Ionization Mode (old)']]
debug_b = debug_b.drop_duplicates(subset=['Feature ID'])

df22 = df22.drop(columns='Ionization Mode (old)')

os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v5')
# df22.to_csv('WW2DW Data Analysis file 4 (Structure level only, all results) v5.csv', 
#             sep=',', float_format="%.15g", encoding='utf-8', index=False)




# Modify the source input file for Safia's visualizations
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v5')

df_file_4 = pd.read_csv('WW2DW Data Analysis file 4 (Structure level only, all results) v5 (formatted).csv')
df_file_2 = pd.read_csv('WW2DW Data Analysis file 2 (substances, raw data all three legs) v5.csv')

# De-duplicate file 2, which has all the converted hazard/authority scores
df_file_2_b = df_file_2.drop_duplicates(subset='DTXCID')
df_file_2_c = df_file_2_b.loc[:, df_file_2_b.columns.str.endswith('DTXCID') | df_file_2_b.columns.str.endswith('_score') | df_file_2_b.columns.str.endswith('_authority')]

df_file_5 = pd.merge(df_file_4, df_file_2_c, how='left', on='DTXCID')
df_file_5_columns = df_file_5.columns.tolist()
# Merge on the converted scores to file 4

df_file_5 = df_file_5.apply(lambda x: pd.to_numeric(x, errors='ignore') if x.dtype == 'object' else x)
#df_file_5.to_csv('WW2DW Data Analysis file 5 (for bar chart visuals) v5.csv', sep=',', float_format="%.15g", encoding='utf-8', index=False)








################## 
# Everything below is getting tallies for the paper using df21 which can also just be read in from the export above
# So start with this step in the future for doing the paper analysis tallies
##################
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v5')
df22 = pd.read_csv('WW2DW Data Analysis file 4 (Structure level only, all results) v5 (formatted).csv')

# Rename columns
df22 = df22.rename(columns={'Structure_Patents':'Patents', 'Structure_Articles':'Articles', 'Structure_PubMed Record Count':'PubMed Records', 'Structure_Sources':'Sources',
                            'Structure_AMOS spectra count':'AMOS spectra', 'Structure_AMOS methods count':'AMOS methods', 
                            'Structure_AMOS fact sheets count':'AMOS fact sheets', 'Structure_Presence in water lists count':'Water lists', 
                            'Structure_total_norm':'Metadata Score', 'MS2 quotient score':'MS2 Score'})


# Merge in specific hazard endpoint data
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/Data analysis/v5')

# Need to pull in and merge the endpoint specific data
df_hazard = pd.read_csv('WW2DW Data Analysis file 2 (substances, raw data all three legs) v5.csv')

# De-duplicate hazard data on structure (keeping highest QAH, and completeness score structures)
df_hazard = df_hazard.sort_values(by=['Feature ID', 'DTXCID', 'Quality-Adjusted Hazard Score', 'Completeness Score'], ascending=[True, True, False, False])
df_hazard = df_hazard.drop_duplicates(subset=['Feature ID', 'DTXCID'], keep='first')


# Identify authority columns and their prefixes
authority_columns = [col for col in df_hazard.columns if col.endswith('_authority')]
prefixes = [col.rsplit('_', 1)[0] for col in authority_columns]

# Select columns with the same prefixes
hazard_columns = [col for col in df_hazard.columns if any(col.startswith(prefix) for prefix in prefixes)]
hazard_columns_b = ['Feature ID', 'DTXCID'] + hazard_columns

df22 = pd.merge(df22, df_hazard[hazard_columns_b], how='left', on=['Feature ID', 'DTXCID'])



# Get total features count
total_feature = df22['Feature ID'].tolist()
df22_col = df22.columns.tolist()
total_feature = list(set(total_feature))

# Get counts for each ionization mode
df22_feature = df22.drop_duplicates(subset=['Feature ID'])
df_feature_ion_counts = df22_feature['Ionization Mode'].value_counts()

# Get counts of total features w/ MS2 data
def has_non_nan(group):
    return group.notna().any()

grouped = df22.groupby('Feature ID')['MS2 Score'].apply(has_non_nan)

total_feature_has_MS2 = grouped.sum()
total_feature_no_MS2 = len(grouped) - total_feature_has_MS2


# Get counts of level 1/SSM chemicals
ID_chemicals = df22['SSM or level 1 name'].tolist()
ID_chemicals = list(set(ID_chemicals)) # Subtract 1 because of nan

# Get counts of level 1/SSM features
ID_results = df22[df22['SSM or level 1 feature'] == True]
ID_feature = ID_results.drop_duplicates(subset=['Feature ID'])
ID_feature_ion_counts = ID_feature['Ionization Mode'].value_counts()

# Get counts of level1/SSM features in different ionization modes
def mode_status(modes):
    unique_modes = set(modes)
    if 'ESI+' in unique_modes and 'ESI-' in unique_modes:
        return 'both'
    elif 'ESI+' in unique_modes:
        return 'ESI+'
    elif 'ESI-' in unique_modes:
        return 'ESI-'
    return 'none'

ID_feature_modes = ID_results.groupby('SSM or level 1 name')['Ionization Mode'].apply(mode_status)
ID_feature_modes_counts = ID_feature_modes.value_counts()

# Get counts of level1/SSM features w/ MS2 data
grouped = ID_results.groupby('Feature ID')['MS2 Score'].apply(has_non_nan)

ID_feature_has_MS2 = grouped.sum()
ID_feature_no_MS2 = len(grouped) - ID_feature_has_MS2

# Get counts of level1/SSM chemicals w/ MS2 data
ID_results['has_MS2'] = ID_results.groupby('Feature ID')['MS2 Score'].transform(lambda group: group.notna().any())


# Distribution plots of counts of structures
plt.figure(figsize=(8, 6))
sns.violinplot(y=df22_feature['structure_count'], color='skyblue')
plt.title('Chemical structures retrieved from database per feature')
plt.ylabel('Number of chemical structures retrieved')
plt.grid(axis='x', linestyle='--', alpha=0.7)
plt.show()

# Define number of bins
num_bins = 10

# Plot histogram
plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(df22_feature['structure_count'], bins=num_bins, color='skyblue', edgecolor='black')

# Calculate percentages
percentages = (counts / counts.sum()) * 100

# Annotate percentages on the histogram
for count, percentage, patch in zip(counts, percentages, patches):
    plt.annotate(f'{percentage:.1f}%', xy=(patch.get_x() + patch.get_width() / 2, count),
                 xytext=(0, 5), textcoords='offset points', ha='center', va='bottom', fontsize=10)

# Add labels and title
plt.title('Chemical structures retrieved from database per feature')
plt.xlabel('Number of chemical structures retrieved')
plt.ylabel('Frequency')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show plot
plt.show()

# Histogram for number of structures per feature: All features
# Define bin size
bin_size = 50

# Calculate bins based on bin size
min_value = df22_feature['structure_count'].min()
max_value = df22_feature['structure_count'].max()
bins = np.arange(min_value, max_value + bin_size, bin_size)

# Plot histogram
plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(df22_feature['structure_count'], bins=bins, color='skyblue', edgecolor='black')

# Calculate percentages
percentages = (counts / counts.sum()) * 100

# Annotate percentages on the histogram
for count, percentage, patch in zip(counts, percentages, patches):
    # Calculate the center of the patch
    x_center = patch.get_x() + patch.get_width() / 2 + 2
    plt.annotate(f'{percentage:.1f}%', xy=(x_center, count),
                 xytext=(0, 5), textcoords='offset points', ha='center', va='bottom', fontsize=10)


# Adjust y-axis limit to provide extra space for annotations
plt.ylim(0, max(counts) * 1.1)

# Add labels and title
plt.title('Chemical structures retrieved from database per feature: All features')
plt.xlabel('Number of chemical structures retrieved')
plt.ylabel('Frequency')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show plot
plt.show()



# Histogram for number of structures per feature: Known features
# Define bin size
bin_size = 50

ID_results_combined_features = ID_results_combined.drop_duplicates(subset='Feature ID') # ID_results_combined initialized below

debug = ID_results_combined_features['structure_count'].tolist()

# Calculate bins based on bin size - Use all features bins
#min_value = ID_results_combined_features['structure_count'].min()
#max_value = ID_results_combined_features['structure_count'].max()
#bins = np.arange(min_value, max_value + bin_size, bin_size)

# Plot histogram
plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(ID_results_combined_features['structure_count'], bins=bins, color='skyblue', edgecolor='black')

# Calculate percentages
percentages = (counts / counts.sum()) * 100

# Annotate percentages on the histogram
for count, percentage, patch in zip(counts, percentages, patches):
    # Calculate the center of the patch
    x_center = patch.get_x() + patch.get_width() / 2 + 2
    plt.annotate(f'{percentage:.1f}%', xy=(x_center, count),
                 xytext=(0, 5), textcoords='offset points', ha='center', va='bottom', fontsize=10)


# Adjust y-axis limit to provide extra space for annotations
plt.ylim(0, max(counts) * 1.1)

# Add labels and title
#plt.title('Chemical structures retrieved from database per feature: $\textit{Known chemical features}$')
plt.title('Chemical structures retrieved from database per feature: $\it{Known\ chemical\ features}$')
plt.xlabel('Number of chemical structures retrieved')
plt.ylabel('Frequency')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show plot
plt.show()



# # Plot point cloud using scatter plot
# jitter = np.random.uniform(-0.2, 0.2, size=len(df22_feature))
# plt.figure(figsize=(10, 6))
# plt.scatter(df22_feature['structure_count'], jitter, alpha=0.4, color='blue', s=15)
# plt.title('Chemical structures retrieved from database per feature')
# plt.xlabel('Number of chemical structures retrieved')
# #plt.ylabel('Jitter')
# plt.grid(axis='x', linestyle='--', alpha=0.7)
# plt.yticks([])  # Hide y-axis ticks as jitter is just for visualization
# plt.show()




################## 
# 3.2 Descriptive Statistics for Retrieved Metadata 
##################

non_ID_results = df22[df22['SSM or level 1 feature'] == False]
non_ID_feature = non_ID_results.drop_duplicates(subset=['Feature ID'])

ID_results_cols = ID_results.columns.tolist()

# Summary stats for each metadata field
ID_results_true = ID_results[ID_results['SSM or level 1 structure'] == True]
ID_results_other = ID_results[ID_results['SSM or level 1 structure'] == False]
ID_results_true['source'] = 'Correct Candidate'
ID_results_other['source'] = 'Incorrect Candidate'
ID_results_combined = pd.concat([ID_results_true, ID_results_other])






# metadata columns
# metadata_columns = ['Structure_Sources', 'Structure_Patents', 'Structure_Articles', 'Structure_PubMed Record Count', 'Structure_AMOS methods count', 
#                     'Structure_AMOS fact sheets count', 'Structure_AMOS spectra count', 'Structure_Presence in water lists count']

# metadata_columns_b = ['Structure_Sources', 'Structure_Patents', 'Structure_Articles', 'Structure_PubMed Record Count', 'Structure_AMOS methods count', 
#                     'Structure_AMOS fact sheets count', 'Structure_AMOS spectra count', 'Structure_Presence in water lists count', 'MS2 raw score', 
#                     'MS2 quotient score', 'Hazard Score', 'Hazard Completeness Score']


metadata_columns = ['Patents', 'Articles', 'PubMed Records', 'Sources', 'AMOS spectra', 'AMOS methods', 
                    'AMOS fact sheets', 'Water lists', 'Metadata Score']

metadata_columns_b = ['Patents', 'Articles', 'PubMed Records', 'Sources', 'AMOS spectra', 'AMOS methods', 
                    'AMOS fact sheets', 'Water lists', 'Metadata Score', 'MS2 raw score', 'MS2 Score', 'Hazard Score', 'Hazard Completeness Score']

# metadata columns: normalized
metadata_columns_norm = ['Structure_Sources_norm', 'Structure_Patents_norm', 'Structure_Articles_norm', 'Structure_PubMed Record Count_norm', 
                         'Structure_AMOS methods count_norm', 'Structure_AMOS fact sheets count_norm', 'Structure_AMOS spectra count_norm', 
                         'Structure_Presence in water lists count_norm']

# Calculate average, minimum, and maximum for each column in the list
ID_results_true_stats_norm = ID_results_true[metadata_columns_norm].agg(['mean', 'min', 'max'])
ID_results_other_stats_norm = ID_results_other[metadata_columns_norm].agg(['mean', 'min', 'max'])


# Calculate average, minimum, and maximum for each column in the list
ID_results_true_stats = ID_results_true[metadata_columns].agg(['mean', 'min', 'max'])
ID_results_other_stats = ID_results_other[metadata_columns].agg(['mean', 'min', 'max'])


# Calculate median and IQR for each column - True candidate
results = {}
for column in metadata_columns_b:
    median = ID_results_true[column].median()
    q1 = ID_results_true[column].quantile(0.25)
    q3 = ID_results_true[column].quantile(0.75)
    iqr = q3 - q1
    mean = ID_results_true[column].mean()
    std_dev = ID_results_true[column].std()
    results[column] = {'Mean': mean, 'Standard Deviation': std_dev, 'Median': median, 'IQR': iqr}

# Convert results to DataFrame for better visualization
ID_results_true_stats_b = pd.DataFrame(results).T

# Calculate median and IQR for each column - Other candidates
results = {}
for column in metadata_columns_b:
    median = ID_results_other[column].median()
    q1 = ID_results_other[column].quantile(0.25)
    q3 = ID_results_other[column].quantile(0.75)
    iqr = q3 - q1
    mean = ID_results_other[column].mean()
    std_dev = ID_results_other[column].std()
    results[column] = {'Mean': mean, 'Standard Deviation': std_dev, 'Median': median, 'IQR': iqr}

# Convert results to DataFrame for better visualization
ID_results_other_stats_b = pd.DataFrame(results).T


# Normalized value calculations
# Calculate median and IQR for each column - True candidate
results = {}
for column in metadata_columns_norm:
    median = ID_results_true[column].median()
    q1 = ID_results_true[column].quantile(0.25)
    q3 = ID_results_true[column].quantile(0.75)
    iqr = q3 - q1
    mean = ID_results_true[column].mean()
    std_dev = ID_results_true[column].std()
    results[column] = {'Mean': mean, 'Standard Deviation': std_dev, 'Median': median, 'IQR': iqr}

# Convert results to DataFrame for better visualization
ID_results_true_stats_norm_b = pd.DataFrame(results).T

# Calculate median and IQR for each column - Other candidates
results = {}
for column in metadata_columns_norm:
    median = ID_results_other[column].median()
    q1 = ID_results_other[column].quantile(0.25)
    q3 = ID_results_other[column].quantile(0.75)
    iqr = q3 - q1
    mean = ID_results_other[column].mean()
    std_dev = ID_results_other[column].std()
    results[column] = {'Mean': mean, 'Standard Deviation': std_dev, 'Median': median, 'IQR': iqr}

# Convert results to DataFrame for better visualization
ID_results_other_stats_norm_b = pd.DataFrame(results).T


all_columns = ['Sources', 'Patents', 'Articles', 'PubMed Records', 'AMOS methods', 
                    'AMOS fact sheets', 'AMOS spectra', 'Water lists', 'MS2 Score', 'Hazard Score']

all_columns_b = ['Patents', 'Articles', 'PubMed Records', 'Sources', 'AMOS spectra', 'AMOS methods',  
                    'AMOS fact sheets',  'Water lists', 'Metadata Score', 'MS2 Score', 'Hazard Score']


# Impute zero values to minimum/sqrt(2) within each feature


# Function to impute zero values
def impute_zeros(group, columns):
    for col in columns:
        # Extract non-zero values
        non_zero_values = group[col][group[col] > 0]
        
        # Check for the condition where no imputation is needed
        if len(non_zero_values) == 1 and group.loc[non_zero_values.index[0], 'SSM or level 1 structure']:
            continue
        
        # Calculate the minimum non-zero value divided by the square root of 2
        if not non_zero_values.empty:
            min_value = non_zero_values.min()
            imputed_value = min_value / np.sqrt(2)
        else:
            imputed_value = 1  # Impute 1 if all values are zero
        
        # Replace zero values with the imputed value
        group[col] = group[col].replace(0, imputed_value)
    return group


ID_results_combined_imputed = ID_results_combined.groupby('Feature ID').apply(impute_zeros, columns=all_columns_b)



# # Non-imputed
# # Melt the DataFrame to long format for use with seaborn
# df_melted = ID_results_combined.melt(id_vars='source', value_vars=all_columns_b, var_name='value_type', value_name='value')

# plot_font_size = 12

# # Define a custom color palette
# custom_palette = {
#     'Correct Candidate': '#1f77b4',  # Blue
#     'Incorrect Candidate': '#ff7f0e'   # Orange
# }
# # Set up the figure
# plt.figure(figsize=(14, 8))

# # Create a strip plot with increased jitter and decreased transparency
# sns.stripplot(data=df_melted, x='value_type', y='value', hue='source', jitter=0.2, dodge=True, palette=custom_palette, alpha=0.7, size=3)

# # Set y-axis to log scale
# plt.yscale('log')

# # Customize the plot
# plt.title('Comparison of chemical results for level 1/SSM chemical features', fontsize=plot_font_size)
# #plt.xlabel('Result type')
# plt.ylabel('Counts or Scores', fontsize=plot_font_size)
# plt.xticks(rotation=45, ha='right', fontsize=plot_font_size)  # Rotate x-axis labels for better readability
# #plt.legend(title='Source', bbox_to_anchor=(1.05, 1), loc='upper left')

# # Place the legend inside the plot
# plt.legend(title='Source', loc='upper right', fontsize=plot_font_size, title_fontsize=plot_font_size+2)

# plt.tight_layout()
# plt.show()






# Example DataFrame setup
# Assuming ID_results_combined_imputed and all_columns_b are defined elsewhere

# Melt the DataFrame to long format for use with seaborn
df_melted = ID_results_combined_imputed.melt(id_vars='source', value_vars=all_columns_b, var_name='value_type', value_name='value')

plot_font_size = 16

# Define a custom color palette
custom_palette = {
    'Correct Candidate': '#1f77b4',  # Blue
    'Incorrect Candidate': '#ff7f0e'   # Orange
}

# Set up the figure
plt.figure(figsize=(14, 10))

# Create a strip plot with increased jitter and decreased transparency
sns.stripplot(data=df_melted, x='value_type', y='value', hue='source', jitter=0.2, dodge=True, palette=custom_palette, alpha=0.7, size=3)

# Set y-axis to log scale
plt.yscale('log')

# Customize the plot
#plt.title('Counts of Metadata, MS$^2$ and Hazard Scores for Known Chemical Features', fontsize=plot_font_size+4)
plt.title('Counts of Metadata, MS$^2$ and Hazard Scores for $\it{Known\ Chemical\ Features}$', fontsize=plot_font_size+4)
plt.xlabel(None)  # Remove x-axis label
plt.ylabel('Counts or Scores', fontsize=plot_font_size)
plt.xticks(rotation=45, ha='right', fontsize=plot_font_size)  # Rotate x-axis labels for better readability
plt.yticks(fontsize=plot_font_size)

# Access x-axis tick labels and modify "MS2 Score" to have a bold superscript "2"
ax = plt.gca()
current_labels = [label.get_text() for label in ax.get_xticklabels()]
new_labels = [label if label != "MS2 Score" else r"MS$\mathbf{^2}$ Score" for label in current_labels]
new_labels = [label if label != "Patents" else "PubChem Patents" for label in new_labels]
new_labels = [label if label != "Articles" else "PubChem Articles" for label in new_labels]
new_labels = [label if label != "PubMed Records" else "PubMed Articles" for label in new_labels]
new_labels = [label if label != "Sources" else "PubChem Sources" for label in new_labels]
new_labels = [label if label != "Water lists" else "Dashboard Water Lists" for label in new_labels]

ax.set_xticklabels(new_labels, fontsize=plot_font_size)

# Access x-axis tick labels and make the last three bold
for label in ax.get_xticklabels()[-3:]:
    label.set_fontweight('bold')

# Configure y-axis ticks
ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=100))
ax.yaxis.set_minor_formatter(plt.NullFormatter())  # Do not label minor ticks

# Set consistent font size for y-axis tick labels
ax.tick_params(axis='y', labelsize=plot_font_size)

# Place the legend inside the plot
plt.legend(loc='upper right', fontsize=plot_font_size, title_fontsize=plot_font_size+2)

plt.tight_layout()
plt.show()







# # Plot distribution of total metadata scores
# plt.figure(figsize=(10, 6))
# sns.boxplot(x='source', y='Structure_total_norm', data=ID_results_combined, palette='pastel')
# plt.title("Comparison of normalized metadata scores for Level 1/SSM chemicals")
# plt.xlabel("Candidate category")
# plt.ylabel('Normalized metadata scores')
# plt.grid(axis='y', linestyle='--', alpha=0.7)
# #plt.yscale('log')
# plt.show()


# Calculate ratio columns
# # Function to calculate ratios
# def calculate_ratios(group):
#     # Find the values where 'is_true' is True
#     true_values = group.loc[group['SSM or level 1 structure'] == True, all_columns_b]
    
#     # If there are no True values, return the original group (to avoid division by zero)
#     if true_values.empty:
#         return group
    
#     # Calculate the ratio for each specified column
#     for column in all_columns_b:
#         ratio_column_name = f'{column}_ratio'
#         group[ratio_column_name] = group[column] / true_values[column].values[0]
    
#     return group

# # Apply the function to each group
# ID_results_combined_ratios = ID_results_combined.groupby('Feature ID').apply(calculate_ratios).reset_index(drop=True)


# Function to calculate reversed ratios
def calculate_reversed_ratios(group):
    # Find the values where 'is_true' is True
    true_values = group.loc[group['SSM or level 1 structure'] == True, all_columns_b]
    
    # If there are no True values, return the original group (to avoid division by zero)
    if true_values.empty:
        return group
    
    # Calculate the reversed ratio for each specified column
    for column in all_columns_b:
        ratio_column_name = f'{column}_ratio'
        group[ratio_column_name] = true_values[column].values[0] / group[column]
    
    return group

# Apply the function to each group
# Reset the index to remove 'Feature ID' from being an index
ID_results_combined_imputed = ID_results_combined_imputed.drop(columns='Feature ID')
ID_results_combined_imputed = ID_results_combined_imputed.reset_index()
ID_results_combined_imputed_ratios = ID_results_combined_imputed.groupby('Feature ID').apply(calculate_reversed_ratios).reset_index(drop=True)

Ratios_to_plot = ID_results_combined_imputed_ratios[ID_results_combined_imputed_ratios['SSM or level 1 structure'] == False]
all_columns_ratios = []
for col in all_columns_b:
    all_columns_ratios.append(str(col) + '_ratio')




# Calculate percentages of each field ratio > 1
# Function to calculate the percentage of non-inf values greater than 1
def percentage_greater_than_one(df, columns):
    percentages = {}
    for column in columns:
        # Filter out infinite values
        non_inf_values = df[column].replace([np.inf, -np.inf], np.nan).dropna()
        
        # Count values greater than 1
        count_greater_than_one = (non_inf_values > 1).sum()
        
        # Calculate percentage
        percentage = (count_greater_than_one / len(non_inf_values)) * 100
        percentages[column] = percentage
    
    return percentages


def percentage_less_than_one(df, columns):
    percentages = {}
    for column in columns:
        # Filter out infinite values
        non_inf_values = df[column].replace([np.inf, -np.inf], np.nan).dropna()
        
        # Count values greater than 1
        count_less_than_one = (non_inf_values < 1).sum()
        
        # Calculate percentage
        percentage = (count_less_than_one / len(non_inf_values)) * 100
        percentages[column] = percentage
    
    return percentages

def percentage_equal_to_one(df, columns):
    percentages = {}
    for column in columns:
        # Filter out infinite values
        non_inf_values = df[column].replace([np.inf, -np.inf], np.nan).dropna()
        
        # Count values greater than 1
        count_equal_to_one = (non_inf_values == 1).sum()
        
        # Calculate percentage
        percentage = (count_equal_to_one / len(non_inf_values)) * 100
        percentages[column] = percentage
    
    return percentages


# Calculate the percentages
ratios_greater_than_one = percentage_greater_than_one(ID_results_combined_imputed_ratios, all_columns_ratios)
ratios_less_than_one = percentage_less_than_one(ID_results_combined_imputed_ratios, all_columns_ratios)
ratios_equal_to_one = percentage_equal_to_one(ID_results_combined_imputed_ratios, all_columns_ratios)


# Imputed
# Melt the DataFrame to long format for use with seaborn
df_melted = Ratios_to_plot.melt(id_vars='source', value_vars=all_columns_ratios, var_name='value_type', value_name='value')

plot_font_size = 16

# Define a custom color palette
custom_palette = {
    'Correct Candidate': '#ff7f0e',
    'Incorrect Candidate': '#1f77b4'
}

# Set up the figure
plt.figure(figsize=(14, 9))

# Create a strip plot with increased jitter and decreased transparency
sns.stripplot(data=df_melted, x='value_type', y='value', hue='source', jitter=0.2, dodge=True, palette=custom_palette, alpha=0.2, size=2, legend=False)

# Remove the suffix '_suffix' from x-axis labels
suffix_to_remove = '_ratio'
current_labels = plt.gca().get_xticklabels()
new_labels = [label.get_text().replace(suffix_to_remove, '') for label in current_labels]

# Update the specific label "MS2 Score" to have a superscript "2"
new_labels = [label if label != "MS2 Score" else r"MS$\mathbf{^2}$ Score" for label in new_labels]
new_labels = [label if label != "Patents" else "PubChem Patents" for label in new_labels]
new_labels = [label if label != "Articles" else "PubChem Articles" for label in new_labels]
new_labels = [label if label != "PubMed Records" else "PubMed Articles" for label in new_labels]
new_labels = [label if label != "Sources" else "PubChem Sources" for label in new_labels]
new_labels = [label if label != "Water lists" else "Dashboard Water Lists" for label in new_labels]

plt.gca().set_xticklabels(new_labels)

# Draw axhline with higher zorder
#.axhline(y=1, color='red', linestyle='--', linewidth=3, zorder=3)
#plt.axhline(y=1, color='red', linestyle='--', linewidth=3, zorder=3)
plt.axhline(y=1, color='black', linestyle='-', linewidth=2, zorder=3)

# Set y-axis to log scale
plt.yscale('log')

# Customize the plot
#plt.title('Ratios of Values for Correct Candidates vs. Incorrect Candidates of Known Chemical Features', fontsize=plot_font_size+4)
plt.title('Ratios of Values for Correct Candidates vs. Incorrect Candidates of $\it{Known\ Chemical\ Features}$', fontsize=plot_font_size+4)
plt.xlabel(None)  # Remove x-axis label
plt.ylabel('Ratios', fontsize=plot_font_size)
plt.yticks(fontsize=plot_font_size) 
plt.xticks(rotation=45, ha='right', fontsize=plot_font_size)  # Rotate x-axis labels for better readability
plt.xticks(fontsize=plot_font_size)

# Extend the y-axis maximum to create more space at the top
# Increase the maximum y limit by a certain factor
current_ylim = plt.ylim()  # Get current y-axis limits
plt.ylim(current_ylim[0] * 50, current_ylim[1] * 5)  # Increase the upper limit by 0.5

current_xlim = plt.xlim()  # Get current y-axis limits
plt.xlim(current_xlim[0]-2, current_xlim[1])  # Increase the upper limit by 0.5

# Access x-axis tick labels and make the last three bold
ax = plt.gca()
for label in ax.get_xticklabels()[-3:]:  
    label.set_fontweight('bold')


# Annotate each strip plot with the percentage
for i, column in enumerate(all_columns_ratios):
    # Calculate the maximum y-value for the current column
    max_value = df_melted[df_melted['value_type'] == column]['value'].replace([np.inf, -np.inf], np.nan).max()
    min_value = df_melted[df_melted['value_type'] == column]['value'].replace([np.inf, -np.inf], np.nan).min()
    # Annotate slightly above the maximum value
    percentage_text_greater = f"{ratios_greater_than_one[column]:.1f}%"
    percentage_text_equal = f"{ratios_equal_to_one[column]:.1f}%"
    percentage_text_lesser = f"{ratios_less_than_one[column]:.1f}%"
    #plt.annotate(percentage_text_greater, xy=(i-0.4, max_value * 2), ha='center', fontsize=plot_font_size, color='black', rotation=90)  # Adjust multiplier as needed
    plt.annotate(percentage_text_greater, xy=(i, max_value * 2), ha='center', fontsize=plot_font_size, color='black')  # Adjust multiplier as needed
    plt.annotate(percentage_text_equal, xy=(i-0.35, 2), ha='center', fontsize=plot_font_size, color='black', rotation=90)  # Adjust multiplier as needed
    #plt.annotate(percentage_text_lesser, xy=(i, min_value), ha='center', fontsize=plot_font_size, color='black')  # Adjust multiplier as needed
    plt.annotate(percentage_text_lesser, xy=(i, 0.01), ha='center', fontsize=plot_font_size, color='black')  # Adjust multiplier as needed

plt.annotate("% > 1", xy=(-1.6, 700000), ha='center', fontsize=plot_font_size + 6, fontweight='bold', color='black')  # Adjust multiplier as needed
plt.annotate("% = 1", xy=(-1.6, 2), ha='center', fontsize=plot_font_size + 6, fontweight='bold', color='black')  # Adjust multiplier as needed
plt.annotate("% < 1", xy=(-1.6, 0.009), ha='center', fontsize=plot_font_size + 6, fontweight='bold', color='black')  # Adjust multiplier as needed

plt.tight_layout()
plt.show()




################## 
# 3.3 Descriptive Statistics for CFM-ID Scores  
##################

# Plot distribution of CFM-ID scores - all features

df22['has_MS2'] = df22.groupby('Feature ID')['MS2 Score'].transform(lambda group: group.notna().any())
df22_has_MS2 = df22[df22['has_MS2'] == True]

# # Distribution plots of counts of structures
# plt.figure(figsize=(8, 6))
# sns.violinplot(y=df22_has_MS2['MS2 raw score'], color='skyblue')
# plt.title('MS2 scores for chemical candidates')
# plt.ylabel('MS2 similarity score')
# plt.grid(axis='x', linestyle='--', alpha=0.7)
# plt.show()

# # # Plot distribution of CFM-ID scores - correct vs. incorrect candidates
# # plt.figure(figsize=(10, 6))
# # sns.boxplot(x='source', y='MS2 raw score', data=ID_results_combined, palette='pastel')
# # plt.title("Comparison of MS2 scores for Level 1/SSM chemicals")
# # plt.xlabel('Candidate category')
# # plt.ylabel('MS2 score')
# # plt.grid(axis='y', linestyle='--', alpha=0.7)
# # #plt.yscale('log')
# # plt.show()

# # Plot point cloud using scatter plot
# jitter = np.random.uniform(-0.2, 0.2, size=len(df22_has_MS2))
# plt.figure(figsize=(10, 6))
# plt.scatter(df22_has_MS2['MS2 raw score'], jitter, alpha=0.2, color='blue', s=4)
# plt.title('MS2 raw scores for candidates with MS2 data')
# plt.xlabel('MS2 raw score')
# #plt.ylabel('Jitter')
# plt.grid(axis='x', linestyle='--', alpha=0.7)
# plt.yticks([])  # Hide y-axis ticks as jitter is just for visualization
# plt.show()

# Histograms of MS2 raw scores: All candidates
# Define bin size
bin_size = 0.2

# Calculate bins based on bin size
min_value = df22_has_MS2['MS2 raw score'].min()
max_value = df22_has_MS2['MS2 raw score'].max()
bins = np.arange(min_value, max_value + bin_size, bin_size)

# Plot histogram
plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(df22_has_MS2['MS2 raw score'], bins=bins, color='skyblue', edgecolor='black')

# Calculate percentages
percentages = (counts / counts.sum()) * 100

# Annotate percentages on the histogram
for count, percentage, patch in zip(counts, percentages, patches):
    # Calculate the center of the patch
    x_center = patch.get_x() + patch.get_width() / 2
    plt.annotate(f'{percentage:.1f}%', xy=(x_center, count),
                 xytext=(0, 5), textcoords='offset points', ha='center', va='bottom', fontsize=10)


# Adjust y-axis limit to provide extra space for annotations
plt.ylim(0, max(counts) * 1.1)

# Add labels and title
plt.title('MS$^2$ Raw Scores for Candidates with MS$^2$ Data: All Features')
plt.xlabel('MS$^2$ raw score')
plt.ylabel('Frequency')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show plot
plt.show()


# Histograms of MS2 raw scores: Candidates for Known Chemical Features
# Define bin size
bin_size = 0.2

ID_results_combined_has_MS2 = ID_results_combined[ID_results_combined['has_MS2'] == True]

debug = ID_results_combined_has_MS2[['Feature ID', 'MS2 Score', 'MS2 raw score']]

# Group by the 'Group' column
debug_grouped = ID_results_combined_has_MS2.groupby('Feature ID')

# Identify groups where all values in the 'Value' column are NaN
debug_nan_MS2 = [name for name, group in debug_grouped if group['MS2 raw score'].isna().all()]

print(debug_nan_MS2)

# Use same bins/ranges as all features
# min_value = ID_results_combined_has_MS2['MS2 raw score'].min()
# max_value = ID_results_combined_has_MS2['MS2 raw score'].max()
# bins = np.arange(min_value, max_value + bin_size, bin_size)

# Plot histogram
plt.figure(figsize=(10, 6))
counts, bins, patches = plt.hist(ID_results_combined_has_MS2['MS2 raw score'], bins=bins, color='skyblue', edgecolor='black')

# Calculate percentages
percentages = (counts / counts.sum()) * 100

# Annotate percentages on the histogram
for count, percentage, patch in zip(counts, percentages, patches):
    # Calculate the center of the patch
    x_center = patch.get_x() + patch.get_width() / 2
    plt.annotate(f'{percentage:.1f}%', xy=(x_center, count),
                 xytext=(0, 5), textcoords='offset points', ha='center', va='bottom', fontsize=10)


# Adjust y-axis limit to provide extra space for annotations
plt.ylim(0, max(counts) * 1.1)

# Add labels and title
#plt.title('MS$^2$ Raw Scores for Candidates with MS$^2$ Data: Known Chemical Features')
plt.title('MS$^2$ Raw Scores for Candidates with MS$^2$ Data: $\it{Known\ Chemical\ Features}$')
plt.xlabel('MS$^2$ raw score')
plt.ylabel('Frequency')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show plot
plt.show()




################## 
# 3.3 Descriptive Statistics for Hazard Scores
##################


# # Calculate average, minimum, and maximum for each column in the list
# ID_results_true_stats_hazard = ID_results_true[hazard_columns].agg(['mean', 'min', 'max'])
# ID_results_other_stats_hazard = ID_results_other[hazard_columns].agg(['mean', 'min', 'max'])


# Calculate median and IQR for each column - True candidate
results = {}
for column in hazard_columns:
    median = ID_results_true[column].median()
    q1 = ID_results_true[column].quantile(0.25)
    q3 = ID_results_true[column].quantile(0.75)
    iqr = q3 - q1
    mean = ID_results_true[column].mean()
    std_dev = ID_results_true[column].std()
    results[column] = {'Mean': mean, 'Standard Deviation': std_dev, 'Median': median, 'IQR': iqr}

# Convert results to DataFrame for better visualization
ID_results_true_stats_hazard_b = pd.DataFrame(results).T

# Calculate median and IQR for each column - Other candidates
results = {}
for column in hazard_columns:
    median = ID_results_other[column].median()
    q1 = ID_results_other[column].quantile(0.25)
    q3 = ID_results_other[column].quantile(0.75)
    iqr = q3 - q1
    mean = ID_results_other[column].mean()
    std_dev = ID_results_other[column].std()
    results[column] = {'Mean': mean, 'Standard Deviation': std_dev, 'Median': median, 'IQR': iqr}

# Convert results to DataFrame for better visualization
ID_results_other_stats_hazard_b = pd.DataFrame(results).T


# Calculate median and IQR for each column - All feature candidates
results = {}
for column in hazard_columns:
    median = df22[column].median()
    q1 = ID_results_true[column].quantile(0.25)
    q3 = ID_results_true[column].quantile(0.75)
    iqr = q3 - q1
    mean = ID_results_true[column].mean()
    std_dev = ID_results_true[column].std()
    results[column] = {'Mean': mean, 'Standard Deviation': std_dev, 'Median': median, 'IQR': iqr}

# Convert results to DataFrame for better visualization
All_results_true_stats_hazard_b = pd.DataFrame(results).T

# Function to round a float to 3 significant figures
def round_to_sig_figs(value, sig_figs=3):
    if value == 0:
        return 0
    else:
        return round(value, sig_figs - int(np.floor(np.log10(abs(value))) + 1))


ID_results_true_stats_hazard_b = ID_results_true_stats_hazard_b.fillna(0)
All_results_true_stats_hazard_b = All_results_true_stats_hazard_b.fillna(0)

# Create the third column with formatted strings
ID_results_true_stats_hazard_b['Mean (SD)'] = ID_results_true_stats_hazard_b.apply(lambda row: f"{round_to_sig_figs(row['Mean']):g} ({round_to_sig_figs(row['Standard Deviation']):g})", axis=1)
ID_results_true_stats_hazard_b['Median (IQR)'] = ID_results_true_stats_hazard_b.apply(lambda row: f"{round_to_sig_figs(row['Median']):g} ({round_to_sig_figs(row['IQR']):g})", axis=1)
ID_results_other_stats_hazard_b['Mean (SD)'] = ID_results_other_stats_hazard_b.apply(lambda row: f"{round_to_sig_figs(row['Mean']):g} ({round_to_sig_figs(row['Standard Deviation']):g})", axis=1)
ID_results_other_stats_hazard_b['Median (IQR)'] = ID_results_other_stats_hazard_b.apply(lambda row: f"{round_to_sig_figs(row['Median']):g} ({round_to_sig_figs(row['IQR']):g})", axis=1)
All_results_true_stats_hazard_b['Mean (SD)'] = All_results_true_stats_hazard_b.apply(lambda row: f"{round_to_sig_figs(row['Mean']):g} ({round_to_sig_figs(row['Standard Deviation']):g})", axis=1)
All_results_true_stats_hazard_b['Median (IQR)'] = All_results_true_stats_hazard_b.apply(lambda row: f"{round_to_sig_figs(row['Median']):g} ({round_to_sig_figs(row['IQR']):g})", axis=1)


# Create scatterplot of QAH vs. Completeness Score - all features
# Create scatter plot using seaborn
plt.figure(figsize=(8, 6))
sns.scatterplot(x='Hazard Score', y='Hazard Completeness Score', data=df22, color='blue', marker='o', s=100, alpha=0.1)

# Customize the plot
plt.title('Hazard Data for All Feature Candidates')
plt.xlabel('Hazard Score')
plt.ylabel('Hazard Completeness Score')
plt.grid(True)
plt.show()


# Create scatterplot of QAH vs. Completeness Score - ID features
# Define colors for True and False
colors = {True: 'blue', False: 'orange'}

# Create scatter plot with different colors based on the 'is_special' column
plt.figure(figsize=(8, 6))
for is_special, color in colors.items():
    subset = ID_results_combined[ID_results_combined['SSM or level 1 structure'] == is_special]
    plt.scatter(subset['Hazard Score'], subset['Hazard Completeness Score'], color=color, label=str(is_special), marker='o', s=50, alpha=0.4)

# Customize the plot
plt.title('Hazard Data for Level 1/SSM Candidates')
plt.xlabel('Hazard Score')
plt.ylabel('Hazard Completeness Score')
plt.legend()
plt.grid(True)

plt.show()

# Create scatterplot of QAH vs. Completeness Score - All features - Correct candidates highlighted
# Define colors for True and False
colors = {True: 'blue', False: 'orange'}
alphas = {True: 0.9, False: 0.1}
markers = {True: '^', False: 'o'}

# Create scatter plot with different colors based on the 'is_special' column
plt.figure(figsize=(8, 6))

# Plot False markers first
subset_false = df22[df22['SSM or level 1 structure'] == False]
plt.scatter(subset_false['Hazard Score'], subset_false['Hazard Completeness Score'], color=colors[False], alpha=alphas[False],
            label='False', marker=markers[False], s=50)

# Plot True markers second
subset_true = df22[df22['SSM or level 1 structure'] == True]
plt.scatter(subset_true['Hazard Score'], subset_true['Hazard Completeness Score'], color=colors[True], alpha=alphas[True],
            label='True', marker=markers[True], s=50)

# Customize the plot
plt.title('Hazard Data for All Features - Correct Candidates Highlighted')
plt.xlabel('Hazard Score')
plt.ylabel('Hazard Completeness Score')
plt.legend()
#plt.grid(True)

plt.show()



# Create scatterplot of QAH vs. Completeness Score - 3 categories: correct, incorrect and other
# Need to run the section of code above for some reason first, or else datapoints don't show up
# Define colors for True and False
colors = {True: 'blue', False: 'orange'}
alphas = {True: 0.9, False: 0.2}
markers = {True: '^', False: 'o'}
marker_size = 50
# Create scatter plot with different colors based on the 'is_special' column
plt.figure(figsize=(8, 6))

# Plot other markers first
subset_other = df22[df22['SSM or level 1 feature'] == False]
plt.scatter(subset_false['Hazard Score'], subset_false['Hazard Completeness Score'], color='grey', alpha=0.1,
            label='Unknown', marker=markers[False], s=marker_size)

# Plot False markers first
subset_false = df22[(df22['SSM or level 1 structure'] == False) & (df22['SSM or level 1 feature'] == True)]
#subset_false = df22[df22['SSM or level 1 structure'] == False]
plt.scatter(subset_false['Hazard Score'], subset_false['Hazard Completeness Score'], color='#ff7f0e', alpha=0.4,
            label='Known - Incorrect Candidate', marker=markers[False], s=marker_size)

# Plot True markers second
subset_true = df22[(df22['SSM or level 1 structure'] == True) & (df22['SSM or level 1 feature'] == True)]
#subset_true = df22[df22['SSM or level 1 structure'] == True]
plt.scatter(subset_true['Hazard Score'], subset_true['Hazard Completeness Score'], color='blue', alpha=0.9,
            label='Known - Correct Candidate', marker=markers[True], s=marker_size)

# Customize the plot
plt.title('Distribution of Hazard Data for All Features')
plt.xlabel('Hazard Score')
plt.ylabel('Hazard Completeness Score')
plt.legend()
plt.show()


################## 
# 3.4 Aggregate Scoring and Visualization  
##################

# Plot A:
cmp='YlOrRd'
font_size = 20
occurrence_scaler = 30
alpha_value = 0.9

df_tri_true = df22[df22['SSM or level 1 structure'] == True]

# Filter by has MS2, and fill nan values with zero
df_tri_true = df_tri_true[df_tri_true['has_MS2'] == True]
df_tri_true['MS2 Score'] = df_tri_true['MS2 Score'].fillna(0)

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df22['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_true = pd.concat([df_tri_true, dummy_df], ignore_index=True)
df_tri_true = pd.concat([df_tri_true, dummy_df_b], ignore_index=True)

debug = df_tri_true[['Feature ID', 'DTXCID', 'Metadata Score', 'MS2 Score', 'Final Occurrence Count', 'Hazard Score']]

df_tri_true = df_tri_true.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_true['Metadata Score'], df_tri_true['MS2 Score'], alpha = alpha_value, s=df_tri_true['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_true['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

#debug = df_tri_true[['Feature ID', 'Metadata Score', 'MS2 Score', 'Final Occurrence Count', 'Hazard Score']]

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel(r"MS$^2$ score", fontsize=font_size)
plt.title('A) Correct Candidates of Known Features', fontsize=font_size)

# Calculate points in each quadrant - Start Q1 and Q3 at minus due to having dummy points
quadrants = {'Q1': -1, 'Q2': 0, 'Q3': -1, 'Q4': 0}

for index, row in df_tri_true.iterrows():
    x = row['Metadata Score']
    y = row['MS2 Score']
    
    if x >= x_midpoint and y >= y_midpoint:
        quadrants['Q1'] += 1
    elif x < x_midpoint and y >= y_midpoint:
        quadrants['Q2'] += 1
    elif x < x_midpoint and y < y_midpoint:
        quadrants['Q3'] += 1
    elif x >= x_midpoint and y < y_midpoint:
        quadrants['Q4'] += 1

# Calculate percentages
# First drop rows that don't have MS2 scores (as these are not plotted)
df_tri_true_has_MS2 = df_tri_true.dropna(subset=['MS2 Score'])
total_points = len(df_tri_true_has_MS2) - 2 # Subtract two dummy points
quadrant_percentages = {key: (value / total_points) * 100 for key, value in quadrants.items()}

# Annotate each quadrant with count and percentage
plt.text(x_midpoint + 0.2, y_midpoint + 0.02, f"$\\bf{{Q1}}$: {quadrants['Q1']} ({quadrant_percentages['Q1']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint - 4.3, y_midpoint + 0.02, f"$\\bf{{Q2}}$: {quadrants['Q2']} ({quadrant_percentages['Q2']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(-0.3, -0.03, f"$\\bf{{Q3}}$: {quadrants['Q3']} ({quadrant_percentages['Q3']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint + 0.2, -0.03, f"$\\bf{{Q4}}$: {quadrants['Q4']} ({quadrant_percentages['Q4']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))




# Plot 8 - True chemicals, other structures
df_tri_false = df22[(df22['SSM or level 1 structure'] == False) & (df22['SSM or level 1 feature'] == True)]

# Filter by has MS2, and fill nan values with zero
df_tri_false = df_tri_false[df_tri_false['has_MS2'] == True]
df_tri_false['MS2 Score'] = df_tri_false['MS2 Score'].fillna(0)


dummy_x = np.nan
dummy_y = np.nan
dummy_color = df_tri_false['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_false = pd.concat([df_tri_false, dummy_df], ignore_index=True)
df_tri_false = pd.concat([df_tri_false, dummy_df_b], ignore_index=True)

df_tri_false = df_tri_false.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_false['Metadata Score'], df_tri_false['MS2 Score'], alpha = alpha_value, s=df_tri_false['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_false['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel(r"MS$^2$ score", fontsize=font_size)
plt.title('B) Incorrect Candidates of Known Features', fontsize=font_size)

# Calculate points in each quadrant - Start Q1 and Q3 at minus due to having dummy points
quadrants = {'Q1': -1, 'Q2': 0, 'Q3': -1, 'Q4': 0}

for index, row in df_tri_false.iterrows():
    x = row['Metadata Score']
    y = row['MS2 Score']
    
    if x >= x_midpoint and y >= y_midpoint:
        quadrants['Q1'] += 1
    elif x < x_midpoint and y >= y_midpoint:
        quadrants['Q2'] += 1
    elif x < x_midpoint and y < y_midpoint:
        quadrants['Q3'] += 1
    elif x >= x_midpoint and y < y_midpoint:
        quadrants['Q4'] += 1

# Calculate percentages
df_tri_false_has_MS2 = df_tri_false.dropna(subset=['MS2 Score'])
total_points = len(df_tri_false_has_MS2) - 2 # Subtract two dummy points
quadrant_percentages = {key: (value / total_points) * 100 for key, value in quadrants.items()}

# Annotate each quadrant with count and percentage
plt.text(x_midpoint + 0.2, y_midpoint + 0.02, f"$\\bf{{Q1}}$: {quadrants['Q1']} ({quadrant_percentages['Q1']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint - 4.3, y_midpoint + 0.02, f"$\\bf{{Q2}}$: {quadrants['Q2']} ({quadrant_percentages['Q2']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(-0.3, -0.03, f"$\\bf{{Q3}}$: {quadrants['Q3']} ({quadrant_percentages['Q3']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint + 0.2, -0.03, f"$\\bf{{Q4}}$: {quadrants['Q4']} ({quadrant_percentages['Q4']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))


# Plot C - Other chemical features
df_tri_other = df22[df22['SSM or level 1 feature'] == False]

# Filter by has MS2, and fill nan values with zero
df_tri_other = df_tri_other[df_tri_other['has_MS2'] == True]
df_tri_other['MS2 Score'] = df_tri_other['MS2 Score'].fillna(0)


dummy_x = np.nan
dummy_y = np.nan
dummy_color = df_tri_other['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_other = pd.concat([df_tri_other, dummy_df], ignore_index=True)
df_tri_other = pd.concat([df_tri_other, dummy_df_b], ignore_index=True)

df_tri_other = df_tri_other.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_other['Metadata Score'], df_tri_other['MS2 Score'], alpha = alpha_value, s=df_tri_other['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_other['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel(r"MS$^2$ score", fontsize=font_size)
plt.title('C) All Candidates of Unknown Features', fontsize=font_size)

# Calculate points in each quadrant - Start Q1 and Q3 at minus due to having dummy points
quadrants = {'Q1': -1, 'Q2': 0, 'Q3': -1, 'Q4': 0}

for index, row in df_tri_other.iterrows():
    x = row['Metadata Score']
    y = row['MS2 Score']
    
    if x >= x_midpoint and y >= y_midpoint:
        quadrants['Q1'] += 1
    elif x < x_midpoint and y >= y_midpoint:
        quadrants['Q2'] += 1
    elif x < x_midpoint and y < y_midpoint:
        quadrants['Q3'] += 1
    elif x >= x_midpoint and y < y_midpoint:
        quadrants['Q4'] += 1

# Calculate percentages
df_tri_other_has_MS2 = df_tri_other.dropna(subset=['MS2 Score'])
total_points = len(df_tri_other_has_MS2) - 2 # Subtract two dummy points
quadrant_percentages = {key: (value / total_points) * 100 for key, value in quadrants.items()}

# Annotate each quadrant with count and percentage
plt.text(x_midpoint + 0.2, y_midpoint + 0.02, f"$\\bf{{Q1}}$: {quadrants['Q1']} ({quadrant_percentages['Q1']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint - 4.3, y_midpoint + 0.02, f"$\\bf{{Q2}}$: {quadrants['Q2']} ({quadrant_percentages['Q2']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(-0.3, -0.03, f"$\\bf{{Q3}}$: {quadrants['Q3']} ({quadrant_percentages['Q3']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint + 0.2, -0.03, f"$\\bf{{Q4}}$: {quadrants['Q4']} ({quadrant_percentages['Q4']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))





# Do some filtering on the plots first to only show those with high occurrence and/or high hazard

# Plot A:
cmp='YlOrRd'
font_size = 20
occurrence_scaler = 30
alpha_value = 0.9
# Plot A - True chemicals, true structures
df_tri_true = df22[df22['SSM or level 1 structure'] == True]
df_tri_true = df_tri_true[(df_tri_true['Hazard Score'] > 4) & (df_tri_true['Final Occurrence Count'] >= 6)]

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df22['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_true = pd.concat([df_tri_true, dummy_df], ignore_index=True)
df_tri_true = pd.concat([df_tri_true, dummy_df_b], ignore_index=True)

df_tri_true = df_tri_true.sort_values(by='Hazard Score', ascending=True)


plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_true['Metadata Score'], df_tri_true['MS2 Score'], alpha = alpha_value, s=df_tri_true['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_true['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

#debug = df_tri_true[['Feature ID', 'Metadata Score', 'MS2 Score', 'Final Occurrence Count', 'Hazard Score']]

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('A) Correct Candidates of Chemical Features', fontsize=font_size)

# Calculate points in each quadrant - Start Q1 and Q3 at minus due to having dummy points
quadrants = {'Q1': -1, 'Q2': 0, 'Q3': -1, 'Q4': 0}

for index, row in df_tri_true.iterrows():
    x = row['Metadata Score']
    y = row['MS2 Score']
    
    if x >= x_midpoint and y >= y_midpoint:
        quadrants['Q1'] += 1
    elif x < x_midpoint and y >= y_midpoint:
        quadrants['Q2'] += 1
    elif x < x_midpoint and y < y_midpoint:
        quadrants['Q3'] += 1
    elif x >= x_midpoint and y < y_midpoint:
        quadrants['Q4'] += 1

# Calculate percentages
# First drop rows that don't have MS2 scores (as these are not plotted)
df_tri_true_has_MS2 = df_tri_true.dropna(subset=['MS2 Score'])
total_points = len(df_tri_true_has_MS2) - 2 # Subtract two dummy points
quadrant_percentages = {key: (value / total_points) * 100 for key, value in quadrants.items()}

# Annotate each quadrant with count and percentage
plt.text(x_midpoint + 0.2, y_midpoint + 0.02, f"$\\bf{{Q1}}$: {quadrants['Q1']} ({quadrant_percentages['Q1']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint - 4.3, y_midpoint + 0.02, f"$\\bf{{Q2}}$: {quadrants['Q2']} ({quadrant_percentages['Q2']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(-0.3, -0.03, f"$\\bf{{Q3}}$: {quadrants['Q3']} ({quadrant_percentages['Q3']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint + 0.2, -0.03, f"$\\bf{{Q4}}$: {quadrants['Q4']} ({quadrant_percentages['Q4']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))



# Plot 8 - True chemicals, other structures
df_tri_false = df22[(df22['SSM or level 1 structure'] == False) & (df22['SSM or level 1 feature'] == True)]
df_tri_false = df_tri_false[(df_tri_false['Hazard Score'] >= 4) & (df_tri_false['Final Occurrence Count'] > 6)]

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df_tri_false['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_false = pd.concat([df_tri_false, dummy_df], ignore_index=True)
df_tri_false = pd.concat([df_tri_false, dummy_df_b], ignore_index=True)

df_tri_false = df_tri_false.sort_values(by='Hazard Score', ascending=True)


debug = df_tri_false[['Metadata Score', 'MS2 Score', 'Hazard Score', 'Final Occurrence Count']]

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_false['Metadata Score'], df_tri_false['MS2 Score'], alpha = alpha_value, s=df_tri_false['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_false['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('B) Incorrect Candidates of Known Features', fontsize=font_size)

# Calculate points in each quadrant - Start Q1 and Q3 at minus due to having dummy points
quadrants = {'Q1': -1, 'Q2': 0, 'Q3': -1, 'Q4': 0}

for index, row in df_tri_false.iterrows():
    x = row['Metadata Score']
    y = row['MS2 Score']
    
    if x >= x_midpoint and y >= y_midpoint:
        quadrants['Q1'] += 1
    elif x < x_midpoint and y >= y_midpoint:
        quadrants['Q2'] += 1
    elif x < x_midpoint and y < y_midpoint:
        quadrants['Q3'] += 1
    elif x >= x_midpoint and y < y_midpoint:
        quadrants['Q4'] += 1

# Calculate percentages
df_tri_false_has_MS2 = df_tri_false.dropna(subset=['MS2 Score'])
total_points = len(df_tri_false_has_MS2) - 2 # Subtract two dummy points
quadrant_percentages = {key: (value / total_points) * 100 for key, value in quadrants.items()}

# Annotate each quadrant with count and percentage
plt.text(x_midpoint + 0.2, y_midpoint + 0.02, f"$\\bf{{Q1}}$: {quadrants['Q1']} ({quadrant_percentages['Q1']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint - 4.3, y_midpoint + 0.02, f"$\\bf{{Q2}}$: {quadrants['Q2']} ({quadrant_percentages['Q2']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(-0.3, -0.03, f"$\\bf{{Q3}}$: {quadrants['Q3']} ({quadrant_percentages['Q3']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint + 0.2, -0.03, f"$\\bf{{Q4}}$: {quadrants['Q4']} ({quadrant_percentages['Q4']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))


# Plot C/D - Other chemical features
df_tri_other = df22[df22['SSM or level 1 feature'] == False]
df_tri_other = df_tri_other[(df_tri_other['Hazard Score'] >= 4) & (df_tri_other['Final Occurrence Count'] > 6)]

# Filter by has MS2, and fill nan values with zero
df_tri_other = df_tri_other[df_tri_other['has_MS2'] == True]
df_tri_other['MS2 Score'] = df_tri_other['MS2 Score'].fillna(0)

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df_tri_other['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_other = pd.concat([df_tri_other, dummy_df], ignore_index=True)
df_tri_other = pd.concat([df_tri_other, dummy_df_b], ignore_index=True)

df_tri_other = df_tri_other.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_other['Metadata Score'], df_tri_other['MS2 Score'], alpha = alpha_value, s=df_tri_other['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_other['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel(r"MS$^2$ score", fontsize=font_size)
plt.title('D) Filtered Candidates of Unknown Features', fontsize=font_size)

# Calculate points in each quadrant - Start Q1 and Q3 at minus due to having dummy points
quadrants = {'Q1': -1, 'Q2': 0, 'Q3': -1, 'Q4': 0}

for index, row in df_tri_other.iterrows():
    x = row['Metadata Score']
    y = row['MS2 Score']
    
    if x >= x_midpoint and y >= y_midpoint:
        quadrants['Q1'] += 1
    elif x < x_midpoint and y >= y_midpoint:
        quadrants['Q2'] += 1
    elif x < x_midpoint and y < y_midpoint:
        quadrants['Q3'] += 1
    elif x >= x_midpoint and y < y_midpoint:
        quadrants['Q4'] += 1

# Calculate percentages
df_tri_other_has_MS2 = df_tri_other.dropna(subset=['MS2 Score'])
total_points = len(df_tri_other_has_MS2) - 2 # Subtract two dummy points
quadrant_percentages = {key: (value / total_points) * 100 for key, value in quadrants.items()}

# Annotate each quadrant with count and percentage
plt.text(x_midpoint + 0.2, y_midpoint + 0.02, f"$\\bf{{Q1}}$: {quadrants['Q1']} ({quadrant_percentages['Q1']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint - 4.3, y_midpoint + 0.02, f"$\\bf{{Q2}}$: {quadrants['Q2']} ({quadrant_percentages['Q2']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(-0.3, -0.03, f"$\\bf{{Q3}}$: {quadrants['Q3']} ({quadrant_percentages['Q3']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint + 0.2, -0.03, f"$\\bf{{Q4}}$: {quadrants['Q4']} ({quadrant_percentages['Q4']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))





# Plot just one Known feature
#df_tri_one_feature = df22[df22['Feature ID'] == 247]
df_tri_one_feature = df22[df22['Feature ID'] == 247]

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df22['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_one_feature = pd.concat([df_tri_one_feature, dummy_df], ignore_index=True)
df_tri_one_feature = pd.concat([df_tri_one_feature, dummy_df_b], ignore_index=True)

df_tri_one_feature = df_tri_one_feature.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
# scatter = plt.scatter(df_tri_one_feature['Metadata Score'], df_tri_one_feature['MS2 Score'], alpha = 0.9, s=150, 
#                       c=df_tri_one_feature['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
scatter = plt.scatter(df_tri_one_feature['Metadata Score'], df_tri_one_feature['MS2 Score'], alpha = alpha_value, s=df_tri_one_feature['Final Occurrence Count'] * occurrence_scaler*2, 
                      c=df_tri_one_feature['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

#debug = df_tri_true[['Feature ID', 'Metadata Score', 'MS2 Score', 'Final Occurrence Count', 'Hazard Score']]

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
#plt.title('A) Correct Candidates of Known Features', fontsize=font_size)



# Calculate points in each quadrant - Start Q1 and Q3 at minus due to having dummy points
quadrants = {'Q1': -1, 'Q2': 0, 'Q3': -1, 'Q4': 0}

for index, row in df_tri_one_feature.iterrows():
    x = row['Metadata Score']
    y = row['MS2 Score']
    
    if x >= x_midpoint and y >= y_midpoint:
        quadrants['Q1'] += 1
    elif x < x_midpoint and y >= y_midpoint:
        quadrants['Q2'] += 1
    elif x < x_midpoint and y < y_midpoint:
        quadrants['Q3'] += 1
    elif x >= x_midpoint and y < y_midpoint:
        quadrants['Q4'] += 1

# Calculate percentages
# First drop rows that don't have MS2 scores (as these are not plotted)
df_tri_one_feature_has_MS2 = df_tri_one_feature.dropna(subset=['MS2 Score'])
total_points = len(df_tri_one_feature_has_MS2) - 2 # Subtract two dummy points
quadrant_percentages = {key: (value / total_points) * 100 for key, value in quadrants.items()}

# Annotate each quadrant with count and percentage
plt.text(x_midpoint + 0.2, y_midpoint + 0.02, f"$\\bf{{Q1}}$: {quadrants['Q1']} ({quadrant_percentages['Q1']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint - 4.3, y_midpoint + 0.02, f"$\\bf{{Q2}}$: {quadrants['Q2']} ({quadrant_percentages['Q2']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(-0.3, -0.03, f"$\\bf{{Q3}}$: {quadrants['Q3']} ({quadrant_percentages['Q3']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))
plt.text(x_midpoint + 0.2, -0.03, f"$\\bf{{Q4}}$: {quadrants['Q4']} ({quadrant_percentages['Q4']:.1f}%)", fontsize=font_size-4, ha='left', va='bottom',
         bbox=dict(facecolor='white', alpha=1, edgecolor='black'))






# Plot a marker size legend
df_tri_other = df22[df22['SSM or level 1 feature'] == False]
df_tri_other = df_tri_other[(df_tri_other['Hazard Score'] >= 4) & (df_tri_other['Final Occurrence Count'] > 6)]

# Filter by has MS2, and fill nan values with zero
df_tri_other = df_tri_other[df_tri_other['has_MS2'] == True]
df_tri_other['MS2 Score'] = df_tri_other['MS2 Score'].fillna(0)

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df_tri_other['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score':[8], 'MS2 Score':[1], 'Final Occurrence Count':0.0000001, 'Hazard Score':[12]})
dummy_df_b = pd.DataFrame({'Metadata Score':[0], 'MS2 Score':[0], 'Final Occurrence Count':0.0000001, 'Hazard Score':[0]})
df_tri_other = pd.concat([df_tri_other, dummy_df], ignore_index=True)
df_tri_other = pd.concat([df_tri_other, dummy_df_b], ignore_index=True)

df_tri_other = df_tri_other.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_other['Metadata Score'], df_tri_other['MS2 Score'], alpha = alpha_value, s=df_tri_other['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_other['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('D) Filtered Candidates of Unknown Features', fontsize=font_size)


# Create the markersize and hazard legend
# Convert the colormap name to a colormap object
cmp_b = cm.get_cmap('YlOrRd')

# Your existing data preparation code
df_tri_other = df22[df22['SSM or level 1 feature'] == False]
df_tri_other = df_tri_other[(df_tri_other['Hazard Score'] >= 4) & (df_tri_other['Final Occurrence Count'] > 6)]
df_tri_other = df_tri_other[df_tri_other['has_MS2'] == True]
df_tri_other['MS2 Score'] = df_tri_other['MS2 Score'].fillna(0)

dummy_df = pd.DataFrame({'Metadata Score': [8], 'MS2 Score': [1], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [12]})
dummy_df_b = pd.DataFrame({'Metadata Score': [0], 'MS2 Score': [0], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [0]})
df_tri_other = pd.concat([df_tri_other, dummy_df], ignore_index=True)
df_tri_other = pd.concat([df_tri_other, dummy_df_b], ignore_index=True)

df_tri_other = df_tri_other.sort_values(by='Hazard Score', ascending=True)

# Scatter plot code
plt.figure(figsize=(10, 8))
plt.xticks(fontsize=font_size-4)
plt.yticks(fontsize=font_size-4)
scatter = plt.scatter(df_tri_other['Metadata Score'], df_tri_other['MS2 Score'], alpha=alpha_value, s=df_tri_other['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_other['Hazard Score'], cmap=cmp_b, edgecolors='black', linewidths=0.5)
cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size) 
cbar.ax.tick_params(labelsize=font_size-4)
plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

# Calculate midpoints for x and y axes
x_midpoint = (plt.xlim()[0] + plt.xlim()[1]) / 2
y_midpoint = (plt.ylim()[0] + plt.ylim()[1]) / 2

# Draw vertical and horizontal lines at the midpoints
plt.axvline(x=x_midpoint, color='black', linestyle='--', linewidth=3)
plt.axhline(y=y_midpoint, color='black', linestyle='--', linewidth=3)

plt.xlabel("Metadata score", fontsize=font_size)
plt.ylabel("MS2 score", fontsize=font_size)
plt.title('D) Filtered Candidates of Unknown Features', fontsize=font_size)

# Add marker size legend outside the plot with reversed order
sizes = np.arange(1, 14)  # Sizes from 1 to 13
colors = [cmp_b(i / 13) for i in sizes]  # Gradient colors using the colormap

# Reverse the order of legend elements and modify the labels
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', markersize=np.sqrt(size * occurrence_scaler), 
                              label=f"{round(size * 7.692)}%", markeredgewidth=0.5, markeredgecolor='black') for size, color in zip(sizes[::-1], colors[::-1])]


# Position the legend on the right side, lower down to avoid overlap with the colorbar
plt.legend(handles=legend_elements, title='Occurrence %', loc='upper left', bbox_to_anchor=(1.3, 0.8), fontsize=font_size-4, title_fontsize=font_size-4, frameon=False)

plt.show()




# Figures for Heather/Graphical Abstract
# Assuming df22 is already defined
df_tri_other = df22[df22['SSM or level 1 feature'] == True]

# Filter by has MS2, and fill nan values with zero
df_tri_other = df_tri_other[df_tri_other['has_MS2'] == True]
df_tri_other['MS2 Score'] = df_tri_other['MS2 Score'].fillna(0)

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df_tri_other['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score': [8], 'MS2 Score': [1], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [12]})
dummy_df_b = pd.DataFrame({'Metadata Score': [0], 'MS2 Score': [0], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [0]})
df_tri_other = pd.concat([df_tri_other, dummy_df], ignore_index=True)
df_tri_other = pd.concat([df_tri_other, dummy_df_b], ignore_index=True)

df_tri_other = df_tri_other.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
scatter = plt.scatter(df_tri_other['Metadata Score'], df_tri_other['MS2 Score'], alpha=alpha_value, 
                      s=df_tri_other['Final Occurrence Count'] * occurrence_scaler, 
                      c=df_tri_other['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)

cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size + 8, labelpad=15)  # Adjust labelpad for spacing


# Remove tick marks and labels from axes
plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, 
                labelbottom=False, labelleft=False)

# Remove tick marks and labels from colorbar
cbar.set_ticks([])
cbar.ax.set_yticklabels([])

plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

plt.xlabel("Metadata Score", fontsize=font_size + 8)
plt.ylabel(r"MS$^2$ Match Score", fontsize=font_size + 8)
# plt.title('C) All Candidates of Unknown Features', fontsize=font_size)

plt.show()




# Figures for Heather/Graphical Abstract - Just amphetamine feature
# Assuming df22 is already defined
df_tri_other = df22[df22['SSM or level 1 feature'] == True]
df_tri_other = df_tri_other[df_tri_other['Feature ID'] == 247]

# Filter by has MS2, and fill nan values with zero
df_tri_other = df_tri_other[df_tri_other['has_MS2'] == True]
df_tri_other['MS2 Score'] = df_tri_other['MS2 Score'].fillna(0)

# dummy_x = np.nan
# dummy_y = np.nan
# dummy_color = df_tri_other['Hazard Score'].max()
# dummy_df = pd.DataFrame({'Metadata Score': [8], 'MS2 Score': [1], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [12]})
# dummy_df_b = pd.DataFrame({'Metadata Score': [0], 'MS2 Score': [0], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [0]})
# df_tri_other = pd.concat([df_tri_other, dummy_df], ignore_index=True)
# df_tri_other = pd.concat([df_tri_other, dummy_df_b], ignore_index=True)

df_tri_other = df_tri_other.sort_values(by='Hazard Score', ascending=True)

plt.figure(figsize=(10, 8))
scatter = plt.scatter(df_tri_other['Metadata Score'], df_tri_other['MS2 Score'], alpha=alpha_value, 
                      s=500, 
                      c=df_tri_other['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=2)

cbar = plt.colorbar(scatter)
cbar.set_label('Hazard Score', size=font_size + 8, labelpad=15)  # Adjust labelpad for spacing


# Remove tick marks and labels from axes
plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, 
                labelbottom=False, labelleft=False)

# Remove tick marks and labels from colorbar
cbar.set_ticks([])
cbar.ax.set_yticklabels([])

plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

plt.xlabel("Metadata Score", fontsize=font_size + 8)
plt.ylabel(r"MS$^2$ Match Score", fontsize=font_size + 8)
# plt.title('C) All Candidates of Unknown Features', fontsize=font_size)

plt.show()





# Assuming df22 is already defined
df_tri_other = df22[df22['SSM or level 1 feature'] == True]

# Filter by has MS2, and fill nan values with zero
df_tri_other = df_tri_other[df_tri_other['has_MS2'] == True]
df_tri_other['MS2 Score'] = df_tri_other['MS2 Score'].fillna(0)

dummy_x = np.nan
dummy_y = np.nan
dummy_color = df_tri_other['Hazard Score'].max()
dummy_df = pd.DataFrame({'Metadata Score': [8], 'MS2 Score': [1], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [12]})
dummy_df_b = pd.DataFrame({'Metadata Score': [0], 'MS2 Score': [0], 'Final Occurrence Count': 0.0000001, 'Hazard Score': [0]})
df_tri_other = pd.concat([df_tri_other, dummy_df], ignore_index=True)
df_tri_other = pd.concat([df_tri_other, dummy_df_b], ignore_index=True)

df_tri_other = df_tri_other.sort_values(by='Hazard Score', ascending=True)

# Separate the rows with 'DTXCID402600'
highlight_df = df_tri_other[df_tri_other['DTXCID'] == 'DTXCID402600']
normal_df = df_tri_other[df_tri_other['DTXCID'] != 'DTXCID402600']

plt.figure(figsize=(10, 8))

# Plot normal markers
plt.scatter(normal_df['Metadata Score'], normal_df['MS2 Score'], alpha=alpha_value, 
            s=normal_df['Final Occurrence Count'] * occurrence_scaler, 
            c=normal_df['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5)

# Plot highlighted markers with extra bold outline
plt.scatter(highlight_df['Metadata Score'], highlight_df['MS2 Score'], alpha=alpha_value, 
            s=highlight_df['Final Occurrence Count'] * occurrence_scaler, 
            c=highlight_df['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=30, label='Highlighted')

# Plot highlighted markers with extra bold outline
plt.scatter(highlight_df['Metadata Score'], highlight_df['MS2 Score'], alpha=alpha_value, 
            s=highlight_df['Final Occurrence Count'] * occurrence_scaler, 
            c=highlight_df['Hazard Score'], cmap=cmp, edgecolors='black', linewidths=0.5, label='Highlighted')

cbar = plt.colorbar()
cbar.set_label('Hazard / Exposure Score', size=font_size + 8, labelpad=15)

# Remove tick marks and labels from axes
plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, 
                labelbottom=False, labelleft=False)

# Remove tick marks and labels from colorbar
cbar.set_ticks([])
cbar.ax.set_yticklabels([])

plt.xlim(-0.5, 8.5)
plt.ylim(-0.05, 1.05)

plt.xlabel("Metadata Score", fontsize=font_size + 8)
plt.ylabel(r"MS$^2$ Match Score", fontsize=font_size + 8)
# plt.title('C) All Candidates of Unknown Features', fontsize=font_size)

plt.show()



################## 
# 3.4 Pair-wise Spearman Correlations
##################

# all_columns_c = ['Patents', 'Articles', 'PubMed Records', 'Sources', 'AMOS spectra', 'AMOS methods',  
#                     'AMOS fact sheets',  'Water lists', 'MS2 Score', 'Hazard Score']

all_columns_c = ['Patents', 'Articles', 'PubMed Records', 'Sources', 'AMOS spectra', 'AMOS methods',  
                    'AMOS fact sheets',  'Water lists', 'Metadata Score', 'MS2 Score', 'Hazard Score']

# Initialize an empty DataFrame to store correlation results
correlation_matrix = pd.DataFrame(index=all_columns_c, columns=all_columns_c)

# df23 = df22.copy()
# df23['MS2 Score'] = df23['MS2 Score'].fillna(0)

# df24 = df23[all_columns_c]
# df24 = df24.fillna(0)

# # Calculate Spearman's correlation matrix
# correlation_matrix = df24.corr(method='spearman')

# # Create a mask for the diagonal
# mask = np.eye(len(correlation_matrix), dtype=bool)

# # Plot the correlation matrix as a heatmap
# plt.figure(figsize=(8, 6))

# # Plot heatmap without annotating the diagonal
# sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap='Blues', cbar=True, linewidths=.5)

# # Manually remove annotations from the diagonal
# for i in range(len(correlation_matrix)):
#     plt.gca().texts[i * (len(correlation_matrix) + 1)].set_text('')

# plt.title("Spearman's Correlation Matrix for Metadata, MS2 and Hazard Data Fields")
# plt.show()



# Assuming df22 and all_columns_c are already defined
df23 = df22.copy()
df23['MS2 Score'] = df23['MS2 Score'].fillna(0)

df24 = df23[all_columns_c]
df24 = df24.fillna(0)

# Calculate Spearman's correlation matrix
correlation_matrix = df24.corr(method='spearman')

# Create a mask for the diagonal
mask = np.eye(len(correlation_matrix), dtype=bool)

# Plot the correlation matrix as a heatmap
plt.figure(figsize=(8, 6))

# Plot heatmap without annotating the diagonal
sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap='Blues', cbar=True, linewidths=.5)

# Manually remove annotations from the diagonal
for i in range(len(correlation_matrix)):
    plt.gca().texts[i * (len(correlation_matrix) + 1)].set_text('')

# Get the current axis
ax = plt.gca()

# Access x-axis tick labels and modify "MS2 Score" to have a bold superscript "2"
current_labels = [label.get_text() for label in ax.get_xticklabels()]
new_labels = [label if label != "MS2 Score" else r"MS$\mathbf{^2}$ Score" for label in current_labels]
new_labels = [label if label != "Patents" else "PubChem Patents" for label in new_labels]
new_labels = [label if label != "Articles" else "PubChem Articles" for label in new_labels]
new_labels = [label if label != "PubMed Records" else "PubMed Articles" for label in new_labels]
new_labels = [label if label != "Sources" else "PubChem Sources" for label in new_labels]
new_labels = [label if label != "Water lists" else "Dashboard Water Lists" for label in new_labels]
ax.set_xticklabels(new_labels)

current_labels = [label.get_text() for label in ax.get_yticklabels()]
new_labels = [label if label != "MS2 Score" else r"MS$\mathbf{^2}$ Score" for label in current_labels]
new_labels = [label if label != "Patents" else "PubChem Patents" for label in new_labels]
new_labels = [label if label != "Articles" else "PubChem Articles" for label in new_labels]
new_labels = [label if label != "PubMed Records" else "PubMed Articles" for label in new_labels]
new_labels = [label if label != "Sources" else "PubChem Sources" for label in new_labels]
new_labels = [label if label != "Water lists" else "Dashboard Water Lists" for label in new_labels]
ax.set_yticklabels(new_labels)

# Make the last three x-axis labels bold
for label in ax.get_xticklabels()[-3:]:
    label.set_fontweight('bold')

# Make the last three y-axis labels bold
for label in ax.get_yticklabels()[-3:]:
    label.set_fontweight('bold')

plt.title("Spearman's Correlation Matrix for Metadata, MS$^2$ and Hazard Scores")
plt.show()






print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))  