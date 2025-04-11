# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 11:16:08 2024

@author: AChao
"""

import requests
import json
import pandas as pd
import time

time_list = []
start_time = time.time()
time_list.append(time.time())

# Read in WW2DW results - Chemicals
df1 = pd.read_excel('MZmine3_cal_tracer_NTA_WebApp_results.xlsx', sheet_name='chemical_results')

# Read in the matched SSM chemicals when ran as tracers
df2 = pd.read_excel('MZmine3_cal_SSM_NTA_WebApp_results.xlsx', sheet_name='Tracers_Summary')

# Mark the correct chemicals
df2['SSM_chemical'] = 'Y'
df2 = df2[['Feature_ID', 'DTXSID', 'SSM_chemical']]
SSM_features = df2['Feature_ID'].tolist()

# Grab just the chemicals corresponding to 63 SSM chemicals
df3 = df1[df1['Feature_ID'].isin(SSM_features)]
df3 = pd.merge(df3, df2, how='left', on=['Feature_ID', 'DTXSID'])

# Export the list of chemicals to pull back ALL metadata from DSSTox API
df3_chemicals = df3['DTXSID'].tolist()


headers = {
    'x-api-key': '18358de8-3355-4511-9ed7-ddb68f013a7c',
    'accept': 'application/json',
    'content-type': 'application/json',
}

params = {
    'projection': 'ntatoolkit',
}


json_data = [
    'DTXSID7020182',
    'DTXSID9020112',
    'DTXSID50382005',
    'DTXSID7020185',
    'DTXSID9020111',
    'DTXSID50382015',
    'DTXSID1034343',
    'DTXSID6091554',
    'DTXSID7021360',
    'DTXSID0020446',
    'DTXSID7026073',
    'DTXSID10274073',
    'DTXSID5049817',
    'DTXSID5020079',
    'DTXSID7064819',
    'DTXSID90896818',
    'DTXSID8020044',
    'DTXSID4039231',
    'DTXSID6091556',
    'DTXSID3021982',
    'DTXSID50931532',
    'DTXSID6022977',
    'DTXSID101022327',
    'DTXSID0049650',
    'DTXSID3020122',
    'DTXSID90892315',
    'DTXSID7025683',
    'DTXSID1064827',
    'DTXSID1027263',
    'DTXSID50884938',
    'DTXSID3042427',
]




# Limit API calls to batches of 1000 DTXSID's
batch_size = 1000
for i in range(0, len(df3_chemicals), batch_size):
    #print(df3_chemicals[i:i+batch_size])
    chemical_batch = df3_chemicals[i:i+batch_size]
    #print(i)
    response = requests.post(
        'https://api-ccte.epa.gov/chemical/detail/search/by-dtxsid/',
        #params=params,
        headers=headers,
        json=chemical_batch,
        #json=json_data,
    )

    api_results = response.json()

    if i == 0:
        df4 = pd.DataFrame(api_results)
    else:
        df4 = pd.concat([df4, pd.DataFrame(api_results)])
    
df4.rename(columns={'dtxsid':'DTXSID'}, inplace=True) 
df5 = pd.merge(df3, df4, how='left', on='DTXSID')

'''
response = requests.post(
    'https://api-ccte.epa.gov/chemical/detail/search/by-dtxsid/',
    #params=params,
    headers=headers,
    json=df3_chemicals,
    #json=json_data,
)


api_results = response.json()

df4 = pd.DataFrame(api_results)
'''

# Note: json_data will not be serialized by requests
# exactly as it was in the original request.
#data = '["DTXSID7020182","DTXSID9020112", "DTXSID50382005","DTXSID7020185","DTXSID9020111", "DTXSID50382015","DTXSID1034343","DTXSID6091554","DTXSID7021360","DTXSID0020446","DTXSID7026073","DTXSID10274073","DTXSID5049817","DTXSID5020079","DTXSID7064819","DTXSID90896818","DTXSID8020044","DTXSID4039231","DTXSID6091556","DTXSID3021982","DTXSID50931532","DTXSID6022977","DTXSID101022327","DTXSID0049650","DTXSID3020122","DTXSID90892315","DTXSID7025683","DTXSID1064827","DTXSID1027263","DTXSID50884938","DTXSID3042427"]'
#response = requests.post('https://api-ccte.epa.gov/chemical/detail/search/by-dtxsid/', params=params, headers=headers, data=data)



print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))  