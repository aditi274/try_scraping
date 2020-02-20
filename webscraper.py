# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 16:17:33 2020

@author: datamining-externe
"""

import readFile2
import requests

'''
Use dataframe from module readFile.py
'''
use_of_df = readFile2.humanMethylation_df
probe_id = readFile2.probe_id
check_probes = readFile2.check_probe_id
use_of_df.fillna("NA", inplace=True)
#get_df_output = readFile2.get_probe

use_of_df = use_of_df[use_of_df['IlmnID'].isin(probe_id)]
'''
Use UCSC_RefGene_Name column for scraping Hugo database
'''
#UGRN = use_of_df.loc[use_of_df["UCSC_RefGene_Name"] == probe_id].iloc[0]
UGRN = use_of_df['UCSC_RefGene_Name'].values[0:]
UGRN = [i.split(';')[0] for i in UGRN]
#print(UGRN)
#print(use_of_df)

def connect_HDB(UGRN):
    
    url = 'http://rest.genenames.org/fetch/symbol/'+str(UGRN)
    request = requests.get(url)
    request = request.text
    print(request)
    
connect_HDB(UGRN)