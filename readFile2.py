# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:58:35 2020

@author: datamining-externe
"""

import pandas as pd

'''
read the csv file, with the selected columns. Read the dataframe in chucksize.

humanMethylation_TP is not a dataframe but pandas.io.parsers.TextFileReader
humanMethylation_df concat all chunks to df, bacuase of output of function.
It is necessary to add parameter ignore index to function concat, because avoiding 
duplicity

'''

humanMethylation_TP = pd.read_csv('HumanMethylation450_CHIP_INFO_LIMITED.csv',chunksize = 1000000, iterator=True, 
                                  error_bad_lines=False, skiprows = 7, 
 usecols = ['IlmnID', 'Chromosome_36', 'Coordinate_36', 'Strand', 'Probe_SNPs', 'UCSC_RefGene_Name',
            'UCSC_RefGene_Group', 'UCSC_CpG_Islands_Name', 'Relation_to_UCSC_CpG_Island', 'Enhancer'])


humanMethylation_df = pd.concat(humanMethylation_TP, ignore_index =True)
humanMethylation_df.reset_index #no double indexes on the right side of the csv files
humanMethylation_df.fillna("NA", inplace=True) # fills empty columns with "NA" instead of NaN

"""
Ask user for probe id(s) and search in dataframe.
.split(',') makes sure that all the values in the input have to be separated 
""" 
probe_id = [item.strip() for item in input("enter probe id(s) and separate with comma: ").split(',')]


def check_probe_id(probe_id):
     
        
    """"
    Check if given probe id is in a valid format
    """
    if any("cg" not in item for item in probe_id):
        print( "Please control your probe id(s), it should start with 'cg'.")
    elif any(len(item) !=10 for item in probe_id):
        print('The probe id(s) should be 10 characters.')
    return True

def get_probe(probe_id):
    
    """"
    Extract data from the row where IlmnID matches the probe_id variable
    If no result is found return "no results found" message.
    """
    
    Ilmn_ID = humanMethylation_df[humanMethylation_df['IlmnID'].isin(probe_id)]
    if Ilmn_ID.empty:
        print("No results found")
    else:
        return Ilmn_ID.to_string()

"""
If probe is in given format, and is in the IlmnID
print the probe and probe information
"""    
   
if check_probe_id(probe_id):
    result = get_probe(probe_id)
    #print(result)
