# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 16:34:43 2019

@author: datamining-externe
"""

import pandas as pd

#This is how you read the .csv file, with selected columns
df_450CHIP = pd.read_csv('HumanMethylation450_CHIP_INFO_LIMITED.csv', error_bad_lines=False, skiprows = 7,
 usecols = ['IlmnID', 'Chromosome_36', 'Coordinate_36', 'Strand', 'Probe_SNPs', 'UCSC_RefGene_Name',
            'UCSC_RefGene_Group', 'UCSC_CpG_Islands_Name', 'Relation_to_UCSC_CpG_Island', 'Enhancer'])

"""
With the reset_index function the index will not be shown double in the
csv file anymore
"""
df_450CHIP.reset_index
"""
With the fillna function, the empty spaces are going to be
replaced with 'NA'.
"""
df_450CHIP.fillna("NA", inplace=True)
#print(df_450CHIP.head(6))


"""
Ask user for probe id and search in data frame
""" 
probe_id = input('Please enter your Probe ID: ')

def check_probe_id(probe_id):
    """"
    Check if given probe id is in a valid format
    """
    if not probe_id.startswith("cg"):
        print( "The probe ID should start with lower case 'cg'.")
    elif len(probe_id) !=10:
        print( "The probe ID should be 10 characters.")
    return True


def get_probe(probe_id):
    """"
    Extract data from the row where IlmnID matches the probe_id variable
    If no result is found return "no results found" message.
    """
    Ilmn_ID = df_450CHIP.loc[df_450CHIP['IlmnID'] == probe_id]
    if Ilmn_ID.empty:
        return 'No results found'
    else:
        return Ilmn_ID.to_string()

"""
If probe is in given format, and is in the IlmnID
print the probe and probe information
"""
if check_probe_id(probe_id):
    result = get_probe(probe_id)
    #print(result)
