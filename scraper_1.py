# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:12:07 2020

@author: datamining-externe
"""

import requests
import xml.etree.ElementTree as ET
import readFile
from bs4 import BeautifulSoup
from urllib.request import urlopen

'''
Import readFile, from readFile using dataframe, in a new variable:
input_HMdf. In this Dataframe UCSC_RefGene_Name will be used to get information from HUGOdb.
df_450CHIP will be saved in the variable input_HMdf (input Human Methylation File).

from readFile, probe_id is also imported since that is the input, and the input
will be coupled to the input of the user and the data that will be fetched from
HUGO database.

After the input of probe_id the check_probe_id function is being imported since
the probe_id should get checked.


TODO: fix error: IndexError: index 0 is out of bounds for axis 0 with size 0
TODO: fix NoneTypeError
'''

input_HMdf = readFile.df_450CHIP
probe_id = readFile.probe_id
check_input = readFile.check_probe_id
input_HMdf.fillna("NA", inplace=True) 

'''
Get the column from the dataframe wich contains the probe_id in the IlmnID column
'''
input_HMdf = input_HMdf.loc[input_HMdf["IlmnID"] == probe_id]

'''
Get the UCSC_RefGene_Name values
'''
UGRN = input_HMdf["UCSC_RefGene_Name"].values[0]

'''
The values can now be used to search HugoDB, by ignoring ";"
'''
UGRN = UGRN.split(";")[0]



def connect_Hugodb(UGRN):
    
    url = 'http://rest.genenames.org/fetch/symbol/'+UGRN
    r = requests.get(url)
    r = r.text
    r = r.replace("\n", "")
    '''
    Make a xml object for later use
    '''
    xml = ET.fromstring(r)
    return xml,  url
    
def hugo_str_atrribs(x, UGRN):
    
    xml, url = connect_Hugodb(UGRN)
    
    try:
        for s in xml.find("result").find("doc").iter("str"):
            if x in s.attrib["name"]:
                x = s.text
    except KeyError:
        # Function does give a keyerror because the value name
        print()
    except AttributeError:
        print("There is no hgnc_id object. Maybe there is no RefGene_name")
    
    
    return x

def hugo_arr_attribs(y, UGRN):
    
    xml, url = connect_Hugodb(UGRN)
    
    try:
        for s in xml.find("result").find("doc").iter("arr"):
            if y in s.attrib["name"]:
                y = s.text
    except KeyError:
        # Function does give a keyerror because the value name
        print()
    except AttributeError:
        print("There is no hgnc_id object. Maybe there is no RefGene_name")
    
    return y

    '''
    Adding new columns to the dataframe: name, Locus_type, Alias_symbol, uniprot_ids
    
    New output is with these new 4 columns. In order to get the output, the objects
    are places in an Arraylist and the arraylist is placed in the dataframe input_HMdf.
    '''
    
    alias_arr = [alias_symbol]
    name_arr = [name]
    locus_arr = [locus_type]
    UP_arr = [uniprot_ids]
    
    input_HMdf['Alias_symbol'] = alias_arr
    input_HMdf['Appr_name'] = name_arr
    input_HMdf['Locus_type'] = locus_arr
    input_HMdf['uniprot_ids'] = UP_arr
    

def connect_uniProt_db():
    
    #uniprot_ids = get_attribs()
    
        
    UPID = input_HMdf["uniprot_ids"].values[0]
    '''
    Specify the URL (url_uniprot)
    '''
    url_uniprot = "http://www.uniprot.org/uniprot/"+ str(UPID)
    
    def fetchdata(url_uniprot):
        try:
            return urlopen(url_uniprot)
        except:
            fetchdata(url_uniprot)
    fetchdata(url_uniprot)
   
    '''
    Query the website
    '''
    page = fetchdata(url_uniprot)
    #print(page)
   
    '''
    Parse the html, store it in Beautiful Soup format
    '''
    bsf = BeautifulSoup(page, "lxml")
    return bsf
    
def parse_uniprot_info():
    
    bsf = connect_uniProt_db()
    
    '''
    Parse the # get Gene ID
    '''
    gene_id = ""
    ext_data = bsf.find('table', class_ = 'databaseTable GENOME')
    if(not(ext_data is None)):
        data = ext_data.find_all(text=True)
        if "GeneID" in data:
            i = data.index("GeneID")
            gene_id = data[i+2]
        #print(gene_id)
   
    '''
    get Pathway
    '''
    getPW = []
    ext_data = bsf.find('table', class_= 'databaseTable PATHWAY')
    if(not(ext_data is None)):
        data = ext_data.find_all(text=True)
        if "Reactome" in data:
            i = data.index("Reactome")
            getPW.append(data[4]+data[5])
            print(getPW)
    
    
    '''
    Molecular Function GO indeannotation
    '''
    molecular_function_go = []
    ext_data = bsf.find('ul', class_ ='noNumbering molecular_function')
    if(not(ext_data is None)):
        for data in ext_data.find_all('li'):
            cells = data.find(lambda tag: tag.name =="a" and (tag.has_attr('onclick')))
            if(not(cells is None)):
                cells = cells.find(text=True).strip()
                molecular_function_go.append(cells)
    #print(molecular_function_go)
    
    '''
    biologoical process GO annotation
    '''
    biological_process_go = []
    ext_data = bsf.find('ul', class_='noNumbering biological_process')
    if(not(ext_data is None)):
        for data in ext_data.find_all('li'):
            cells = data.find(lambda tag: tag.name == "a" and (tag.has_attr("onclick")))
            if(not(cells is None)):
                cells = cells.find(text=True).strip()
                biological_process_go.append(cells)
        #print(biological_process_go)
        
    '''
    cellular Component GO annotation
    '''
    cellular_component_go = []
    ext_data = bsf.find('div', id='table-go_annotation')
    if(not(ext_data is None)):
        lt = ext_data.find('ul', class_='noNumbering subcellLocations')
        for li in lt.find_all('li'):
            cell = li.find('ul')
            if(not(cell is None)):
                cell_loc = cell.find(text=True)
                cellular_component_go.append(cell_loc)
        #print(cellular_component_go)
        
    '''
    cellular Component Keywords
    '''
    cellular_component_kw = []
    ext_data = bsf.find('div', class_='section', id='subcellular_location')
    if(not(ext_data is None)):
        header = ext_data.find('h4')
        if(not(header is None)):
            head_data = header.find_all(text=True)
            if 'Keywords - Cellular component' in head_data:
                data = header.next_sibling
                vals = data.find_all('a')
                val = [v.find(text=True) for v in vals]
                v = list(filter(lambda x :  x != ',', vals))
                val = []
                for kw in v:
                    kws = str(kw)
                    if kws[:18]=='<a href="/keywords':
                        val.append(kw.find(text=True))
                cellular_component_kw = val
                
                                    
    '''
    keywords - Molecular function and Biological process
    '''
    molecular_function_kw = []
    biological_process_kw = []
    ext_data = bsf.find('table', class_='databaseTable')
    if(not(ext_data is None)):
        ext_data = ext_data.find_all('tr')
        for row in ext_data:
            data = row.find_all('td')
            head = data[0].find(text=True)
            vals = data[1].find_all('a')
            val = [v.find(text=True) for v in vals]
            v = list(filter(lambda x : x != ', ', vals))
            val = []
            for kw in v:
                kws = str(kw)
                if kws[:18]=='<a href="/keywords':
                    val.append(kw.find(text=True))
            if(head == 'Molecular function'):
                molecular_function_kw = val
            if(head == 'Biological process'):
                biological_process_kw = val
                
    '''
    Get disease summary
    '''
    disease_info = []
    ext_data = bsf.find('div', class_='section', id='pathology_and_biotech')
    if(not(ext_data is None)):
        disTxt = ext_data.find('div', class_='diseaseAnnotation')
        for info in disTxt:
            data = info.find_all(text=True)
            disease_info.append(data[0])
    
    '''
    Get text tissue specificity
    '''
    tis_specificity = []
    ext_data = bsf.find('div', class_='section', id = 'expression')
    if(not(ext_data is None)):
        tisText = ext_data.find_all('div', class_='annotation')
        for info in tisText:
            data = info.find_all(text=True)
            tis_specificity.append(data[0])
    
    
    

    
    '''
    Adding new columns to the dataframe: name, Locus_type, Alias_symbol, uniprot_ids
    
    New output is with these new 7 columns. In order to get the output, the objects
    are placed in an Arraylist and the arraylist is placed in the dataframe input_HMdf.
    Next to the HUGO database output columns
    '''            
    geneid_arr = [gene_id]
    getPW_arr = [getPW]
    mol_func_arr = [molecular_function_go]
    mfkw_arr = [molecular_function_kw]
    bio_proc_arr = [biological_process_go]
    bpkw_arr = [biological_process_kw]
    cell_comp_arr = [cellular_component_go]
    cckw_arr = [cellular_component_kw]
    disInfo_arr = [disease_info]
    tisSpec_arr = [tis_specificity]
    
    #print(len(gid_arr, getPw_arr))
    
    
    '''
    Add new columns in the dataframe
    '''
    input_HMdf['GeneID'] = geneid_arr
    input_HMdf['Pathway'] = getPW_arr
    input_HMdf['Molecular_Function'] = mol_func_arr
    input_HMdf['Mol_func_keywords'] = mfkw_arr
    input_HMdf['Biological_Process'] = bio_proc_arr
    input_HMdf['Bio_proc_keywords'] = bpkw_arr
    input_HMdf['Cellular component'] = cell_comp_arr
    input_HMdf['Cel_comp_keywords'] = cckw_arr
    input_HMdf['Disease_involvement'] = disInfo_arr
    input_HMdf['Tissue Specificity'] = tisSpec_arr
       
connect_Hugodb(UGRN)
hugo_str_atrribs('name', UGRN)
hugo_str_atrribs('locus_type', UGRN)
hugo_arr_attribs('alias_symbol', UGRN)
hugo_arr_attribs('uniprot_ids', UGRN)
connect_uniProt_db()
parse_uniprot_info()


print(input_HMdf)