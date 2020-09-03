# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 11:55:11 2019

@author: smith.j
"""

import pandas as pd
from pathlib import Path
import os
import requests, sys
import io
from gprofiler import GProfiler
from fake_useragent import UserAgent
import jwrangle

def mapAnyID_gPro(IdList, splitstr = ['nosplit'], geneProductType = 'protein', gConvertOrganism = 'hsapiens', gConvertTarget = 'ENTREZGENE', writetopath = 'NO', writeTargetsAsList = 'NO'):
    '''
    Accepts a list of individual IDs, or a list of grouped IDs, and returns a dataframe of remapped IDs. In addition to the original list, 
    these remapped IDs are represented by 3 columns:
    target_all = all the mapped IDs found for each ID or group of IDs
    target_primary = a single ID selected based firstly on swissprot membership then, if unavailable, trembl membership. If there is a group of IDs, then
    only the first will be chosen from the sorted group (more detail under 'Limitations'). This is done, so that string comparisons can be made between 
    experiments by consistently named genes.
    target_name = the gene name matching the target primary.
    
    Input:
    IdList = a list of uniprot IDs or a list of groups of uniprot IDs (i.e. example of the latter being 'Majority protein group IDs')
    IdType = a string representing the type of biomolecule being searched, i.e. 'protein', 'rna'
    splitstr = a list of charcters to split a group of IDs upon. If = 'nosplit' then no split will be done (default)
    gConvertTarget ='ENTREZGENE'                                   #see gProfiler webiste or PyPi page for options
    gConvertOrganism = 'hsapiens'                                  #see gProfiler webiste or PyPi page for info
    writetopath = a list where save[0] = cwd in pathlib format, and save[1] = a subfolder

    Return: A dictionary including
    id_map = A dataframe of the original list with the 3 additional columns
    query_map = A dataframe of the gProfiler query and results
    
    Notes:
    Writetopath is designed to synergise with jwrangle.importMixedFiles() such that a dictionary object resenbling the output of 
    jweb.mapAnyID_gPro() can easily be retrieved from the filesystem.
    
    Limitations:
    'None' is used where no string can be found. This may affect the assignment of neuraminidase for Influenza A virus.
    If a SYMBOL type is selected and a putative or predicted gene is converted for primary as signified by 'LOC' then the non-LOC gene names will be given sorting
    preference when selecting a primary name. Where LOC status is hidden, i.e. by gene acc. then this sorting will not occur. Overall, however, sorting is an 
    imprecise user-intervention that I've intended only for use in comparisons that depend on set theory (i.e. we are forced to select only one label for each protein group
    therefore sort is a means of doing so consistently between MQ data). Yes, how incredibly complicated.
    '''

    gpquery = jwrangle.SplitList(IdList, splitstr, replace = '$&', remove_nonstring_values = True, drop_empty_strings = True)
    gpquery = [i.replace('CON__','') for i in gpquery]

    ua = UserAgent()
    # gplist = GProfiler(ua.chrome).gconvert(query = gpquery, organism = gConvertOrganism, target = gConvertTarget, region_query=False, numeric_ns=None)        ##gprofiler v0.3.5
    gplist = GProfiler(ua.chrome).convert(query = gpquery, organism = gConvertOrganism, target_namespace = gConvertTarget) 
    gp_df = pd.DataFrame(gplist).replace(regex=[gConvertTarget+':'], value='')
        
    gPro_all = []
    gPro_primary = []
    gPro_name = []
    
    if geneProductType != 'protein':
        for entry in IdList:
            ## Check if list iterables are ID groups or single
            for i in splitstr:
                entry = entry.replace(i,'$&')
                mplist = entry.split('$&')
            ## replace np.nan so strings can be indexed for searching
            # hits = gp_df[gp_df[1].isin(mplist)].replace(np.nan, 'na')                                                                                            ##gprofiler v0.3.5
            # hits = gp_df[gp_df['incoming'].isin(mplist)].replace(np.nan, 'na')                                                                                   ##gprofiler v0.3.5
            hits = gp_df[gp_df['incoming'].isin(mplist)].replace('None', 'None') #np.nan > 'None', na' > 'None' (update to v1.0)


            if len(hits) > 0:
                ## sort hits for primary selection
                # primaryhits_list = sorted(hits[3].tolist(), key= lambda x:('LOC' in x, x))                                                                       ##gprofiler v0.3.5
                primaryhits_list = sorted(hits['converted'].tolist(), key= lambda x:('LOC' in x, x))  

                gPro_primary.append(primaryhits_list[0])            
                # gPro_name.append(hits.loc[hits[3]==primaryhits_list[0], 5].values[0])                                                                            ##gprofiler v0.3.5
                gPro_name.append(hits.loc[hits['converted']==primaryhits_list[0], 'description'].values[0])

                
                ## pool entries for target_all
                # listhits = hits[3].tolist()                                                                                                                      ##gprofiler v0.3.5
                listhits = hits['converted'].tolist()

                trimlisthits = []
                for i in listhits:
                    if type(i) == str and i not in trimlisthits:
                        trimlisthits.append(i)
                trimlisthits.sort()
                if 'na' in trimlisthits and len(trimlisthits) > 1:
                    trimnahits = [i for i in trimlisthits if i != 'None']  #np.nan > 'None', na' > 'None' (update to v1.0)
                    gPro_all.append(';'.join(trimnahits))
                else:    
                    gPro_all.append(';'.join(trimlisthits))
        
            elif len(hits) == 0:
                gPro_all.append('None')  #np.nan > 'None', na' > 'None' (update to v1.0)
                gPro_primary.append('None')  #np.nan > 'None', na' > 'None' (update to v1.0)
                gPro_name.append('None')  #np.nan > 'None', na' > 'None' (update to v1.0)
        
        newdf = pd.DataFrame({'Query': IdList, gConvertTarget + '_gPro all': gPro_all, gConvertTarget + '_gPro primary': gPro_primary, 
                              gConvertTarget + '_gPro name': gPro_name})
 
    else:
        gPro_protein_status = []
        
        for entry in IdList:
            ## Check if list iterables are ID groups or single
            for i in splitstr:
                entry = entry.replace(i,'$&')
                mplist = entry.split('$&')
            ## replace np.nan so strings can be indexed for searching
            # hits = gp_df[gp_df[1].isin(mplist)].replace(np.nan, 'na')                                                                                             ##gprofiler v0.3.5
            # hits = gp_df[gp_df['incoming'].isin(mplist)].replace(np.nan, 'na')                                                                                    ##gprofiler v0.3.5
            hits = gp_df[gp_df['incoming'].isin(mplist)].replace('None', 'None')  #np.nan > 'None', na' > 'None' (update to v1.0)

    
            ## pool entries for target_all
            if len(hits) > 0:
                # listhits = hits[3].tolist()
                listhits = hits['converted'].tolist()
                trimlisthits = []
                for i in listhits:
                    if type(i) == str and i not in trimlisthits:
                        trimlisthits.append(i)
                trimlisthits.sort()
                if 'na' in trimlisthits and len(trimlisthits) > 1:
                    trimnahits = [i for i in trimlisthits if i != 'None']  #np.nan > 'None', na' > 'None' (update to v1.0)
                    gPro_all.append(';'.join(trimnahits))
                else:    
                    gPro_all.append(';'.join(trimlisthits))
        
            elif len(hits) == 0:
                gPro_all.append('None')  #np.nan > 'None', na' > 'None' (update to v1.0)
            
            ## filter entries for swissprot vs sptrembl origin
            # sphits = hits[hits[6].str.contains('UNIPROTSWISSPROT')]                                                                                               ##gprofiler v0.3.5
            sphits = hits[hits['namespaces'].str.contains('UNIPROTSWISSPROT')]
            
            # otherhits = hits[(hits[6].str.contains('UNIPROTSPTREMBL')) | (hits[6].str.contains('UNIPROT_GN'))]                                                    ##gprofiler v0.3.5                                                                                                              ##gprofiler v0.3.5
            otherhits = hits[(hits['namespaces'].str.contains('UNIPROTSPTREMBL')) | (hits['namespaces'].str.contains('UNIPROT_GN'))]
            
            ## check for swissprot values first, sort then choose the first
            if len(sphits) > 0:
                # primaryhits_list = sorted(sphits[3].tolist(), key= lambda x:('LOC' in x, x))
                primaryhits_list = sorted(sphits['converted'].tolist(), key= lambda x:('LOC' in x, x))

                gPro_primary.append(primaryhits_list[0])            
                
                # gPro_name.append(hits.loc[hits[3]==primaryhits_list[0], 5].values[0])
                gPro_name.append(hits.loc[hits['converted']==primaryhits_list[0], 'description'].values[0])
                
                
                gPro_protein_status.append('SWISSPROT')
            ## if no swissprot entires then check for trembl values second, sort then choose the first
            elif len(sphits) == 0:
                if len(otherhits) > 0:
                    # primaryhits_list = sorted(otherhits[3].tolist(), key= lambda x:('LOC' in x, x))                                                              ##gprofiler v0.3.5
                    primaryhits_list = sorted(otherhits['converted'].tolist(), key= lambda x:('LOC' in x, x))
                    
                    gPro_primary.append(primaryhits_list[0])            
                    
                    # gPro_name.append(hits.loc[hits[3]==primaryhits_list[0], 5].values[0])                                                                        ##gprofiler v0.3.5
                    gPro_name.append(hits.loc[hits['converted']==primaryhits_list[0], 'description'].values[0])

                    gPro_protein_status.append('TREMBL')
            ## if no entries can be found at all add 'na'
                else:
                    gPro_primary.append('None')  #np.nan > 'None', na' > 'None' (update to v1.0)
                    gPro_name.append('None')  #np.nan > 'None', na' > 'None' (update to v1.0)
                    gPro_protein_status.append('None')  #np.nan > 'None', na' > 'None' (update to v1.0)
    
        newdf = pd.DataFrame({'Query': IdList, gConvertTarget + '_gPro all': gPro_all, gConvertTarget + '_gPro primary': gPro_primary, 
                              gConvertTarget + '_gPro name': gPro_name, 'UNIPROT_gPro status': gPro_protein_status})
    
    # gp_df = pd.DataFrame(gplist).replace(np.nan,'na')                                                                                                            ##gprofiler v0.3.5
    gp_df = pd.DataFrame(gplist).replace('None','None')  #np.nan > 'None', na' > 'None' (update to v1.0)

    if writetopath != 'NO':

        if writeTargetsAsList == 'NO':
            p = Path(writetopath[0] / 'downloads' / writetopath[1])
            if p.exists() == False:
                Path(writetopath[0] / 'downloads' / writetopath[1]).mkdir(parents=True, exist_ok=True)

            newdf.to_csv(writetopath[0] / 'downloads' / writetopath[1] / 'id_map.csv', index = False)
            gp_df.to_csv(writetopath[0] / 'downloads' / writetopath[1] / 'query_map.csv', index = False)
        
        else:
            p = Path(writetopath[0] / 'outputs')
            if p.exists() == False:
                Path(writetopath[0] / 'outputs').mkdir(parents=True, exist_ok=True)
            
            listname = writetopath[-1] + '_targetList' + '.txt'
            with open(writetopath[0] / 'outputs' / listname, 'w') as f:

                for item in gPro_primary:
                    f.write("%s " % item)
    
    return {'id_map':newdf, 'query_map':gp_df}


def fetchQuickGO(QG_goId, QG_geneProductType = 'protein', QG_taxonId = '9606', QG_geneProductSubset = ['Swiss-Prot', 'TrEMBL'], 
                 gConvertTarget='ENTREZGENE', gConvertOrganism = 'hsapiens', writetopath = 'NO'):
    '''
    Sends a list of GO codes to the QuickGo server and parses the response element into a dictionary (df_dict) where keys = GO query and value = html response
    as a DataFrame. Also takes each DataFrame in df_dict and builds a second dictionary (df_lists) where keys = GO query, values = lists of unique values for 
    the specified gene product type. Finally, sends each value in df_lists to GProfiler and converts each ID into a specified convention (i.e. ENTREZGENE)
    then creates a third dictionary (gConvert) where keys = GO query and value = converted IDs. Useful for fetching up-to-date GO elements: to preserve list 
    orders for dataframe annotation see the function gConvertFromDictionary().
    
    Based on QuickGo programmatic access and gProfiler, https://pypi.org/project/gprofiler-official/
    QG_goId = ["GO:0003723"]                                       #queries QuickGo webservice for GO terms in a list
    QG_geneProductType = 'protein'                                 #choose from 'miRNA', 'complex' or 'protein'
    QG_taxonId = '9606'                                            #9606 = human
    QG_geneProductSubset = ['Swiss-Prot', 'TrEMBL']                #default = 'None' (for rna and complexes search).
                                                                    for protein searches; choose from "Swiss-Prot" or "TrEMBL" or get both
    gConvertTarget ='ENTREZGENE'                                   #see gProfiler webiste or PyPi page for options
    gConvertOrganism = 'hsapiens'                                  #see gProfiler webiste or PyPi page for info
    
    Returns: A dictionary whereby
    query_map: df_dict =  a dictionary where keys = GO query and value = html response as a DataFrame
    QG_geneProductType + '_lists': df_lists = a dictionary where keys = GO query, values = lists of unique values for the specified gene product type (i.e. rna, complexes or protein)
    gConvertTarget + '_lists': gConvert = a dictionary where keys = GO query and value = converted IDs
    Yes that's right, a dictionary within a dictionary....absolute madness.
    
    Limitations:
    API won't allow more than 10000 records per automated query. Must do it manually where records >10000.
    Some options are hardcoded for the sake of simplicity. See Source to modify.
    A manual QuickGO download can retrieve the gene symbol- I'm not sure how this can be achived with the REST API (not present in this function
    '''
    df_dict = {}
    for value in QG_goId:
        
        if QG_geneProductSubset != ['None']:
        
            for entry in QG_geneProductSubset:
                # Build URL for protein subgroups
                requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?" + \
                                                                                                        "goId=" + value +\
                                                                                                        "&taxonId=" + QG_taxonId +\
                                                                                                        "&taxonUsage=exact&goUsage=descendants" +\
                                                                                                        "&geneProductType=" + QG_geneProductType +\
                                                                                                        "&geneProductSubset=" + entry
                response = requests.get(requestURL, headers={ "Accept" : "text/gpad"})
        
                if not response.ok:
                  response.raise_for_status()
                  sys.exit()
                responseBody = response.text
        
                quick_go_df = pd.read_csv(io.StringIO(responseBody), skiprows = 9, header = None, sep='\t', names= \
                               ['GENE_PRODUCT_DB', 'GENE_PRODUCT_ID', 'QUALIFIER_QGO', 'GO_TERM', 'REFERENCE', 'ECO_ID', 'WITH/FROM',\
                                'INTERACTING_TAXON','DATE','ASSIGNED_BY','ANNOTATION_EXTENSION_QGO','GO_EVIDENCE_CODE'])
                
                quick_go_df['GO_EVIDENCE_CODE'] = quick_go_df['GO_EVIDENCE_CODE'].apply(lambda x: x.replace('goEvidence=',''))
                quick_go_df = quick_go_df[['GENE_PRODUCT_DB', 'GENE_PRODUCT_ID', 'QUALIFIER_QGO', 'GO_TERM', 'REFERENCE', 'ECO_ID', 
                                           'WITH/FROM', 'ASSIGNED_BY','ANNOTATION_EXTENSION_QGO','GO_EVIDENCE_CODE']]
                
                manualURL = requestURL.replace('services/annotation/downloadSearch?', 'annotations?')
                
                if len(quick_go_df) == 10000:
                    print(value.replace(':','') + '_' + entry + ' has exceeded 10000 record limit, a manual download will be required: \n' + manualURL + '\n')
                else:
                    record_no = str(len(quick_go_df))
                    print(value.replace(':','') + '_' + entry + ' has ' + record_no + ' records: \n' + manualURL + '\n')
                          
                df_dict[value.replace(':','') + '_' + entry] = quick_go_df
                
        else:
            # Build URL where no subgroups are wanted....I can't be bothered avoiding the code duplication.
            requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?" + \
                                                                                                    "goId=" + value +\
                                                                                                    "&taxonId=" + QG_taxonId +\
                                                                                                    "&taxonUsage=exact&goUsage=descendants" +\
                                                                                                    "&geneProductType=" + QG_geneProductType
            response = requests.get(requestURL, headers={ "Accept" : "text/gpad"})
    
            if not response.ok:
              response.raise_for_status()
              sys.exit()
            responseBody = response.text
    
            quick_go_df = pd.read_csv(io.StringIO(responseBody), skiprows = 9, header = None, sep='\t', names= \
                               ['GENE_PRODUCT_DB', 'GENE_PRODUCT_ID', 'QUALIFIER_QGO', 'GO_TERM', 'REFERENCE', 'ECO_ID', 'WITH/FROM',\
                                'INTERACTING_TAXON', 'ASSIGNED_BY','ANNOTATION_EXTENSION_QGO','GO_EVIDENCE_CODE'])
            
            quick_go_df['GO_EVIDENCE_CODE'] = quick_go_df['GO_EVIDENCE_CODE'].apply(lambda x: x.replace('goEvidence=',''))
            quick_go_df = quick_go_df[['GENE_PRODUCT_DB', 'GENE_PRODUCT_ID', 'QUALIFIER_QGO', 'GO_TERM', 'REFERENCE', 'ECO_ID', 
                                       'WITH/FROM', 'DATE','ASSIGNED_BY','ANNOTATION_EXTENSION_QGO','GO_EVIDENCE_CODE']]
            
            manualURL = requestURL.replace('services/annotation/downloadSearch?', 'annotations?')
     
            if len(quick_go_df) == 10000:
                manualURL = requestURL.replace('%3A', ':')
                manualURL = manualURL.replace('services/annotation/downloadSearch?', 'annotations?')
                print(value.replace(':','') + '_' + QG_geneProductType + ' has exceeded 10000 record limit, a manual download will be required: \n' + manualURL + '\n')
            else:
                record_no = str(len(quick_go_df))
                print(value.replace(':','') + '_' + QG_geneProductType + ' has ' + record_no + ' records: \n' + manualURL + '\n')

            df_dict[value.replace(':','') + '_' + QG_geneProductType] = quick_go_df

    df_dict_merge = {}
    for key, value in df_dict.items():
        gPro_query = value['GENE_PRODUCT_ID'].tolist()
        gPro_dict = mapAnyID_gPro(gPro_query, splitstr = ['nosplit'], geneProductType = QG_geneProductType, gConvertOrganism = gConvertOrganism, gConvertTarget = gConvertTarget)
        gPro_df = gPro_dict['id_map'].rename(columns={'Query': 'GENE_PRODUCT_ID'}).drop_duplicates()
        gPro_df['GENE_PRODUCT_ID'] = gPro_df['GENE_PRODUCT_ID'].astype(str)
        value['GENE_PRODUCT_ID'] = value['GENE_PRODUCT_ID'].astype(str)
        query_map = pd.merge(value, gPro_df, on='GENE_PRODUCT_ID', how='left' )
        df_dict_merge[key] = query_map

    cwd  = Path(os.getcwd())
    p = Path(cwd / 'downloads' / 'QuickGo')
    if p.exists() == False:
        Path(cwd / 'downloads' / 'QuickGo').mkdir(parents=True, exist_ok=True)

    if writetopath != 'NO':
        for key, value in df_dict_merge.items():
            fname = key + '_' + gConvertOrganism + '_' + QG_geneProductType + '-' + gConvertTarget + '.csv'
            value.to_csv(cwd / 'downloads' / 'QuickGo' / fname, index = False)
            
    return df_dict_merge

def fetchQuickGO_stats(QG_goId, QG_geneProductType = 'protein', QG_taxonId = '9606', QG_geneProductSubset = ['Swiss-Prot', 'TrEMBL']):
    '''
    Sends a list of GO codes to the QuickGo server and parses the response element into a dictionary (df_dict) where keys = GO query and value = html response
    as a dictionary of values or dataframes. The value holds all QuickGo stats related to the selected GO ID.
    
    Based on QuickGo programmatic access and gProfiler, https://pypi.org/project/gprofiler-official/
    QG_goId = ["GO:0003723"]                                       #queries QuickGo webservice for GO terms in a list
    QG_geneProductType = 'protein'                                 #choose from 'miRNA', 'complex' or 'protein'
    QG_taxonId = '9606'                                            #9606 = human
    QG_geneProductSubset = ['Swiss-Prot', 'TrEMBL']                #default = 'None' (for rna and complexes search).
                                                                    for protein searches; choose from "Swiss-Prot" or "TrEMBL" or get both
    gConvertTarget ='ENTREZGENE'                                   #see gProfiler webiste or PyPi page for options
    gConvertOrganism = 'hsapiens'                                  #see gProfiler webiste or PyPi page for info
    
    Returns: A dictionary whereby
    query_map: df_dict =  a dictionary where keys = GO query and value = html response as a DataFrame
    QG_geneProductType + '_lists': df_lists = a dictionary where keys = GO query, values = lists of unique values for the specified gene product type (i.e. rna, complexes or protein)
    gConvertTarget + '_lists': gConvert = a dictionary where keys = GO query and value = converted IDs
    Yes that's right, a dictionary within a dictionary....absolute madness.
    
    Limitations:
    Some options are hardcoded for the sake of simplicity. See Source to modify.
    '''

    df_dict = {}
    for value in QG_goId:
        
        if QG_geneProductSubset != ['None']:
        
            for entry in QG_geneProductSubset:
                # Build URL for protein subgroups
                requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/stats?" + \
                                                                                                        "goId=" + value +\
                                                                                                        "&taxonId=" + QG_taxonId +\
                                                                                                        "&taxonUsage=exact&goUsage=descendants" +\
                                                                                                        "&geneProductType=" + QG_geneProductType +\
                                                                                                        "&geneProductSubset=" + entry

                r = requests.get(requestURL)
                
                if not r.ok:
                  r.raise_for_status()
                  sys.exit()
                
                responseBody = r.json()
                
                if responseBody['numberOfHits'] > 0:
                
                    annotations_Hits = responseBody['results'][0]['totalHits']
                    genesProducts_Hits = responseBody['results'][1]['totalHits']
                    
                    annotation_Refs = pd.DataFrame.from_dict(responseBody['results'][0]['types'][0]['values'])
                    genesProduct_Refs = pd.DataFrame.from_dict(responseBody['results'][1]['types'][0]['values'])
                    
                    annotation_AssignedBy = pd.DataFrame.from_dict(responseBody['results'][0]['types'][1]['values'])
                    genesProduct_AssignedBy = pd.DataFrame.from_dict(responseBody['results'][1]['types'][1]['values'])
                    
                    annotation_HitStats = pd.DataFrame.from_dict(responseBody['results'][0]['types'][2]['values'])
                    genesProduct_HitStats = pd.DataFrame.from_dict(responseBody['results'][1]['types'][2]['values'])
                    
                    annotation_Taxon = pd.DataFrame.from_dict(responseBody['results'][0]['types'][3]['values'])
                    genesProduct_Taxon = pd.DataFrame.from_dict(responseBody['results'][1]['types'][3]['values'])
                    
                    annotation_ECOID = pd.DataFrame.from_dict(responseBody['results'][0]['types'][4]['values'])
                    genesProduct_ECOID = pd.DataFrame.from_dict(responseBody['results'][1]['types'][4]['values'])
                    
                    annotation_Aspect = pd.DataFrame.from_dict(responseBody['results'][0]['types'][5]['values'])
                    genesProduct_Aspect = pd.DataFrame.from_dict(responseBody['results'][1]['types'][5]['values'])
                    
                    GOstats =   {'annotations_Hits':annotations_Hits, 'genesProducts_Hits':genesProducts_Hits,
                                 'annotation_Refs':annotation_Refs, 'genesProduct_Refs':genesProduct_Refs,
                                 'annotation_AssignedBy':annotation_AssignedBy, 'genesProduct_AssignedBy':genesProduct_AssignedBy,
                                 'annotation_HitStats':annotation_HitStats,'genesProduct_HitStats':genesProduct_HitStats,
                                 'annotation_Taxon':annotation_Taxon, 'genesProduct_Taxon':genesProduct_Taxon,
                                 'annotation_ECOID':annotation_ECOID, 'genesProduct_ECOID':genesProduct_ECOID,
                                 'annotation_Aspect':annotation_Aspect, 'genesProduct_Aspect':genesProduct_Aspect}
                    
                    df_dict[value.replace(':','') + '_' + entry] = GOstats
                
                elif responseBody['numberOfHits'] == 0:
                    df_dict[value.replace(':','') + '_' + entry] = {'annotations_Hits':0, 'genesProducts_Hits':0}

        else:
            # Build URL where no subgroups are wanted....I can't be bothered avoiding the code duplication.
            requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/stats?" + \
                                                                                                    "goId=" + value +\
                                                                                                    "&taxonId=" + QG_taxonId +\
                                                                                                    "&taxonUsage=exact&goUsage=descendants" +\
                                                                                                    "&geneProductType=" + QG_geneProductType
            r = requests.get(requestURL)
            
            if not r.ok:
              r.raise_for_status()
              sys.exit()
            
            responseBody = r.json()
            
            if responseBody['numberOfHits'] > 0:
            
                annotations_Hits = responseBody['results'][0]['totalHits']
                genesProducts_Hits = responseBody['results'][1]['totalHits']
                
                annotation_Refs = pd.DataFrame.from_dict(responseBody['results'][0]['types'][0]['values'])
                genesProduct_Refs = pd.DataFrame.from_dict(responseBody['results'][1]['types'][0]['values'])
                
                annotation_AssignedBy = pd.DataFrame.from_dict(responseBody['results'][0]['types'][1]['values'])
                genesProduct_AssignedBy = pd.DataFrame.from_dict(responseBody['results'][1]['types'][1]['values'])
                
                annotation_HitStats = pd.DataFrame.from_dict(responseBody['results'][0]['types'][2]['values'])
                genesProduct_HitStats = pd.DataFrame.from_dict(responseBody['results'][1]['types'][2]['values'])
                
                annotation_Taxon = pd.DataFrame.from_dict(responseBody['results'][0]['types'][3]['values'])
                genesProduct_Taxon = pd.DataFrame.from_dict(responseBody['results'][1]['types'][3]['values'])
                
                annotation_ECOID = pd.DataFrame.from_dict(responseBody['results'][0]['types'][4]['values'])
                genesProduct_ECOID = pd.DataFrame.from_dict(responseBody['results'][1]['types'][4]['values'])
                
                annotation_Aspect = pd.DataFrame.from_dict(responseBody['results'][0]['types'][5]['values'])
                genesProduct_Aspect = pd.DataFrame.from_dict(responseBody['results'][1]['types'][5]['values'])
                
                GOstats =   {'annotations_Hits':annotations_Hits, 'genesProducts_Hits':genesProducts_Hits,
                             'annotation_Refs':annotation_Refs, 'genesProduct_Refs':genesProduct_Refs,
                             'annotation_AssignedBy':annotation_AssignedBy, 'genesProduct_AssignedBy':genesProduct_AssignedBy,
                             'annotation_HitStats':annotation_HitStats,'genesProduct_HitStats':genesProduct_HitStats,
                             'annotation_Taxon':annotation_Taxon, 'genesProduct_Taxon':genesProduct_Taxon,
                             'annotation_ECOID':annotation_ECOID, 'genesProduct_ECOID':genesProduct_ECOID,
                             'annotation_Aspect':annotation_Aspect, 'genesProduct_Aspect':genesProduct_Aspect}
                    
                df_dict[value.replace(':','') + '_' + QG_geneProductType] = GOstats
            
            elif responseBody['numberOfHits'] == 0:
                df_dict[value.replace(':','') + '_' + QG_geneProductType] = {'annotations_Hits':0, 'genesProducts_Hits':0}

    return df_dict                                                                                            


def filterQuickGoDict(QuickGoDict, evidence_code = 'Experimental'):
    '''
    Reads a QuickGo Dictionary of DataFrames and provides preset lists of evidence codes for filtering.
    
    Evidence Codes:
    Experimental: ['GO_ExpEvidence','EXP','IDA','IPI','IMP','IGI','IEP']
    HighThroughput: ['GO_HTPEvidence','HTP','HDA','HMP','HGI','HEP']
    Computational: ['GO_CompEvidence','ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA']
    Author: ['GO_AuthEvidence','TAS','NAS']
    Curated: ['GO_CurEvidence','IC','ND']
    InferredElectronic: ['GO_eEvidence','IEA']
    
    Input:
    QuickGoDict: a QuickGo Dictionary, such as that which can be generated by jweb.fetchQuickGO()
    i.e. 'enables' = GO Molecular Function, 'involved_in' = GO Biological Process, 'part_of' = cellular component
    evidence_code: a string, i.e. 'Curated', which pointing to the evidence codes to filter on. Otherwise a custom list of codes.
        
    Return:
    QckGO_filterDict: A Dictionary of filtered QuickGoD DataFrames (as values) based on keys with a common GO code. 
    
    Limitations:
    Previous versions used an input 'qualifier' to sort by relationship. This has been dropped due to inconsistent use of qualifiers
    by various GO datbases. Sort in pandas if qualifiers are needed.
    '''
    # GO_evidence_presets:
    if evidence_code == 'Experimental':
        evidence = ['EXP','IDA','IPI','IMP','IGI','IEP']
    elif evidence_code == 'HighThroughput':
        evidence = ['HTP','HDA','HMP','HGI','HEP']
    elif evidence_code == 'Computational':
        evidence = ['ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA']
    elif evidence_code == 'Author':
        evidence = ['TAS','NAS']
    elif evidence_code == 'Curated':
        evidence = ['IC','ND']
    elif evidence_code == 'InferredElectronic':
        evidence = ['IEA']
    elif evidence_code == 'All':
        evidence = ['EXP','IDA','IPI','IMP','IGI','IEP',
                    'HTP','HDA','HMP','HGI','HEP',
                    'ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA',
                    'TAS','NAS',
                    'IC','ND',
                    'IEA']
    else:
        evidence = evidence_code
    
    QckGO_filterDict = {}
    for key, value in QuickGoDict.items():
        QckGO_filterDict[key] = value[value['GO_EVIDENCE_CODE'].str.contains('|'.join(evidence))]
    
    return QckGO_filterDict

def getQuickGO_stats(GO_Stats_Dict):
    '''
    Takes a dictionary of quickgo statistical dictionaries such as created by jweb.fetchQuickGO_stats() and tabulates the annotation and
    gene hit statistics.
    
    Input:
    GO_Stats_Dict = a statistical dictionary, can be protein, trembl, swissprot, miRNA, complex.
    
    Returns:
    A Dataframe of record numbers.
    '''
    Stats_Dict_list = []
    Stats_Dict_sp_list = []
    Stats_Dict_trembl_list = []

    for key, value in GO_Stats_Dict.items():
        newdict = {'go ID': key.split('_')[0].replace('GO', 'GO:'), 
        key.split('_')[1] + '_annotations': value['annotations_Hits'],
        key.split('_')[1] + '_genesProducts': value['genesProducts_Hits']}

        if 'Swiss-Prot' in key:
            Stats_Dict_sp_list.append(newdict)
        elif 'TrEMBL' in key:
            Stats_Dict_trembl_list.append(newdict)
        else:
            Stats_Dict_list.append(newdict)

    df_mergeList = []
    meta_list = [Stats_Dict_list, Stats_Dict_sp_list, Stats_Dict_trembl_list]
    for item in meta_list:
        if len(item) > 0:
            df_mergeList.append(item)

    if len(df_mergeList) == 2:
        df_zero = pd.DataFrame.from_dict(df_mergeList[0])
        df_one = pd.DataFrame.from_dict(df_mergeList[1])
        combined_stats_df = pd.merge(df_zero, df_one, on = 'go ID', how = 'outer')
    elif len(df_mergeList) == 3:
        df_zero = pd.DataFrame.from_dict(df_mergeList[0])
        df_one = pd.DataFrame.from_dict(df_mergeList[1])
        df_two = pd.DataFrame.from_dict(df_mergeList[2])
        merge_one = pd.merge(df_zero, df_one, on = 'go ID', how = 'outer')
        combined_stats_df = pd.merge(merge_one, df_two, on = 'go ID', how = 'outer')

    combined_stats_df = combined_stats_df.fillna('not searched')
    
    return combined_stats_df


def mapQuickGO(df_dict, QG_geneProductType = 'protein', gConvertOrganism = 'hsapiens', gConvertTarget = 'ENTREZGENE', writetopath = 'NO'):
    '''
    Similiar to jweb.fetchQuickGO() except instead of accepting an id list, will instead accept a quickgo dataframe (as dicitonary). 
    Useful for tables that have been manually downloaded from the QuickGo site. 
    
    df_dict = As for fetchQuickGO() this should be a dictionary. Assign a similar key name to maintain convention.
    
    Return:
    A dictonary with the searched columns added to each value, a written file if needed.
    
    Limitations:
    The manual TSV downloads from quickGO will have different column string headers to the json fetched ones.
    Also, the column types will vary based on user selection.
    This function will trim the column headers to match those retrieved via gaf request in jweb.fetchQuickGO()
    '''
   
    df_dict_new = {}
    for key, value in df_dict.items():
        cols = list(value.columns)
        value.columns = [i.replace(' ','_') for i in cols]
        df_dict_new[key] = value
    
    df_dict_merge = {}
    for key, value in df_dict_new.items():
        
        value = value[['GENE_PRODUCT_DB', 'GENE_PRODUCT_ID', 'QUALIFIER_QGO', 'GO_TERM', 'REFERENCE', 'ECO_ID', 
                       'WITH/FROM', 'ASSIGNED_BY','ANNOTATION_EXTENSION_QGO','GO_EVIDENCE_CODE']]
        
        gPro_query = value['GENE_PRODUCT_ID'].tolist()
        gPro_dict = mapAnyID_gPro(gPro_query, splitstr = ['nosplit'], geneProductType = QG_geneProductType, gConvertOrganism = gConvertOrganism, gConvertTarget = gConvertTarget)
        gPro_df = gPro_dict['id_map'].rename(columns={'Query': 'GENE_PRODUCT_ID'}).drop_duplicates()
        gPro_df['GENE_PRODUCT_ID'] = gPro_df['GENE_PRODUCT_ID'].astype(str)
        value['GENE_PRODUCT_ID'] = value['GENE_PRODUCT_ID'].astype(str)
        query_map = pd.merge(value, gPro_df, on='GENE_PRODUCT_ID', how='left' )
        df_dict_merge[key] = query_map
    
    cwd  = Path(os.getcwd())
    p = Path(cwd / 'downloads' / 'QuickGo')
    if p.exists() == False:
        Path(cwd / 'downloads' / 'QuickGo').mkdir(parents=True, exist_ok=True)

    if writetopath != 'NO':
        for key, value in df_dict_merge.items():
            fname = key + '_' + gConvertOrganism + '_' + QG_geneProductType + '-' + gConvertTarget + '.csv'
            value.to_csv(cwd / 'downloads' / 'QuickGo' / fname, index = False)
    
    return df_dict_merge, gPro_dict


