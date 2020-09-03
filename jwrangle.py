# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 15:10:14 2019

@author: smith.j
"""
from pathlib import Path
import copy
import numpy as np
import pandas as pd
from functools import partial, reduce
import os


#### NEED TO make compatible WITH GROUPED LISTS i.e. maj prot id?
def AnnotateDataFrameCtrls(df_search, ctrl_dictionary, search_match = 'Gene names', dict_match = 'ENTREZGENE_gPro primary', none_col = 'none col'):
    '''Adds new columns named for ctrl_dictionary[key] to the input DataFrame. Column values comprise booleans that 
    indicate whether a gene name, by row entry, has been identified in the value (list) associated with condition_dictionary[key].
    
    Input:
    df_search = the dataframe to be searched
    ctrl_dictionary = a dictionary where keys = name of list, and values = list of gene names. See function ImportCtrlLists()
    search_match = the df_search column name for the values being matched, default 'Gene names'
    dict_match = the ctrl_dictionary[DataFrame] column name for the values being matched, default 'ENTREZGENE_gPro primary'
    
    Returns:
    df_ann = the annotated dataframe described above
    sub_dfs = a dictionary of individual dataframes for each of key in ctrl_dictionary plus one for entries that had no observation
    ctrl_colnames = a list of column names in df_ann for each of key in ctrl_dictionary plus one for entries that had no observation '''

    ctrl_dict_match = {}
    for key, value in ctrl_dictionary.items():
        match_list = value[dict_match].tolist()
        match_list_trim = [i for i in match_list if i != 'na']
        ctrl_dict_match[key] = df_search[search_match].isin(match_list_trim)
    
    ctrl_hits_df = pd.DataFrame.from_dict(ctrl_dict_match)
    
    ctrl_nohits_df = ctrl_hits_df[~ctrl_hits_df.any(1)]                                                                       #find rows that are all false
    ctrl_nohits_df[none_col] = True                                                                           #add zero annotation column to annotation columns
    ctrl_colnames = list(ctrl_nohits_df.columns.values)                                                                       #extract colnames for annotation columns to build sub dfs
    
    df_ann = df_search.merge(ctrl_hits_df, how='outer', left_index=True, right_index=True)                                    #merge ctrl columns to volcano matrix
    ctrl_nohits_df_merge = ctrl_nohits_df[[none_col]]                                                         #extract zero annotation matrix, merge to volcano matrix, fill np.nan
    df_ann = df_ann.merge(ctrl_nohits_df_merge, how='outer', left_index=True, right_index=True)
    
    df_ann[none_col] = df_ann[none_col].replace(np.nan, False)  ####this line probably source of error message?
    
    sub_dfs = {}
    for name in ctrl_colnames:
        sub_dfs[name] = df_ann[df_ann[name]==True]
    
    return {'ann_df':df_ann, 'ann_subdf':sub_dfs}

def concatGO_DataFrameDict(GO_DataFrameDict):
    '''
    A very specific function for concatenating a dictionary of GO Dataframes where key = GO########_swissprot/trembl
    
    Input:
    GO_DataFrameDict: A dictionary of dataframes
    
    Returns
    GO_DataFrameDict: A dictionary of dataframes concatenated on shared GO code.
    '''
    
    key_list = []
    for key, value in GO_DataFrameDict.items():
        key_list.append(key.split('_')[0])
    
    concatdict = {}
    for i in key_list:
        concatgroup = []
        for key, value in GO_DataFrameDict.items():
            if i in key:
                concatgroup.append(value)
        concatdict[i] = concatgroup
    
    GO_DataFrameDict_concat = {}
    for key, value in concatdict.items():
        GO_DataFrameDict_concat[key] = pd.concat(value, sort = True).drop_duplicates().reset_index() ## Drop duplicates necessary?
    
    return GO_DataFrameDict_concat

def importMixedFiles(path, dropSuffix = 'NO'):
    '''
    Imports all dataframes or lists of genes in a specified folder and converts these to a dictionary where
    keys = filename, and values = DataFrame or List.
    dropSuffix = string. whether to include the filename suffix in the dictionary key or not.
    
    Input:
    path = a Pathlib style path command i.e. Path(os.getcwd()) / 'downloads' 
    
    Returns:
    A dictionary as described above
    
    Limitations:
    Hardcoded to discriminate .csv DataFrames, tab delimited .txt DataFrames, and space separated .txt lists'''

    files = os.listdir(path)
    
    file_dict = {}
       
    for f in files:
        if '.csv' in f:
            if dropSuffix != 'NO':
                file_dict[Path(f).stem] = pd.read_csv(path / f)
            else:                
                file_dict[f] = pd.read_csv(path / f)
            #### tried to sniff out delimiter with mixed success....
#            with open(path / f, 'r') as csvfile:
#                dialect = csv.Sniffer().sniff(csvfile.read(1024))
#                file_dict[f] = pd.read_csv(path / f, delimiter = dialect.delimiter)
                
        elif '.tsv' in f:
            if dropSuffix != 'NO':
                file_dict[Path(f).stem] = pd.read_csv(path / f, delimiter = '\t')
            else:
                file_dict[f] = pd.read_csv(path / f, delimiter = '\t')
                
        elif '.txt' in f:
            with open(path / f, 'r') as txtfile:
                check = txtfile.read()
                if '\t' in check:
                    if dropSuffix != 'NO':
                        file_dict[Path(f).stem] = pd.read_csv(path / f, delimiter = '\t')
                    else:
                        file_dict[f] = pd.read_csv(path / f, delimiter = '\t')
                
                elif ' ' in check:
                    if dropSuffix != 'NO':
                        file_dict[Path(f).stem] = check.split(' ')
                    else:
                        file_dict[f] = check.split(' ')
                        
                elif '\n' in check:
                    if dropSuffix != 'NO':
                        file_dict[Path(f).stem] = check.split('\n')
                    else:
                        file_dict[f] = check.split('\n')

        else:
            file_dict[f] = 'No Import'

    return file_dict

def ListsToDataFrameSets(dict_list):
    '''
    Takes a dictionary where keys = labels, values = list of elements, and creates a dataframe of booleans corresponding to their
    sets value. Will throw out elements with len <2 or np.nan. Will accept lists with recurring elements. Useful processing function
    prior to upset.plot. Also useful for generating value counts and sets for plotting.
    
    Returns:
    knitdf = The dataframe described above where 'gene' = the merged upon column.
    
    Limitations:
    Presently hardcoded for strings.
    '''
    
    dict_list_clean = {}
    for key, value in dict_list.items():
        clean_list = []
        for entry in value:
            if ';' in entry and len(entry) > 2:
                clean_list.append(entry.split(';')[0])
            elif len(entry) < 2 or entry == np.nan:
                pass
            else:
                clean_list.append(entry)
        dict_list_clean[key] = clean_list
    
    dict_df = {}
    for key, value in dict_list_clean.items():
        unique_list = list(set(value))
        unique_list_bool = len(unique_list)*[True]
        dict_df[key] = pd.DataFrame({'gene':unique_list, key:unique_list_bool})
    
    my_reduce = partial(pd.merge, on='gene', how='outer') 
    knitdf = reduce(my_reduce, dict_df.values())   
    knitdf.reset_index(inplace=True, drop= True)
    knitdf.replace(np.nan, False, inplace = True)

    return knitdf

def Log2_ByPrefix(pGroups, prefix, sub = np.nan):
    '''
    Applies log2 transformation on columns with a given prefix in the DataFrame pGroups.
    
    sub = a value to substitute for the value -1*(np.inf) generated when log is applied to 0.
    '''
    
    log2_df = pGroups.apply(lambda x: np.log2(x) if prefix in x.name else x).replace(-1*(np.inf), sub)
    
    return log2_df

#def Log2_ByColumn(df, col):
#    '''A transform function utlising masking and broadcasting.
#    
#    Input: 
#    df = dataframe, col = column of values
#    
#    Return: 
#    The input dataframe with the specified column log2 transformed.
#    
#    Limitations: 
#    Not sure if this tolerates np.nan'''
#    
#    # Generate boolean series where != 0
#    mask = (df[col] != 0)
#    # Generate data frame of valid values 
#    df_valid = df[mask]
#    # Apply transform to valid values, modify original column
#    df.loc[mask, col] = np.log2(df_valid[col])
#    
#    return df
#
#def Log2_ByPrefix(pGroups, prefix):
#    '''Identifies prefix columns in a dataframe and log2 transforms them
#    
#    Input: 
#    pGroups = protein groups dataframe
#    
#    Return: 
#    pGroups dataframe log2 transformed prefix columns'''
#    
#    # identify prefix columns
#    prefix_cols = []
#    for entry in list(pGroups):
#        if prefix in entry:
#            prefix_cols.append(entry)
#    
#    # Apply log2 transfrom to all columns with given prefix
#    for col in prefix_cols:
#        Log2_ByColumn(pGroups, col)
#    
#    return pGroups

def MQ_getExptDict(metadata):
    '''Creates a dictionary where keys = groups from MQ group specific parameters, and values = all the sample names
    that were searched in that group. Because samples belonging to a specific experiment are generally searched as separate
    Group specific parameters in MaxQuant, this dictionary will typically represent the experiments being tested'''
    
    expt_dict = {}
    MQ_groups = metadata['MQgroups'].unique().tolist()
    
    for i in MQ_groups:
       b = metadata['MQgroups']==i
       c = metadata['sample'][b.values].tolist()
       expt_dict[i] = c
      
    return expt_dict

def MQ_getCondDict(metadata):
    '''Creates a dictionary where keys = groups from MQ group specific parameters, and values = a list where:
        list[0] = condition, and list [1:] = all corresponding sample names that were searched in from that condition.'''
    
    cond_dict = {}
    MQ_groups = metadata['MQgroups'].unique().tolist()
    
    for i in MQ_groups:
       b = metadata['MQgroups']==i
       c = metadata['condition'][b.values].tolist()
       cond_dict[i] = list(set(c))
       
    con_sample_dict = {}
    
    for key, value in cond_dict.items():
        combined_con_s = []
        for elem in value:
            trim = metadata[metadata['condition'] == elem]
            con_s = [elem] + trim['sample'].tolist()
            combined_con_s.append(con_s)
        con_sample_dict[key] = combined_con_s
    
    return con_sample_dict

def MQ_getInstrumentQC_data(QC_path, getMajorityProteinID_lists = 'YES'):
    '''
    General function set up for expansion
    '''  
      
    return_dict = {}
    
    if getMajorityProteinID_lists == 'YES':
        QC_groups_filter = {}
        QC_groups = importMixedFiles(QC_path, dropSuffix = 'no')
        for key in QC_groups:
            if '.' not in key:
                QC_groups_filter[key] = QC_groups[key]
        
        pGroups_dict = {}
        for folder in QC_groups_filter:
            inst_dict = importMixedFiles(QC_path / folder, dropSuffix = 'yes')
            for file, value in inst_dict.items():
                if 'roups' in file:
                    pGroups_dict[folder + '_' + file] = value
                    
        pGroups_TotalProtein_df = {}
        for key, value in pGroups_dict.items():
            pGroups_TotalProtein_df[key] = list(set(SplitList(value['Majority protein IDs'].tolist(), [';', '-'], replace = '$&', remove_nonstring_values = True, drop_empty_strings = True)))
    
        return_dict['MajProtLists'] = pGroups_TotalProtein_df
    
    return return_dict

def MQ_getThreePassFilter(pGroups, custom_exclusion = False):
    '''Input: pGroups = MQ proteinGroups dataframe
    
    Return: 
    MQ proteinGroups dataframe minus 'Only identified by site', 'Reverse', 'Potential contaminant'.
    
    Limitations: 
    Contaminants hardcoded for MQ proteinGroups tables'''
    
    pGroups = pGroups[pGroups['Only identified by site'] != '+']
    pGroups = pGroups[pGroups['Reverse'] != '+']
    pGroups = pGroups[pGroups['Potential contaminant'] != '+']
    pGroups = pGroups[~pGroups['Majority protein IDs'].str.contains('REV__')]
    
    if type(custom_exclusion)==list:
        pGroups = pGroups[~pGroups['Majority protein IDs'].isin(custom_exclusion)]

    return pGroups

def MQ_poolDataByCondition(pGroup_clean, metadata_sorted, prefix_list = ['Intensity', 'Sequence coverage']):
    '''
    Pools results from each experimental condition into dataframes by prefix and then plots the IQR box and whiskers plot for the group.
    
    Input:
    pGroup_clean = An MQ proteinGroups dataframe that has been cleaned of contaminants.
    metadata_sorted = A metadata dataframe sorted to reflect the order desired for visualisation.
    prefix_list = A list of prefixes representing the desired measurements from pGroups. Default = ['Intensity', 'Sequence coverage']
    
    Returns:
    pool_df = A dictionary of dataframes comprised of the values for each group and prefix. Useful for extra statistical functions.
    '''

    conditions_list = []
    for i in metadata_sorted['condition']:
        if i not in conditions_list:
            conditions_list.append(i)
    
    metrics = prefix_list
    sub_dfs = {}
    for condition in conditions_list:
        sub_list = [i for i in list(pGroup_clean) if condition in i]
        for entry in metrics:
            metrics_samples = [i for i in sub_list if entry in i]
            sub_dfs[entry + ': ' + condition] = pGroup_clean[metrics_samples]
    
    groups_pool_dfs ={}
    for key, value in sub_dfs.items():
        pool_series = []
        
        for sample in sub_dfs[key].columns:
            pool_series = pool_series + sub_dfs[key].loc[:,sample].tolist()
        
        groups_pool_dfs[key] = pd.Series(pool_series)
        
    pool_df = pd.DataFrame(groups_pool_dfs)

    return pool_df

def MQ_poolMulti(pGroups, metadata, melt_list, group = 'condition'):
    '''
    Takes all entries in melt_list and stacks all such values in one column, with adjacent columns denoting 
    Sample ID (names that share the same condition variable) and measurement type (from melt list). Thus a 
    long format dataframe in which sample identity and different measures are maintained. Useful where multiple 
    variables need to be pooled into 3 columns. Returns a dictionary where key = group variables, and value = dataframe.
    
    Input:
    pGroups = 
    metadata = sorted in the desired order for plotting
    melt_list = a list of the measures to be stacked i.e. ['Intensity', 'LFQ intensity']
    group = name of the metadata column to split the dictionary up by
    
    Returns:
    A dictionary of dataframes as described above.
    
    Notes: 
    A variant MQ_poolDataByCondition. Can replace it, cannot be replaced by it.
    '''
    
    conditions_list = []
    for i in metadata[group]:
        if i not in conditions_list:
            conditions_list.append(i)
    
    metrics = melt_list
    sub_dfs = {}
    for condition in conditions_list:
        sub_list = [i for i in list(pGroups) if condition in i]
        for entry in metrics:
            metrics_samples = [i for i in sub_list if entry in i]
            melted = pGroups[metrics_samples].melt(value_name = 'Measure', var_name = 'Sample')
            melted['Sample'] = melted['Sample'].apply(lambda x: x.replace(entry + ' ', ''))
            melted['Value'] = len(melted['Sample'])*[entry]
            sub_dfs[entry + ': ' + condition] = melted
    
    cond_dict = {}
    for condition in conditions_list:
        subdfs_list = []
        for key, value in sub_dfs.items():
            if condition in key:
                subdfs_list.append(value)
        cond_df_pool = pd.concat(subdfs_list, axis = 0)
        cond_dict[condition] = cond_df_pool

    return cond_dict

def MQ_getSliceByPrefix(pGroups, metadata, prefix, group = 'condition', add_col = None):
    '''Creates a Dictionary of DataFrames by group.
    
    Input: 
    pGroups = proteinGroups dataframe
    metadata = metadata
    prefix = the common substring defining the specific column for sample name, i.e: 
    'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
    group = the group from metadata by which sample frequencies will be summed, default = 'condition'
    add_col = a list of columns to bring from the input dataframe, i.e. 'Gene names'. default = None
    
    Return: 
    A Dictionary of DataFrames where key = group label and value = the corresponding slice of pGroups.
    
    Limitations: 
    Hardcoded to the 'sample' label in metadata'''

    # Adjust prefix to select for prefix + sample instances
    prefix = prefix + ' '
    
    condition_list = list(metadata[group].unique())
    
    subDF = {}
    for i in condition_list:
        s_list = [prefix + i for i in metadata.loc[metadata[group]==i, 'sample'].tolist()]
        if add_col != None:
            subDF[i] = pGroups[s_list + add_col]
        else:
            subDF[i] = pGroups[s_list]

    return subDF

def MQ_writeMetadata(pGroups, evidence, expt_df, expt_ID, cwd):
    '''
    Writes expt_df as the metadata dataframe then creates a new pGroups and evidence file that is consistent with the metadata, thus
    making any necessary changes to experiment names in the process.
    
    pGroups = original MaxQuant proteinGroups file
    evidence = original MaxQuant evidence file
    expt_df = the metadata dataframe
    expt_ID = a prefix identifier for file naming
    cwd = the current working directory as a Path object
    
    Returns:
    A dictionary whereby
    metadata = metadata dataframe
    pGroups = the metalabeled proteinGroups file
    evidence = the metalabeled evidence file
    
    Output:
    /metadata/metadata.csv
    /MaxQuant/expt_ID + '_proteinGroups_metalabeled.txt'
    /MaxQuant/expt_ID + '_evidence_metalabeled.txt'
    '''
    
    pGroups_copy = copy.deepcopy(pGroups)
    experiment_names = expt_df['experiment'].tolist()
    colnames = list(pGroups.columns.values)
    
    colnames_dict = dict(zip(expt_df['experiment'].tolist(), expt_df['sample'].tolist()))
    
    colnames_rename = []
    for label in colnames:
        label_list = label.split(' ')
    
        # Generate boolean that checks for any entry in label_list being in experiment names
        bool_check = [elem in label_list for elem in experiment_names]
        # Return the index position of True
        bool_index = [elem for elem, x in enumerate(bool_check) if x]
        # Generate a check for boolean True
        result = any(bool_check)

        if result:
            # Find the sample name for the experiment name that was found
            sample = colnames_dict[experiment_names[bool_index[0]]]
            # Replace old name with new one
            newlabel = label.replace(experiment_names[bool_index[0]],sample)
            colnames_rename.append(newlabel)
        else:
            colnames_rename.append(label)

    pGroups_copy.columns = colnames_rename
    
    p = Path(cwd / 'metadata')
    if p.exists() == False:
        Path(cwd / 'metadata').mkdir()
    
    metastr = expt_ID + '_metadata.csv'
    expt_df.to_csv(cwd / 'metadata' / metastr)
    print('metadata.csv created in metadata folder')
    
    pgstr = expt_ID + '_proteinGroups_metalabeled.txt'
    pGroups_copy.to_csv(cwd / 'MaxQuant' / pgstr, sep = '\t', index=False)
    print(expt_ID + '_proteinGroups_metalabeled.txt created in MaxQuant folder')
    
    #### Create a copy of evidence for analysis
    sample_dict = {}
    for pos, value in enumerate(experiment_names):
        sample_dict[value] = expt_df['sample'].tolist()[pos]
        
    evidence['Experiment'] = evidence['Experiment'].apply(lambda name: sample_dict[name] if name in sample_dict else name)
    
    evstr = expt_ID + '_evidence_metalabeled.txt'
    evidence.to_csv(cwd / 'MaxQuant' / evstr, sep = '\t', index=False)
    print(expt_ID + '_evidence_metalabeled.txt created in MaxQuant folder')
    
    return {'metadata':expt_df, 'pGroups':pGroups_copy, 'evidence': evidence}

def SearchString(df, column, stringlist = ['test']):
    '''
    Scans a given dataframe and specified column for instances of each substring in a list and returns a dataframe with the relevant rows.
    
    Input:
    df: A DataFrame
    column: A column label
    stringlist = a list of substrings to search for
    
    Return
    new_df = a single Dataframe comprised of all relevant rows.
    '''
    df_str = df[df[column].apply(type) == str]

    poslist = []
    for i in stringlist:
        for pos, label in enumerate(df_str[column].tolist()):
            if i in  label:
                poslist.append(pos)
    
    rows = []
    for i in poslist:
        rows.append(df_str.iloc[i,:])

    new_df = pd.DataFrame(rows)
    
    return new_df

def SplitList(ls, split_str, replace = '$&', remove_nonstring_values = False, drop_empty_strings = True):
    '''
    Iterates on a list and splits each item on every substring or character listed in split_str.
    Good for lists of ids that use a substring or character as a separator, i.e. Major Protein Groups.
    
    Input:
    ls = list of string items
    split_str = list of substrings or characters to be split on and removed
    replace = a placeholder substring used in the function. Only needs changing if the default is inappropriate. default = '$&'
    remove_nonstring_values = whether to remove non-string values from the final list. default = False.
    drop_empty_strings = whether to drop empty values. default = True.
    
    Returns:
    sub_list = all the items on the list individually split according to split_str
    '''
    ls_rep = []    
    for i in ls:
        if type(i) == str:
            for k in split_str:
                i = i.replace(k, replace)
            
            ls_rep.append(i)
        
        else:
            if remove_nonstring_values == False:
                ls_rep.append(i)
    
    sub_list = []
    for i in ls_rep:
        if type(i) == str:
            splitList = i.split(replace)
            for k in splitList:
                sub_list.append(k)
                
        else:
            if remove_nonstring_values == False:
                sub_list.append(i)

    if drop_empty_strings == True:
        sub_list = [i for i in sub_list if i !='']

    return sub_list


