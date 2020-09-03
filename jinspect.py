# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 12:26:37 2019

@author: smith.j
"""
import pandas as pd
import numpy as np
import padua
import jwrangle


def getIsoelectricPoints(protein_list, pI_DataFrame, Avg_pI_col_name = 'Avg_pI'):
    '''
    Scans a precalculated pI database and returns a dataframe of pI results derived from a number of algorithms. Typically I would
    use the 'Avg_pI' column whose average encompasses all but the most bizarre algorithms (see url). Database from 
    http://isoelectricpointdb.org/ and edited to include column edit: Uniprot ID. This added column is the header.split('|')[1] protein label.
    
    Input:
    A list of proteins- a Majority protein IDs list must be split up first i.e. with jinspect.SplitList()
    A pI dataframe- recommend an edited version of that derived from isoelectricpointdb.org.
    
    Returns:
    A DataFrame (maj_prot_pI_df) of pIs where every row represents an entry in the protein list or Majority protein ID list. In the case of the latter
    len(maj_prot_pI_df) may exceed that of the input dataframe.
    
    Limitations:
    Hard-coded to rename 'Avg_pI' from the isoelectricpointdb dataframe
    '''
    unique_protein_list = list(set(protein_list))
        
    maj_prot_pI = []
    for entry in unique_protein_list:
        search = pI_DataFrame[pI_DataFrame['edit: Uniprot ID'] == entry]
        if type(search) == pd.core.frame.DataFrame:
            maj_prot_pI.append(search)

    maj_prot_pI_df = pd.concat(maj_prot_pI)
    maj_prot_pI_df = maj_prot_pI_df.rename(columns={'Avg_pI': Avg_pI_col_name})
    
    return maj_prot_pI_df

def getIsoelectricPointsByClassDict(ClassDict, dfisoelectric_point_database):
    '''
    A version of getIsoelectricPoints that accepts dictionaries created by jtest.MQ_applyClassifyRBP() and returns a dictionary of lists of pIs.
    '''
    ClassDict_pI = {}
    for key, value in ClassDict.items():
        if 'annStatus' not in key and len(value) > 0:
            search_list = jwrangle.SplitList(value['Majority protein IDs'].tolist(), [';','-'], replace = '$&', remove_nonstring_values = True, drop_empty_strings = True)
            pI_df = getIsoelectricPoints(search_list, dfisoelectric_point_database)
            ClassDict_pI[key + '_pI'] = pI_df['Avg_pI'].tolist()
    
    II_list = []
    for key, value in ClassDict_pI.items():
        if 'II' in key and len(value) > 0:
            II_list = II_list + value
            
    ClassDict_pI['II_pI'] = II_list
    
    return ClassDict_pI

def getSumByAnnotation(group_sum_dict, comparison_dictionary, prefix, transform = 'Log2', matchTO = 'Gene names', none_col = 'none col'):

    '''  
    Takes a dictionary of the format produced by jinspect.MQ_getSumByGroup() and sums all values that match each group in a comparison_dictionary. Will produce a normal 
    DataFrame (for browsing) and a melted DataFrame for plotting. Each DataFrame has columns representing the original values summed and those same columns as log transformations. 
    
    Will provide the log2tranformed sum intensity of all proteins belonging to the categories defined by a comparison dictionary.
    
    Input:
    group_sum_dict = the output of function analysis.SumMeasureByGroup() wherein the sum intensity of each gene identified in an MQ search has been calculated for each condition group.
    comparison_dictionary: A dictionary of categories with gene names to compare against where keys = category name (i.e. GO:0003723) and values are a list of ENTREZ gene names.
    prefix: the measure that was tested, i.e. "LFQ Intensity"
    transform: the transformation to be applied to the summed intensities, options 'Log10' or 'Log2'
    matchTO: row identifier to be used in conjunction with the comparison dictionary. i.e. if comparison dictionary is comprised of ENTREZ gene names then the default MQ 'Gene names' should also be selected.
    none_col: Label to apply to members that fall in no category
    
    Returns: A dictionary:
    log_melt_df = a DataFrame containing all groups and their log2(sum Intensity) by comparison_dictionary[key] in long format
    group_SumByKey = a dictionary of DataFrames named by group with the sum intensity and log sum intensity for each comparison_dictionary[key].
    group_SumByGene = a dictionary of DataFrames named by group with the sum intensity and log sum intensity for each gene of comparison_dictionary[key]
    
    Dependent on jwrangle.AnnotateDataFrameCtrls()
    '''
    group_ann_df = {}
    group_ann_sub_df = {}
    for key, value in group_sum_dict.items():
        value_annotated = jwrangle.AnnotateDataFrameCtrls(value, comparison_dictionary, search_match = matchTO, dict_match = matchTO, none_col = none_col)
        group_ann_df[key] = value_annotated['ann_df']
        group_ann_sub_df[key] = value_annotated['ann_subdf']
    
    group_SumByGene = {} #return
    for key, value in group_ann_sub_df.items():
        
        ctrl_sum_df = {}
        for ctrl, df in value.items():
            df.set_index(matchTO, inplace = True)
            df_trim = df.loc[:,[i for i in list(df) if prefix + ' ' in i]]
            df_trim.loc[ctrl]= df_trim.sum(axis = 0)
            ctrl_sum_df[ctrl] = df_trim
    
        group_SumByGene[key] = ctrl_sum_df
    
    # for each df in group_sum_df we transpose the sum data and label the row with ctrl key
    group_sum_dfplot = {} #return
    for key, value in group_SumByGene.items():
        
        ctrl_sum_dfplot = []
        for ctrl, df in value.items():
            sumonly = df.loc[ctrl, :]
            ctrl_sum_dfplot.append(sumonly)
            
        group_sum_dfplot[key] = pd.DataFrame(ctrl_sum_dfplot)
    
    # Log10 transform values
    log_group_sum_dfplot = {}
    log_group_sum_dfplot_unlabeled = {}

    if transform == 'Log10':
        for key, value in group_sum_dfplot.items():
            rename_cols = [transform + '(sum ' + i + ')' for i in list(value.columns)]
            log_group_sum_dfplot_unlabeled[key] = value.apply(np.log10)
            log_group_sum_dfplot[key] = value.apply(np.log2)
            log_group_sum_dfplot[key].columns = rename_cols
    
    elif transform == 'Log2':
        for key, value in group_sum_dfplot.items():
            rename_cols = [transform + '(sum ' + i + ')' for i in list(value.columns)]
            log_group_sum_dfplot_unlabeled[key] = value.apply(np.log2)
            log_group_sum_dfplot[key] = value.apply(np.log2)
            log_group_sum_dfplot[key].columns = rename_cols

    group_SumByKey = {}
    for key, value in group_sum_dfplot.items():
        group_SumByKey[key] = pd.merge(log_group_sum_dfplot[key], value, right_index=True, left_index=True).reset_index().replace(-1*(np.inf), np.nan)
        group_SumByKey[key] = group_SumByKey[key].replace(0,np.nan)
        
    log_melt_dflist = []
    for key, value in log_group_sum_dfplot_unlabeled.items():
        log_melt = log_group_sum_dfplot_unlabeled[key].reset_index()
        rename_cols = [i.replace(prefix + ' ','') for i in list(log_melt)]
        log_melt.columns = rename_cols
        log_melt = pd.melt(log_melt, id_vars = ['index']) #return
        log_melt.rename(columns = {'index':'Category', 'variable':'Condition', 'value':transform + '(sum ' + prefix + ')'}, inplace=True)
        log_melt['Group'] = key
        log_melt_dflist.append(log_melt)
    log_melt_df = pd.concat(log_melt_dflist).replace(-1*(np.inf), np.nan)
    
    return {'log_melt_df':log_melt_df, 'group_SumByKey':group_SumByKey, 'group_SumByGene':group_SumByGene}

def MQ_dropDuplicateIDs(pGroups_map, metadata, prefix = 'Intensity', ID = 'ENTREZGENE_gPro primary', pool = 'measure', drop_ID = 'None', keep_PoolCalcs = False, applyLFQ_filter = ['Intensity', 'iBAQ']):
    '''
    Recommended: dropping dupes and unmapped genes may lead to loss of contaminant rows.
    
    Takes a pGroups DataFrame and
    1. calculates the frequency and intensity mean of 'pool' (which should typically represent all samples in the MQ)
    2. identifies duplicate rows based on 'ID' then selects one row to keep based on pool frequency and pool mean prefix (in that order of priority)
    3. discards remaining duplicates
    4. drops all ID rows labeled with 'drop_ID'
    5. (optional, toggle off with a non-list string or variable) will filter all measures listed in LFQ filter such that only those with an LFQ value will be kept. This allows us to use MQ to apply unique peptide limits to Intensity and iBAQ scores.
    6. Will return a Dictionary; df_keep is the filtered dataframe. df_droprows is a dataframe of deleted rows. df_ droprows dropped intensity values (you have to get those back yourself)
    
    Input:
    pGroups_map = dataframe. A proteinGroups file, most typically with a remapped gene column
    metadata = dataframe.
    prefix = string. A measure, i.e. 'Intensity', used in steps 1-3 for second priority sorting.
    ID = string. The column to search for duplicates, i.e. 'Gene names' will have many dupes, 'Majority protein IDs' will have very few, if any.
    pool = string. The metadata column whose string applies to all samples to be pooled. default = 'measure', which is shared by every sample in the MQ search.
    drop_ID = string/other. An ID string that should entirely removed from the dataframe. default = 'None' (used by jweb.mapAnyID_gPro() when no gene can be found). switch off with any non-string value.
    keep_PoolCalcs = boolean. If True, the results of step 1 will be added to the keep_df.
    applyLFQ_filter = list/other. A list of prefix measures, i.e. Intensity, iBAQ, that should be discarded wherever a companion LFQ value cannot be found. switch off with any non-list value.
    
    Returns: A Dictionary as described above.
    keep_DF = DataFrame
    df_droprows = DataFrame
    
    Note:
    Because the group mean name will often be set on a metadata value shared by all samples, i.e. 'measure' it can often seem a little nonsensical. i.e. if the
    common shared value is 'Intensity' and 'Peptide' is chosen as the prefix then you'll end up with a column 'Peptide mean: Intensity'. If this bothes you then
    add another metadata column, otherwise, accept it as an otherwise meaningless quirk.
    '''
    # annotate pgroups with the condition frequencies and means
    pGroups_freq = MQ_getFrequencyByGroup(pGroups_map, metadata, prefix, group = pool)
    pGroups_freq_means = MQ_getMeasuredMeansByGroup(pGroups_freq, metadata, prefix, group = pool)
    
    # Select all duplicate rows based on multiple column names in list
    duplicateRowsDF = pGroups_freq_means[pGroups_freq_means.duplicated([ID], keep = False)]
    
    # get unique ids
    unique_list = list(duplicateRowsDF[ID].unique())
    
    # get freq and means colnames   
    sort_list = []
    for k in list(metadata[pool].unique()):
        freqstr = prefix + ' Freq: ' + k
        meanstr = prefix + ' mean: ' + k
        sort_list.append(freqstr)
        sort_list.append(meanstr)
    
    sort_list_two = sorted(sort_list)
    
    # iterate on subDFs and sort highest to lowest freq then intensity
    keep_index = []
    drop_index = []
    for i in unique_list:
        subDF = duplicateRowsDF[duplicateRowsDF[ID] == i]
        subDF_sort = subDF.sort_values(by = sort_list_two, ascending = False)
        keep_index.append(int(subDF_sort.iloc[0:1].index.values))
        for i in subDF_sort.iloc[1:].index.values:
            drop_index.append(i)
    
    if keep_PoolCalcs == False:
        df_keep = pGroups_map.loc[~pGroups_map.index.isin(drop_index)]
        df_droprows = pGroups_map.loc[pGroups_map.index.isin(drop_index)]
        
    elif keep_PoolCalcs == True: 
        df_keep = pGroups_freq_means.loc[~pGroups_freq_means.index.isin(drop_index)]
        df_droprows = pGroups_freq_means.loc[pGroups_freq_means.index.isin(drop_index)]
    
    if type(drop_ID) == str:
        df_drop_ID = df_keep[df_keep[ID] == drop_ID]
        df_keep = df_keep[~(df_keep[ID] == drop_ID)]
        df_droprows = pd.concat([df_droprows, df_drop_ID])
    
    if type(applyLFQ_filter) == list:
        # derive the LFQ and int columns, sort both so order matches
        lfq_cols = sorted([i for i in df_keep if 'LFQ intensity' in i])
        
        # extract the subdf for each of the measures
        lfq_df = df_keep[lfq_cols]
        
        # create a boolean filter where an LFQ intensity = True
        lfq_filter = lfq_df > 1
        
        LFQ_filter_dict = {}
        for measure in applyLFQ_filter:
        # Filter Intensity values
            int_cols = sorted([i.replace('LFQ intensity', measure) for i in lfq_cols])
            int_df = df_keep[int_cols]
    
            replacementNA = 'no'
            if 0 in int_df.values:
                replacementNA = 'yes'
    
            # rename the columns
            int_df.columns = lfq_cols
    
            # apply boolean filter & reinstate original column names
            int_df_filter = int_df[lfq_filter]
            int_df_filter.columns = int_cols
            
            # return missing values to previous value
            if replacementNA == 'yes':
                int_df_filter = int_df_filter.fillna(0)
                
            LFQ_filter_dict[measure] = int_df_filter
        
        # just drop columns and merge back    
        for key, value in LFQ_filter_dict.items():
            df_dropcols = df_keep.drop(columns = list(value.columns))
            df_keep = pd.merge(df_dropcols, value, how='outer', left_index=True, right_index=True)

    return {'df_keep': df_keep, 'df_droprows': df_droprows}

def MQ_filterValidValues(pGroup, prefix, group_filter, sample_list, rship = 'equal or more', impute = False, keep_cols = ['Majority protein IDs', 'Gene names', 'Protein names']):
    '''Filters rows in a pGroups dataframe down to those with > or < x for each specified condition. Imputation from normal distribution where desired.
    NB: where 0 or more is required set the minimum valid as -1.
    
    Input: 
    pGroup = proteinGroups dataframe that includes a count column generated for the condition groups of interest, i.e. see ConditionCount()
    prefix = the prefix defining the specific column for sample name, i.e: 
    'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
    group_filter = a dictionary where keys = condition and value = minimum valid value count inclusive
    sample_list = a list of the samples belonging to the group. Can be derived from keys belonging to the function ExptDictFromMetadata().
    impute = True or False, will impute missing values from a normal distribution
    rship = string. 'more than' or 'less than' to define the filter type. default = 'equal or more'
    keep_cols = a list of columns to keep in the filtered dataframe. Default = ['Majority protein IDs', 'Gene names', 'Protein names']
    
    Returns: 
    A dataframe filtered by row, and trimmed for columns of interest including 'Majority protein IDs', 'Gene names', 'Protein names'.
    '''
    
    count_col = prefix + ' Freq: '
    
    rename_group_filter = {}
    for key, value in group_filter.items():
        rename_group_filter[count_col + key] = value
        
    df = pGroup.copy()
    
    # filter from dictionary with multiple conditions
    # code from https://stackoverflow.com/questions/51618985/how-to-filter-a-pandas-dataframe-using-multiple-conditions
    # another option https://stackoverflow.com/questions/39231176/python-dictionary-iteration-for-data-frame-filtering-with-multiple-conditions
    if rship == 'equal or more':    
        df = df.loc[np.logical_and.reduce(list(map(lambda x: df[x]>=rename_group_filter[x], rename_group_filter.keys())))]
    elif rship == 'equal or less':
        df = df.loc[np.logical_and.reduce(list(map(lambda x: df[x]<=rename_group_filter[x], rename_group_filter.keys())))]

    col_list = [prefix + ' ' + i for i in sample_list]
    dfnames = df.loc[:,keep_cols]
    
    if impute == True:
        df2 = df.loc[:,col_list].replace(0,np.nan)
        dfimputed_values, dfimputed_bool = padua.imputation.gaussian(df2, width=0.3, downshift=-1.8)   ## Extract this padua function so that installing the extra package can be avoided.
        df3 = dfnames.merge(dfimputed_values, how='outer', left_index=True, right_index=True)
    
    else:
        df3 = dfnames.merge(df.loc[:,col_list], how='outer', left_index=True, right_index=True)
        df3.replace(0,np.nan, inplace = True)
    
    return df3

def MQ_getContaminants(pGroups, metadata, gene_col = 'Gene names'):
    '''Creates a trimmed dataframe from proteinGroups, comprising the contaminants, protein IDs, gene names and iBAQ scores for each sample.
    Plots elements from this dataframe including one protein or contaminant ID and one Gene name (if available). Plot type is controlled by layout.
    
    Input:
    pGroups = proteinGroups file prior to log2 transformation
    metadata = metadata
    gene_col = a string name for the column from which gene IDs will be extracted. Default = 'Gene names'
    
    Returns:
    The trimmed log transformed conatminant dataframe for all samples
    Limitations: 
    Hardcoded for iBAQ scores.
    '''
    expt_dict = jwrangle.MQ_getExptDict(metadata)
    
    prefix_cols = [i for i in list(pGroups.columns.values) if 'iBAQ ' in i]
    pGroups_contaminants = pGroups[pGroups['Potential contaminant']=='+']
    pGroups_contaminants = pGroups_contaminants.loc[:,['Majority protein IDs', gene_col] + prefix_cols]
    rename_cols = [i.replace('iBAQ ','') for i in list(pGroups_contaminants.columns.values) if 'iBAQ ' in i]
    pGroups_contaminants.columns = ['Majority protein IDs', gene_col] + rename_cols
    
    mpID = [i.replace(';',':') for i in pGroups_contaminants['Majority protein IDs'].astype(str).tolist()]
    mpID_1 = []
    for i in mpID:
        if len(i.split(':')) > 1:
            mpID_1.append(i.split(':')[1])
        else:
            mpID_1.append(i.split(':')[0])
       
    gN = [i.replace(';',':') for i in pGroups_contaminants[gene_col].astype(str).tolist()]
    gN_1 = []
    for i in gN:
        if len(i.split(':')) > 1:
            gN_1.append(i.split(':')[1])
        else:
            gN_1.append(i.split(':')[0])
    
    mpID_gN = []
    for index, value in enumerate(mpID_1):
        mpID_gN.append(value + ': ' + gN_1[index])
            
    pGroups_contaminants['Protein ID: Gene'] = mpID_gN
    pGroups_contaminants = pGroups_contaminants.drop(['Majority protein IDs', gene_col], axis=1)
    pGroups_contaminants = pGroups_contaminants.set_index('Protein ID: Gene')
    
    # Log2 transform IBAQ data
    pGroups_contaminants_log = (np.log2(pGroups_contaminants)).replace(-np.inf, 0)

    return pGroups_contaminants_log

##########
###########
#### CHANGE the name (or change function), too close to MQ_getFrequencybysample
def MQ_getFrequencyByGroup(pGroups, metadata, prefix, group = 'condition'):
    '''Adds a column called "prefix Freq.: condition" to the input dataframe that tallies the number of non np.nan values (or numbers !=0)
    recorded for each condition. This column can be used to filter dataframes by valid values per group.
    
    Input: 
    pGroups = proteinGroups dataframe
    metadata = metadata
    prefix = the common substring defining the specific column for sample name, i.e: 
    'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
    group = the group from metadata by which sample frequencies will be summed, default = 'condition'
    
    Return: 
    A modified pGroups with cumulative frequency columns for each gene by group.
    
    Limitations: 
    Hardcoded to the 'sample' label in metadata'''

    # Adjust prefix to select for prefix + sample instances
    prefix = prefix + ' '
    
    df = pGroups.replace(0, np.nan)
    condition_list = list(metadata[group].unique())
    
    freq_cols = []
    for i in condition_list:
        s_list = [prefix + i for i in metadata.loc[metadata[group]==i, 'sample'].tolist()]
        df[prefix + 'Freq: ' + i] = df[s_list].count(axis='columns')
        freq_cols.append(prefix + 'Freq: ' + i)

    df_freq = df[freq_cols]

    df_concat = pd.merge(pGroups, df_freq, right_index=True, left_index=True, how='outer')

    return df_concat

def MQ_getFrequencyBySample(pmap, metadata_sorted, freqList = ['Intensity', 'iBAQ', 'LFQ intensity'], measure = False):
    '''
    Extracts the cumulative frequency of protein IDs for each sample (based on non-zero measures such as iBAQ, Intensity, and LFQ). 
    Cumulative frequency is a a count based on presence/absence.
    
    Input:
    pmap = An MQ proteinGroups dataframe that has been cleaned of contaminants and may(optional) be annotated with boolean values for columns in freqList.
    metadata_sorted = A metadata dataframe sorted to reflect the order desired in the barplot.
    freqList = a list. Can be a substring shared by samples i.e 'Intensity', 'iBAQ', 'LFQ intensity'
    Can be a list of columns of boolean values used as a mask to count sample membership.
    measure = a string or boolean. Use 'False' if freqList already conatins values such as counting 'Intensity', 'iBAQ', etc. 
    Otherwise choose a measure to which boolean masking will be applied for each sample- i.e. 'iBAQ'.
    
    Returns:
    The metadata_sorted dataframe with freqlist count totals added as a column.
    
    Limitations:
    Plotting function depends on BarPlotByGroup() function.
    '''
    
    metadata_sorted = metadata_sorted.reset_index(drop = True) 
 
    sample_list = metadata_sorted['sample'].tolist()
    
    dict_count = {}
    
    for entry in freqList:
        count_list = []
        for sample in sample_list:
            if measure == False:
                colname = entry + ' ' + sample
                pmap[colname] = pmap[colname].replace(np.nan, 0)
                count = pmap[colname].astype(bool).sum(axis=0)
                count_list.append(count)
                
            else:
                colname = measure + ' ' + sample
                mask = pmap[pmap[entry]][colname]
                mask = mask.replace(np.nan, 0)
                count = mask.astype(bool).sum(axis=0)
                count_list.append(count)
                    
        dict_count[entry] = count_list
    
    count_df = pd.DataFrame.from_dict(dict_count)
    metadata_ann = pd.concat([metadata_sorted, count_df], axis=1)

    return metadata_ann

# The below can be achieved by .groupby()
def MQ_getMeans(annotated_metadata, meanList, group = 'condition'):
    '''
    Takes an metadata-like dataframe that contains counts, or measurements, and returns then means of any defined group as a new dataframe.
    
    Input:
    annotated_metadata: A metadata-like Dataframe with experimental data,such as created by jinspect.getSumBySample() or jinspect.getFrequencyBySample()
    meanList: list of strings representing the columns to extract means from
    group: string representing the metadata grouping to limit meanlist grouping to.
    
    Returns:
    A new dataframe of means.
    '''    
    
    groupList = []
    for i in annotated_metadata[group].tolist():
        if i not in groupList:
            groupList.append(i)
    
    meandict = {}
    for i in groupList:
        subdf = annotated_metadata[annotated_metadata[group] == i]
        meanlist = []
        for k in meanList:
            meanlist.append(int(np.mean(subdf[k].tolist())))
        meandict[i] = meanlist
            
    meanDF = pd.DataFrame.from_dict(meandict, orient='index')
    meanDF.columns = meanList

    return meanDF

def MQ_getMeasuredMeansByGroup(pGroups, metadata, prefix, log_base = False, group = 'condition'):
    '''Adds a column called "prefix Freq.: condition" to the input dataframe that tallies the number of non np.nan values (or numbers !=0)
    recorded for each condition. This column can be used to filter dataframes by valid values per group.
    
    Input: 
    pGroups = proteinGroups dataframe
    metadata = metadata
    prefix = the common substring defining the specific column for sample name, i.e: 
    'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
    group = the group from metadata by which sample frequencies will be summed, default = 'condition'
    
    Return: 
    A modified pGroups with cumulative frequency columns for each gene by group.
    
    Limitations: 
    Hardcoded to the 'sample' label in metadata'''

    print('WARNING: jinspect.MQ_getMeasuredMeansByGroup() has not been checked')

#    #Identify and drop prefix columns that are not sample specific
#    prefix_nosample = [i for i in  pGroups if prefix == i]
#    df = pGroups.replace(0, np.nan).drop(prefix_nosample, axis = 1)
    df = pGroups.replace(0, np.nan)

    # find the groups from metadata
    condition_list = list(metadata[group].unique())
    
    # Adjust prefix to select for prefix + sample instances
    sample_prefix = prefix + ' '
    
    if log_base == False:
        means_cols = []
        for i in condition_list:
            s_list = [sample_prefix + i for i in metadata.loc[metadata[group]==i, 'sample'].tolist()]
            df[prefix + ' mean: ' + i] = df[s_list].mean(axis = 1)
            means_cols.append(prefix + ' mean: ' + i)

    elif log_base == 2:
        means_cols = []
        for i in condition_list:
            s_list = [sample_prefix + i for i in metadata.loc[metadata[group]==i, 'sample'].tolist()]            
            df[prefix + ' linear mean: ' + i] = df[s_list].apply(lambda x: log_base**x).mean(axis = 1)
            df[prefix + ' log2 mean: ' + i] = df[prefix + ' linear mean: ' + i].apply(lambda x: np.log2(x)).replace(-1*(np.inf), np.nan)

            means_cols.append(prefix + ' linear mean: ' + i)
            means_cols.append(prefix + ' log2 mean: ' + i)

    elif log_base == 10:
        means_cols = []
        for i in condition_list:
            s_list = [sample_prefix + i for i in metadata.loc[metadata[group]==i, 'sample'].tolist()]            
            df[prefix + ' linear mean: ' + i] = df[s_list].apply(lambda x: log_base**x).mean(axis = 1)
            df[prefix + ' log10 mean: ' + i] = df[prefix + ' linear mean: ' + i].apply(lambda x: np.log10(x)).replace(-1*(np.inf), np.nan)

            means_cols.append(prefix + ' linear mean: ' + i)
            means_cols.append(prefix + ' log10 mean: ' + i)

    df_mean = df[means_cols]

    df_concat = pd.merge(pGroups, df_mean, right_index=True, left_index=True, how='outer')

    return df_concat

#### Deprecated for the above (introduces antilog corrections)
# def MQ_getMeasuredMeansByGroup(pGroups, metadata, prefix, group = 'condition'):
#     '''Adds a column called "prefix Freq.: condition" to the input dataframe that tallies the number of non np.nan values (or numbers !=0)
#     recorded for each condition. This column can be used to filter dataframes by valid values per group.
    
#     Input: 
#     pGroups = proteinGroups dataframe
#     metadata = metadata
#     prefix = the common substring defining the specific column for sample name, i.e: 
#     'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
#     group = the group from metadata by which sample frequencies will be summed, default = 'condition'
    
#     Return: 
#     A modified pGroups with cumulative frequency columns for each gene by group.
    
#     Limitations: 
#     Hardcoded to the 'sample' label in metadata'''

#     # Adjust prefix to select for prefix + sample instances
#     prefix = prefix + ' '
    
#     df = pGroups.replace(0, np.nan)
#     condition_list = list(metadata[group].unique())
    
#     means_cols = []
#     for i in condition_list:
#         s_list = [prefix + i for i in metadata.loc[metadata[group]==i, 'sample'].tolist()]
#         df[prefix + 'mean: ' + i] = df[s_list].mean(axis = 1)
#         means_cols.append(prefix + 'mean: ' + i)

#     df_mean = df[means_cols]

#     df_concat = pd.merge(pGroups, df_mean, right_index=True, left_index=True, how='outer')

#     return df_concat


def MQ_getMeasuredMeans(pgroup, metadata, measure = 'iBAQ', log_base = False, keep_col = ['Majority protein IDs'], output = 'dict_df'):
    '''Plots pairwise comparisons of the iBAQ means between each condition per MQgroup.
    
    Input:
    pGroup = a proteinGroups DataFrame, preferably log2 transformed and without imputation
    metadata
    measure = the measure defining the specific column for sample name, i.e: 'LFQ intensity', iBAQ', 'Intensity'
    log_base = int/float/bool. If the input dataframe contains log tramsformed values they will have to be returned to linear values before a means can be calculated. 
    Where log_base = False, the means method will assume linear values are the input. Otherwise log_base should equal the log_base previously applied for transformation (i.e. commonly = 2)
    output = 'dict_df' (a dictionary of dataframes trimmed by measure where key = condition. OR ann_df (the input dataframe with mean columns added)

    Return:
    Means columns annotated to the original dataframe
    MQgroups_dict = a dictionary of DataFrames by MQGroups, trimmed by measure with genes as index
    '''
    print('WARNING: jinspect.MQ_getMeasuredMeans() has been checked, by manual calculation, once only')

    expt_dict = jwrangle.MQ_getExptDict(metadata)
    cond_dict = jwrangle.MQ_getCondDict(metadata)

    MQgroups_dict = {}
    for key, value in expt_dict.items():
        cols = [measure + ' ' + i for i in value]
        MQgroups_dict[key] = pgroup.loc[:,keep_col + cols]

    if log_base == False:
        all_dict = {}
        for key, value in MQgroups_dict.items():
            for group in cond_dict[key]:
                grplist = [measure + ' ' + i for i in group[1:]]
                MQgroups_dict[key][grplist] = MQgroups_dict[key][grplist].replace(0, np.nan)
                MQgroups_dict[key][measure + ' ' + group[0] + ' mean'] = MQgroups_dict[key][grplist].mean(axis = 1)
                all_dict[measure + ' ' + group[0] + ' mean'] = MQgroups_dict[key][grplist].mean(axis = 1)
    
    elif log_base == 2:
        all_dict = {}
        for key, value in MQgroups_dict.items():
            for group in cond_dict[key]:
                grplist = [measure + ' ' + i for i in group[1:]]
                MQgroups_dict[key][grplist] = MQgroups_dict[key][grplist].replace(0, np.nan)
                
                MQgroups_dict[key][measure + ' ' + group[0] + ' linear mean'] = MQgroups_dict[key][grplist].apply(lambda x: log_base**x).mean(axis = 1)
                MQgroups_dict[key][measure + ' ' + group[0] + ' log2 mean'] = MQgroups_dict[key][measure + ' ' + group[0] + ' linear mean'].apply(lambda x: np.log2(x)).replace(-1*(np.inf), np.nan)

                all_dict[measure + ' ' + group[0] + ' linear mean'] = MQgroups_dict[key][grplist].apply(lambda x: log_base**x).mean(axis = 1)
                all_dict[measure + ' ' + group[0] + ' log2 mean'] = MQgroups_dict[key][measure + ' ' + group[0] + ' linear mean'].apply(lambda x: np.log2(x)).replace(-1*(np.inf), np.nan)


    elif log_base == 10:
        all_dict = {}
        for key, value in MQgroups_dict.items():
            for group in cond_dict[key]:
                grplist = [measure + ' ' + i for i in group[1:]]
                MQgroups_dict[key][grplist] = MQgroups_dict[key][grplist].replace(0, np.nan)
                
                MQgroups_dict[key][measure + ' ' + group[0] + ' linear mean'] = MQgroups_dict[key][grplist].apply(lambda x: log_base**x).mean(axis = 1)
                MQgroups_dict[key][measure + ' ' + group[0] + ' log2 mean'] = MQgroups_dict[key][measure + ' ' + group[0] + ' linear mean'].apply(lambda x: np.log10(x)).replace(-1*(np.inf), np.nan)
 
                all_dict[measure + ' ' + group[0] + ' linear mean'] = MQgroups_dict[key][grplist].apply(lambda x: log_base**x).mean(axis = 1)
                all_dict[measure + ' ' + group[0] + ' log2 mean'] = MQgroups_dict[key][measure + ' ' + group[0] + ' linear mean'].apply(lambda x: np.log10(x)).replace(-1*(np.inf), np.nan)
   
    DF_all = pd.DataFrame.from_dict(all_dict)
    ann_df = pd.merge(pgroup, DF_all, right_index=True, left_index=True, how = 'outer')
    
    if output == 'ann_df':
        means_out = ann_df
    elif output =='dict_df':
        means_out = MQgroups_dict
    
    return means_out



def MQ_getMissedCleavages(evidence, metadata, drop_contaminants = False):
    '''Calculates the missed cleavage rate per sample from the evidence file.
    
    Input:
    evidence = evidence file
    cond_dict = dictionary where keys = groups from MQ group specific parameters, and values = a list where
    list[0] = condition, and list [1:] = all corresponding sample names that were searched in from that condition.
    
    See function CondDictFromMetadata()
    category = 'expt', 'group', sample'. This option defines the grouping of the output plot.
    Returns:
        
    The missed cleavage dataframe annotated for 'expt', 'group', 'sample'. 
    
    Limitations:
    Relies on cond_dict and, if using CondDictFromMetadata(), by extension a metadata file (not utilised in this function)'''
    
    cond_dict = jwrangle.MQ_getCondDict(metadata)

    if drop_contaminants == True:
        evidence = evidence[evidence['Potential contaminant']!='+']
    
    # Tally and calculate % missed cleavages for each sample from the evidence file
    expt_list = list(set(evidence['Experiment']))                                       #get sample names from evidence, we use this to filter the dataframe
       
    mc_dict = {}                                                                        #Build a dictionary of lists for pandas (nb. if you don't create a list, i.e. have one value, then index=[0] must be passed in dataframe creation)
    for expt in expt_list:
        sample_record = evidence[evidence['Experiment']==expt]
        mc_sum = int(sample_record['Missed cleavages'].sum())
        mc_len = int(len(sample_record))
        mc_perc = round((mc_sum/mc_len)*100)                                            #Divide # observed missed cleavages / # peptides detected on a per sample basis
        mc_dict[expt] = [mc_perc]
    
    #### messed up!!! the group and expt assignment ####
    
    missed_cleavages_df=pd.DataFrame(mc_dict)                                                    #Convert dictionary to dataframe
    
    missed_cleavages_df_red=missed_cleavages_df.T                                                              #Transpose
    missed_cleavages_df_red.reset_index(level=0, inplace=True) 
    missed_cleavages_df_red.rename(columns={'index': 'sample', 0:'% Missed Cleavages'}, inplace=True) #rename columns

    rename = []
    mc_expt = []
    for name in missed_cleavages_df_red['sample'].tolist():
        for key, value in cond_dict.items():
            for group in value:
                if name in group[1:]:
                    mc_expt.append(key)
                    rename.append(group[0])
        
    missed_cleavages_df_red['group'] = rename
    missed_cleavages_df_red['expt'] = mc_expt
    
    missed_cleavages_df_red.sort_values(by=['expt', 'group', 'sample'], inplace = True)
    
    return missed_cleavages_df_red


def MQ_getSumByGroup(pgroups, prefix, metadata, MatchON = 'Gene names'):
    '''
    Can't remember what it does....I think it....
    Will take in the expt_dict and cond_dict variables created from our metadata and calculate the sum values for each experimental group for a given measure.
    
    Input:
    pgroup: A MaxQuant Dataframe
    expt_dict: Output from preprocess.ExptDictFromMetadata(metadata)
    cond_dict: Output fom preprocess.CondDictFromMetadata(metadata)
    prefix: the measurement type, represented as a prefix, in the MaxQuant table. i.e. 'Intensity'
    index: Specifies which column of the MQ file to generate an index from. This index is necessary to trim the MQ dataframe down to a matrix that comprises only index 
    and experimental values in order to sum across rows. 
    
    Returns
    A dictionary of dataframes where keys = experimental groups and values = a dataframe with the specified index where the sum are calculated for each condition group.
    '''
    expt_dict = jwrangle.MQ_getExptDict(metadata)
    cond_dict = jwrangle.MQ_getCondDict(metadata)

    int_col = [i for i in list(pgroups) if prefix + ' ' in i]
    groups_filtered_intensity = pgroups.loc[:,int_col +[MatchON]]
    
    MQgroups_df = {}
    for key, value in expt_dict.items():
        cols = [prefix + ' ' + i for i in value]
        MQgroups_df[key] = groups_filtered_intensity.loc[:,[MatchON] + cols]
        MQgroups_df[key] = MQgroups_df[key].set_index(MatchON)
    
    #for key, value in MQgroups_df.items():
    group_sum_dict = {}
    for condition, group in cond_dict.items():
        means_dict = {}
        for ls in group:
            grplist = [prefix + ' ' + i for i in ls[1:]]
            means_dict[prefix + ' ' + ls[0]] = MQgroups_df[condition].loc[:, grplist].sum(axis = 1)
        means_df = pd.DataFrame(means_dict).dropna(how = 'all')
        means_df = means_df.replace(np.nan, 5)
        means_df.reset_index(inplace = True)
        means_df = means_df.rename(columns={'index': MatchON})
        group_sum_dict[condition] = means_df

    return group_sum_dict



def MQ_getSumBySample(pmap, metadata_sorted, freqList = ['Peptides'], measure = False):
    '''
    Calulates the sum a measure for each sample (based on non-zero measures such as iBAQ, Intensity, and LFQ). 
    
    Input:
    pmap = An MQ proteinGroups dataframe that has been cleaned of contaminants and may(optional) be annotated with boolean values for columns in freqList.
    metadata_sorted = A metadata dataframe sorted to reflect the order desired in the barplot.
    freqList = a list. Can be a substring shared by samples i.e 'Intensity', 'iBAQ', 'LFQ intensity'
    Can be a list of columns of boolean values used as a mask to count sample membership.
    measure = a string or boolean. Use 'False' if freqList already conatins values such as counting 'Intensity', 'iBAQ', etc. 
    Otherwise choose a measure to which boolean masking will be applied for each sample- i.e. 'iBAQ'.
    
    Returns:
    The metadata_sorted dataframe with peptide totals added as a column.
    
    Limitations:
    Uses a concatenation method on metadata_sorted so the index must have been previously ordered.
    '''
    
    metadata_sorted = metadata_sorted.reset_index(drop = True) 
    
    sample_list = metadata_sorted['sample'].tolist()
    
    dict_count = {}
    
    for entry in freqList:
        count_list = []
        for sample in sample_list:
            if measure == False:
                colname = entry + ' ' + sample
                pmap[colname] = pmap[colname].replace(np.nan, 0)
                count = pmap[colname].sum(axis=0)
                count_list.append(count)
                
            else:
                colname = measure + ' ' + sample
                mask = pmap[pmap[entry]][colname]
                mask = mask.replace(np.nan, 0)
                count = mask.sum(axis=0)
                count_list.append(count)
                    
        dict_count[entry] = count_list
    
    count_df = pd.DataFrame.from_dict(dict_count)
    metadata_ann = pd.concat([metadata_sorted, count_df], axis=1) # should be a merge? why is it even like this?

    return metadata_ann

