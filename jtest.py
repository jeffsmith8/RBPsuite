# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 12:37:36 2019

@author: smith.j
"""
import pandas as pd
import jinspect
import jwrangle
import math
import scipy
import statsmodels.stats.multitest as multi
import qvalue
from scipy.cluster import hierarchy

def labelSigStatus(pGroups, pval_col, foldchange_col, foldchange_cutoff = 1.58, pval_cutoff = 0.05):
    '''Annotates a DataFrame with the enrichment status of each row
    
    Input:
    pGroups: A Dataframe that has been statistically tested
    pval_col: A column name (string) referring to the significance value to be 
    foldchange_col: A column name (string) referring to the foldchange value
    foldchange_cutoff: A cut off value that will determine enrichment, depletion or na. Default is 1.58 (eq. Log2(3))
    pval_cutoff: A cut off value that will determine consideration for enrichment, depletion or na. Default is 0.05
    
    Return:
    The original DataFrame with a column annotated for enrichment status'''
    
    richDF = pGroups.copy()
    richlist = []

    for (idx, row) in richDF.iterrows():
        if row[pval_col] <= 0.05 and row[foldchange_col] >= foldchange_cutoff:
            richlist.append('enriched')
        elif row[pval_col] <= 0.05 and row[foldchange_col] <= -1* foldchange_cutoff:
            richlist.append('depleted')
        else:
            richlist.append('na')
    
    richDF[pval_col + ': enrichment'] = richlist
    
    return richDF

def applyIndependentTTest(pGroup, groupA, groupB, prefix, metadata, nan_pol = 'propagate'):
    '''Performs an independent T-Test of values from a dataframe, calulcates 2 independent FDR corrections (BH and q-val), 
    the means per group and the Log2 fold changes. See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
    
    Input:
    pGroup: A dataframe comprising only those values being tested. Columns containing strings are fine.
    groupA: A group (string from condition column in metadata). Typically the ctrl group, i.e. left side of a volcano plot
    groupB: A group (string from condition column in metadata). Typically the expt group, i.e. right side of a volcano plot
    prefix: A string prefix that indicates the measurement type for each sample to be tested. Typically assigned by MaxQuant, i.e. "LFQ intensity"
    metadata: variable name for metadata dataframe.
    
    Returns:
    The input dataframe annotated with the following columns; 't-statistic', 'p-val', 'BH FDR p-val', 'q-val',
    'mean', 'Log2FC', 'significant'
    
    Limitations:
    Independent t-Test is limited to 2 samples only.'''
    
    pGroup_tTest = pGroup.copy()
    
    matrixA = pGroup_tTest[[prefix + ' ' + i for i in metadata[metadata['condition']==groupA]['sample'].tolist()]]
    matrixB = pGroup_tTest[[prefix + ' ' + i for i in metadata[metadata['condition']==groupB]['sample'].tolist()]]
    
    # calculate stats
    pGroup_tTest['t-statistic'], pGroup_tTest['p-val'] = scipy.stats.ttest_ind(matrixA, matrixB, axis=1, nan_policy = nan_pol)

    ## The FDR corrections were failing where p-val = np.nan was present. Thus for each we will extract a sub df to run each correction then merge back onto the master df
    # BH FDR correction
    bhval_df = pGroup_tTest[pGroup_tTest['p-val']>0].copy()
    a, bhval_df['BH FDR p-val'] = multi.fdrcorrection(bhval_df['p-val'], alpha = 0.01)

    # merge bh into original df
    pGroup_tTest = pd.merge(pGroup_tTest, bhval_df.loc[:, 'BH FDR p-val'], how = 'left', left_index=True, right_index=True)

    # Q-Val FDR correction
    qval_df = pGroup_tTest[pGroup_tTest['p-val']>0].copy()
    qval_df['q-val'] = qvalue.estimate(qval_df['p-val'])
    
    # merge qvals into original df
    pGroup_tTest = pd.merge(pGroup_tTest, qval_df.loc[:, 'q-val'], how = 'left', left_index=True, right_index=True)
    
    pGroup_tTest['-Log(BH FDR p-val)'] = pGroup_tTest['BH FDR p-val'].apply(lambda x: -1*math.log10(x) if x > 0 else 0)
    pGroup_tTest['-Log(q-val)'] = pGroup_tTest['q-val'].apply(lambda x:-1*math.log10(x) if x > 0 else 0)

    pGroup_tTest['sig (q-val)'] = ['significant' if i <= 0.05 else 'not significant' for i in pGroup_tTest['q-val'].tolist()]
    pGroup_tTest['sig (BH FDR p-val)'] = ['significant' if i <= 0.05 else 'not significant' for i in pGroup_tTest['BH FDR p-val'].tolist()]

    groupA_mean = matrixA.mean(axis = 1)
    groupB_mean = matrixB.mean(axis = 1)

    # calculate the log 2 fold changes in both directions
    Log2FC_B = groupB_mean - groupA_mean
    Log2FC_A = groupA_mean - groupB_mean

    pGroup_tTest[groupA + ' mean'] = groupA_mean
    pGroup_tTest[groupB + ' mean'] = groupB_mean
    pGroup_tTest['Log2FC ' + groupB] = Log2FC_B
    pGroup_tTest['Log2FC ' + groupA] = Log2FC_A
        
    return pGroup_tTest

def applyIndependentTTest_ImputeCheck(preImpute, postImpute, groupA, groupB, prefix, metadata, nan_pol = 'propagate'):
    '''
    Performs the same independent T-Test and FDR correction approach as jtest.IndependentTTest() but will test both pre and post imputation dataframes.
    By testing both dataframes, two additional columns 'Impute Dependency (BH)' and 'Impute Dependency (q-val)' will be generated. These columns indicate whether the 
    finding of a significant value was dependent on the imputation method used. Useful to colour volcano plots and thence discriminate between real and imputed data.
    
    Input:
    preImpute: A dataframe with missing values.
    postImpute: A dataframe with missing values imputed by any means.
    groupA: A group (string from condition column in metadata). Typically the ctrl group, i.e. left side of a volcano plot
    groupB: A group (string from condition column in metadata). Typically the expt group, i.e. right side of a volcano plot
    prefix: A string prefix that indicates the measurement type for each sample to be tested. Typically assigned by MaxQuant, i.e. "LFQ intensity"
    metadata: variable name for metadata dataframe.

    Returns:
    the postImpute dataframe with standard columns expected of jtest.IndependentTTest() plus 'Impute Dependency (BH)' and 'Impute Dependency (q-val)'.
    '''

    ttest_pre = applyIndependentTTest(preImpute, groupA, groupB, prefix, metadata, nan_pol = 'propagate')
    ttest_post = applyIndependentTTest(postImpute, groupA, groupB, 'LFQ intensity', metadata, nan_pol = 'propagate')
    
    testcompare_qval = pd.merge(ttest_pre[['t-statistic', 'sig (q-val)']], ttest_post[['sig (q-val)']], how='outer', left_index=True, right_index=True)
    
    compare_qval = (ttest_pre['sig (q-val)'] == ttest_post['sig (q-val)'])
    testcompare_qval['compare'] = compare_qval
    
    impute_dependency_qval = []
    for index, value in testcompare_qval['t-statistic'].isnull().iteritems():
        if value == True and testcompare_qval.loc[index, 'compare'] == False:
            impute_dependency_qval.append('dependent')
        elif value == False and testcompare_qval.loc[index, 'compare'] == True:
            impute_dependency_qval.append('n/a')
        elif value == False and testcompare_qval.loc[index, 'compare'] == False:
            impute_dependency_qval.append('change')
        else:
            impute_dependency_qval.append('dependent')
            
    testcompare_BH = pd.merge(ttest_pre[['t-statistic', 'sig (q-val)']], ttest_post[['sig (q-val)']], how='outer', left_index=True, right_index=True)
    
    compare_BH = (ttest_pre['sig (BH FDR p-val)'] == ttest_post['sig (BH FDR p-val)'])
    testcompare_BH['compare'] = compare_BH
    
    impute_dependency_BH = []
    for index, value in testcompare_BH['t-statistic'].isnull().iteritems():
        if value == True and testcompare_BH.loc[index, 'compare'] == False:
            impute_dependency_BH.append('dependent')
        elif value == False and testcompare_BH.loc[index, 'compare'] == True:
            impute_dependency_BH.append('n/a')
        elif value == False and testcompare_BH.loc[index, 'compare'] == False:
            impute_dependency_BH.append('change')
        else:
            impute_dependency_BH.append('dependent')

    ttest_post['Impute Dependency (BH)'] = impute_dependency_BH
    ttest_post['Impute Dependency (q-val)'] = impute_dependency_qval

    return ttest_post

def getDistanceMatrix(pGroups, sample_list, prefix, method = 'ward'):
    '''Calculates a 2D numpy array distance matrix used for dendrogram plotting based on selected dataframe columns. 
        Also removes any rows composed solely of 0s. See:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
        
        Input: 
        pGroups = proteinGroups dataframe
        sample_list = sample names (not column names). This list can be derived from metadata.
        prefix = the prefix defining the specific column for sample name, i.e: 
        'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
        NB: prefix can be left empty if unnecessary
        method: default = 'ward'. See Scipy docs for other options.
        
        Return: 
        A plot object. Also plots when called.
        
        Limitations: 
        Not set to return an object (just plot). Colors hardcoded and method set to 'ward' by default. '''
    
    # Adjust prefix to select for prefix + sample instances
    prefix = prefix + ' '
    
    # Trim dataframe to desired columns and transpose 
    group = [prefix + i for i in sample_list]
    pGroups_trim = pGroups.loc[:,['Gene names'] + group]
    pGroups_trim = pGroups_trim.loc[(pGroups_trim.sum(axis=1) != 0)]           #remove rows sum = 0
    pGroups_trim = pGroups_trim.set_index('Gene names', drop = True)
    pGroups_trim = pGroups_trim.T
    
    # Create distance matrix
    Z = hierarchy.linkage(pGroups_trim, method)
    
    # Plot dendrogram
    return pGroups_trim, Z

def MQ_getFoldChange(pGroup, groupA, groupB, prefix, metadata):
    '''
    Scans a proteinGroup matrix for the specified groups calculates mean and fold change blah blah blah
    
    make it better so only pGroup and metadata are necessary blah blah blah?
    
    pGroup: A dataframe comprising only those values being tested. Columns containing strings are fine.
    groupA: A group (string from condition column in metadata). Typically the ctrl group, i.e. left side of a volcano plot
    groupB: A group (string from condition column in metadata). Typically the expt group, i.e. right side of a volcano plot
    prefix: A string prefix that indicates the measurement type for each sample to be tested. Typically assigned by MaxQuant, i.e. "LFQ intensity"
    metadata: variable name for metadata dataframe.

    '''
    pGroup_fc = pGroup.copy()
    
    matrixA = pGroup_fc[[prefix + ' ' + i for i in metadata[metadata['condition']==groupA]['sample'].tolist()]]
    matrixB = pGroup_fc[[prefix + ' ' + i for i in metadata[metadata['condition']==groupB]['sample'].tolist()]]

    groupA_mean = matrixA.mean(axis = 1)
    groupB_mean = matrixB.mean(axis = 1)

    # calculate the log 2 fold changes in both directions
    Log2FC_B = groupB_mean - groupA_mean
    Log2FC_A = groupA_mean - groupB_mean

    pGroup_fc[groupA + ' mean'] = groupA_mean
    pGroup_fc[groupB + ' mean'] = groupB_mean
    pGroup_fc['Log2FC ' + groupB] = Log2FC_B
    pGroup_fc['Log2FC ' + groupA] = Log2FC_A

    return pGroup_fc

def MQ_applyClassifyRBP(pGroup_ConCount, measure, metadata, MQgroup, condition_NC, condition_PC, valids, add_cols = []):
    '''
    Takes a pGroups df that has been modified by jinspect.MQ_getFrequencyByGroup() and scores all entries for RBP Status based on non-crosslinked and crosslinked measurements.
    
    RBP Classes are mutually exclusive as follows:
    'I' = RBP Purified; proteins exclusively found in crosslinked group AND in #samples => .
    'IIa' = Purified_RBP candidate; proteins exclusively found in crosslinked group AND found in 1 < #samples < valids.
    'IIb' = Purified_RBP candidate; proteins that meet criteria 'I' but appear in exactly 1 non-crosslinked sample.
    'IIc' = QE RBP candidates; proteins identifed by significance (T-test + q-value) and fold change enrichment (> 3x fc). No imputation is used, dependent on MQ LFQ normalisation and intensity assignment.
    'NC' = Not Classified; proteins that do not meet the QE enrichment criteria set by 'IIc'. Includes significantly depleted proteins. No imputation is used, dependent on MQ LFQ normalisation and intensity assignment.
    'ND' = No Data; proteins that do not meet any of the above criteria

    Input:
    pGroup_ExpCount: A pGroups dataframe that has been annotated by preprocess.ConditionCount(). Must contain LFQ values.
    measure: A string representing the intensity unit to be used for determining protein presence. Choose 'iBAQ' or 'Intensity' if unsure.
    metadata: variable name for metadata dataframe.
    MQgroup: A string name for the non-crosslinked vs crosslinked condition group to be tested, found in metadata.
    condition_NC: A string name representing all samples of the non-crosslinked group, found in metadata
    condition_PC: A string name representing all samples of the crosslinked group, found in metadata
    valids: Integer. Minimum number of samples in which a protein must appear to be considered for class 'I'.
    add_cols: A list of columns to be included in the output dataframes. The default columns for jinspect.MQ_filterValidValues are always present.

    groupA: A group (string from condition column in metadata). Typically the ctrl group, i.e. left side of a volcano plot
    groupB: A group (string from condition column in metadata). Typically the expt group, i.e. right side of a volcano plot
    prefix: A string prefix that indicates the measurement type for each sample to be tested. Typically assigned by MaxQuant, i.e. "LFQ intensity"
    metadata: variable name for metadata dataframe.
    
    Returns:
    'RBP_Dict': A dictionary of dataframes for each of the above classes and two master datframes annotated by class.
    A printed message to console that confirms if all proteins in the original datframe were accounted for.
    '''
    
    expt_dict = jwrangle.MQ_getExptDict(metadata)
    
    NC_freq_col =  measure + ' Freq: ' + condition_NC
    PC_freq_col =  measure + ' Freq: ' + condition_PC
    pGroup_ExpCount = pGroup_ConCount[(pGroup_ConCount[NC_freq_col] > 0) | (pGroup_ConCount[PC_freq_col] > 0)]
    
    #### Category I RBP "purified"
    ### Define 'purified' RBPs as those found in crosslinked group AND have no trace in any sample of the non-crosslinked group.
    ## Find crosslinked proteins that meet the min. valids requirement
    pc_df_valid = jinspect.MQ_filterValidValues(pGroup_ExpCount, measure, {condition_NC:0, condition_PC:valids}, expt_dict[MQgroup], impute = False, 
                                                keep_cols = ['Majority protein IDs', 'Gene names', 'Protein names'] + add_cols)

    ## Find all proteins appearing in the non-crosslinked group.
    nc_df = jinspect.MQ_filterValidValues(pGroup_ExpCount, measure, {condition_NC:1, condition_PC:0}, expt_dict[MQgroup], impute = False,
                                          keep_cols = ['Majority protein IDs', 'Gene names', 'Protein names'] + add_cols)

    ## Extract Category Ia purified RBPs as those which pc_df_valid - nc_df
    pure_df = pc_df_valid[~pc_df_valid.index.isin(nc_df.index.tolist())].copy()
    pure_df['RBP Class'] = 'I'
    pure_df['RBP subClass'] = ''

    #### Category II RBP "candidate"
    ### a) Proteins observed 1 < detection < min. valids AND not identified in any non-crosslinked sample
    ## Find crosslinked proteins that do not meet the min. valids requirement
    pc_df_invalid = pGroup_ExpCount[~pGroup_ExpCount.index.isin(pc_df_valid.index.tolist())]

    ## Extract Category II purified RBP ('candidates')
    pc_df_one = pGroup_ExpCount[pGroup_ExpCount[measure + ' Freq: ' + condition_PC] == 1]
    pure_cand_df = pc_df_invalid[~pc_df_invalid.index.isin(nc_df.index.tolist())]
    pure_cand_df_trim = pure_cand_df[~pure_cand_df.index.isin(pc_df_one.index.tolist())].copy()
    pure_cand_df_trim['RBP Class'] = 'II' 
    pure_cand_df_trim['RBP subClass'] = 'a'

    ### b) Proteins which DO meet the min. valids requirement AND appear in 1 non-crosslinked sample
    nc_df_one = pGroup_ExpCount[pGroup_ExpCount[measure + ' Freq: ' + condition_NC] == 1]
    
    ## Find union of pc_df_valid and nc_df_one
    almostpure_cand_df = pc_df_valid[pc_df_valid.index.isin(nc_df_one.index.tolist())].copy()
    almostpure_cand_df['RBP Class'] = 'II'
    almostpure_cand_df['RBP subClass'] = 'b'

    ### c) Proteins which are identified by quantitative enrichment without imputation.
    ## Derive DataFrame for QE
    qe_df = jinspect.MQ_filterValidValues(pGroup_ExpCount, 'LFQ intensity', {condition_NC:2, condition_PC:0}, expt_dict[MQgroup], impute = False,
                                          keep_cols = ['Majority protein IDs', 'Gene names', 'Protein names'] + add_cols)

    ## Apply tTest to all members of nc_df
    tTest_df = applyIndependentTTest(qe_df, condition_NC, condition_PC, 'LFQ intensity', metadata, nan_pol = 'omit')
    
    ## Extract Category II enriched RBP ('candidates') based on q-val significane + FC > 3
    qe_cand_df = tTest_df[(tTest_df['sig (q-val)']=='significant') & (tTest_df['Log2FC ' + condition_PC] >= 1.585)].copy()
    qe_cand_df['RBP Class'] = 'II'
    qe_cand_df['RBP subClass'] = 'c'
    
    ## Extract N.S. (not significant) proteins
    qe_ns_df = tTest_df[tTest_df['Log2FC ' + condition_PC] < 1.585].copy()
    qe_ns_df['RBP Class'] = 'ND'
    qe_ns_df['RBP subClass'] = ''

    ### Category N.D proteins
    ## These proteins do not fit in any of the above categories  
    catI_catII_list = pure_df.index.tolist() + pure_cand_df_trim.index.tolist() + almostpure_cand_df.index.tolist() + qe_cand_df.index.tolist() + qe_ns_df.index.tolist()
    nd_df = pGroup_ExpCount[~pGroup_ExpCount.index.isin(catI_catII_list)]
    nd_df = jinspect.MQ_filterValidValues(nd_df, 'LFQ intensity', {condition_NC:0, condition_PC:0}, expt_dict[MQgroup], impute = False,
                                          keep_cols = ['Majority protein IDs', 'Gene names', 'Protein names'] + add_cols)
    
    ## Calculate fold changes and annotate nd_df with any tTesting results from qe_df
    nd_df = MQ_getFoldChange(nd_df, condition_NC, condition_PC, 'LFQ intensity', metadata)
    cols_to_use = tTest_df.columns.difference(nd_df.columns)
    nd_df_merge = pd.merge(nd_df, tTest_df[cols_to_use], how = 'left', left_index=True, right_index=True).copy()
    nd_df_merge['RBP Class'] = 'NC'
    nd_df_merge['RBP subClass'] = ''

    ## Annotate pGroup_ExpCount with Category column for CatI, CatII, N.D.
    status_df_list = [pure_df.loc[:,['RBP Class','RBP subClass']], pure_cand_df_trim.loc[:,['RBP Class','RBP subClass']],
                      almostpure_cand_df.loc[:,['RBP Class','RBP subClass']], qe_cand_df.loc[:,['RBP Class','RBP subClass']], 
                      qe_ns_df.loc[:,['RBP Class','RBP subClass']], nd_df_merge.loc[:,['RBP Class','RBP subClass']]]
    
    status_df = pd.concat(status_df_list, axis=0, join='outer')
    pGroup_ExpCount_AnnStatus = pd.merge(pGroup_ExpCount, status_df, how = 'left', left_index=True, right_index=True)
    pGroup_ExpCount_AnnStatus['MQgroup'] = MQgroup

    ## Put relevant dfs into dictionary
    RBP_Dict = {'Input_df_annStatus':pGroup_ExpCount_AnnStatus,
                'Summary_df_annStatus': pGroup_ExpCount_AnnStatus.loc[:,add_cols + ['Gene names', 'RBP Class', 'RBP subClass', 'MQgroup']],
                'I': pure_df,
                'IIa': pure_cand_df_trim,
                'IIb': almostpure_cand_df,
                'IIc': qe_cand_df,
                'NC': qe_ns_df,
                'ND': nd_df_merge}

    ## Finally check all proteins have been accounted for
    if len(RBP_Dict['Input_df_annStatus']) == len(pGroup_ExpCount):
        print('All proteins classified = True')
    else:
        print('All proteins classified = False')
        
    return RBP_Dict
