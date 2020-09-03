# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 10:27:57 2019

@author: smith.j
"""

from matplotlib import pyplot as plt
import pandas as pd
import math
import seaborn as sns
from scipy.cluster import hierarchy
import jtest
from pylab import cm
from matplotlib import colors
import itertools
import numpy as np
import altair as alt
from scipy import stats
from matplotlib.patches import Rectangle
import jwrangle


CommonPalettesAsHex = {'Set1_qual':['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999'],
                       'Set2_qual':['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'],
                       'Set3_qual':['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f'],
                       'Pastel1_qual':['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec', '#f2f2f2'],
                       'Pastel2_qual':['#b3e2cd', '#fdcdac', '#cbd5e8', '#f4cae4', '#e6f5c9', '#fff2ae', '#f1e2cc', '#cccccc'],
                       'Paired_qual':['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928'],
                       'Dark2_qual':['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666'],
                       'RBP_Class':['#0173b2','#9cd7f7','#a3ed95','#fa564d']}


def getMplPaletteAsHex(qual_palette):
    '''
    Extracts the hex color codes, as a list, from a matplotlib colormap (qualitative map only)
    
    Input:
    qual_palette = list. name of a qualitative matplotlib palette
    
    Return:
    the same list as hex colour codes.
    '''
    
    cmap = cm.get_cmap(qual_palette)    # PiYG
    set_pal = []

    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        set_pal.append(colors.rgb2hex(rgb))
    
    return set_pal

def MQ_showDendrogramQC_mplplot(pGroup_clean, prefix, metadata, title, grid = 'YES', fsize = (8, 8)):
    '''Creates a grid of dendrograms where each plot = MQ group, and its members = samples.
    
    Input: 
    pGroup_clean = proteinGroups file with contaminants removed
    expt_dict = dictionary where keys = groups from MQ group specific parameters, and values = all the sample names that were searched in that group 
    see function ExptDictFromMetadata()
    prefix = the prefix defining the specific column for sample name, i.e: 
    'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
    metadata = metadata
    fsize = figure size as tuple, default 8x8
    
    Limitations: 
    Hierarchy algorithm is hardcoed to 'ward'. see main description. When used for QC the dendrograms are based on proteinGroups
    output and therefore will represent untransformed values, without imputation, and inclusive of values = 0. Thus the dendrogram will also reflect 
    relationships in missing data. The plots are intended for early QC only and should not be used in presentations (z scoring on transformed data 
    without missing values is necessary for proper data representation 
    
    Algorithm: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html'''
    
    expt_dict = jwrangle.MQ_getExptDict(metadata)
    
    if grid == 'NO':
        for key, value in expt_dict.items():
            
            fig, axes = plt.subplots(1,1,figsize=(6,len(value)/2))
            pGroupsTrim, Z = jtest.getDistanceMatrix(pGroup_clean, value, prefix, method = 'ward')
            hierarchy.dendrogram(Z, labels=[i.replace(prefix + ' ','') for i in pGroupsTrim.index], 
                                            leaf_rotation=0, orientation="left", above_threshold_color='grey',
                                            leaf_font_size = 12)
            
            plt.title(key + ': ' + prefix, size = 16)
            plt.tight_layout()
            plt.close
    
    elif grid == 'YES':
        splt_sum = len(metadata['MQgroups'].unique())
               
    #    fig, axes = plt.subplots(math.ceil(splt_sum/2),2, figsize=(8, splt_sum))
        fig, axes = plt.subplots(math.ceil(splt_sum/2),2, figsize=fsize)
        count = 1
        
        for key, value in expt_dict.items():
            pGroupsTrim, Z = jtest.getDistanceMatrix(pGroup_clean, value, prefix, method = 'ward')
            plt.subplot(math.ceil(splt_sum/2), 2, count)
            hierarchy.dendrogram(Z, labels=[i.replace(prefix + ' ','') for i in pGroupsTrim.index], leaf_rotation=0, orientation="left", above_threshold_color='grey')
            count += 1
        
        fig.suptitle(title + prefix, size = 16)
        fig.tight_layout()
        fig.subplots_adjust(top=0.90)
        plt.close


def MQ_getContaminants_sbplot(contaminants, metadata, width = 1, length = 1, layout = 'single'):
    '''Accepts a trimmed dataframe from wrangle.MQ_getContaminants(), comprising the contaminants, protein IDs, gene names and iBAQ scores for each sample.
    Plots elements from this dataframe including one protein or contaminant ID and one Gene name (if available). Plot type is controlled by layout.
    
    Input:
    metadata = metadata
    width = a multiplier to size width of plot (multiplier functions relative to number of labels), default = 1
    length = multiplier to size width of plot (multiplier functions relative to number of labels), default = 1
    layout = 'single' provides one large heatmap for all samples in the MQ run
    'individual' provides individual heatmaps for each MQgroup where contaminant iBAQ > 5
    'grid' = all individual plots, without the cutoff, in a single plot grid
    
    Returns:
    Nothing
    
    Limitations: 
    Hardcoded for iBAQ scores.
    Grid option isn't great'''
    
    expt_dict = jwrangle.MQ_getExptDict(metadata)
    
    ### one cont map
    if layout == 'single':
        plt.figure('Contaminant Heatmap- iBAQ')
        fig, ax = plt.subplots(figsize=(width*len(contaminants)/3,length*len(contaminants)/4))
        ax = sns.heatmap(contaminants,cmap="YlGnBu")
        ax.set_title('Contaminants: iBAQ', size = 16)
        plt.tight_layout()
        plt.close
    
    ## individual maps with ibaq >5 cutoff
    if layout == 'individual':
        for key, value in expt_dict.items():
            df_plot = contaminants.loc[:,value]
            df_plot_trim = df_plot[(df_plot > 5).any(1)]
            wide = len(list(df_plot_trim))*width
            long = len(df_plot_trim)*length
            
            plt.figure(key)
            fig, ax = plt.subplots(figsize=(wide,long/4))
            ax = sns.heatmap(df_plot_trim.loc[:,value],cmap="YlGnBu")
            ax.set_title('Contaminants: ' + key, size = 16)
            plt.tight_layout()
            plt.close

    #grid by individual
    if layout == 'grid':
        wide = len(metadata['MQgroups'].unique())*width
        long = len(contaminants)*length
        fig, axes = plt.subplots(math.ceil(wide/2),2, figsize=(wide*1.5, long))
        count = 1
        
        for key, value in expt_dict.items():
            plt.subplot(math.ceil(wide/2), 2, count)
            sns.heatmap(contaminants.loc[:,value], cmap="YlGnBu")
            count += 1
        
        fig.tight_layout()
        plt.close
    
    plt.close

def MQ_showImputationHistogram_mpl(postimpute, imputed, prefix, sample_list, title):
    '''Creates a grid of histogram overlays where each plot = sample and its histograms = total and imputed dataframe value counts per protein
    
    Input: 
    postimpute = A complete dataframe with imputed values
    Imputed = A dataframe of values that were introduced by imputation
    prefix = the prefix defining the specific column for sample name, i.e: 
    'LFQ intensity', 'peptides', iBAQ', 'Intensity' , 'Unique peptides' etc.
    sample_list = the samples being analysed, typically a value derived from the dictionary generated by ExptDictFromMetadata()
    title = the parent title to all plots in the grid'''
    
    col_list = [prefix + ' ' + i for i in sample_list]
    
    splt_sum = len(sample_list)
    fig, axes = plt.subplots(math.ceil(splt_sum/2),2, figsize=(9, splt_sum*1.5))
    ax0 = plt.subplot(math.ceil(splt_sum/2),2,1)
    count = 1

    for sample in col_list:
        sns.set()
        plt.subplot(math.ceil(splt_sum/2), 2, count, sharey=ax0)
        
        plt.hist(postimpute[sample].dropna(), color='#0173b2', bins=14, range=(10,32))
        plt.hist(imputed[sample].dropna(), color='#de8f05', bins=14, range=(10,32))
        
        handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ['#0173b2','#de8f05']]
        labels= ["Total","Imputed"]
        plt.legend(handles, labels)
        
        plt.title(sample.replace(prefix + ' ',''))
        plt.ylabel('Count')
        plt.xlabel('Log2 LFQ Intensity')

        count += 1
    
#    fig.suptitle(title + prefix, size = 16)
    fig.suptitle(title + prefix, y=1+(1/(4*splt_sum)), size = 16)
    fig.tight_layout()
#    fig.subplots_adjust(top=0.90)
    plt.close

def showPearsonRegression_altair(df, mark_color = '#0173b2'):
    '''
    A ready-to use scatterplot function that will;
    1. Remove all (0,0) data points
    2. Shift all data points with 1 zero to a local minimum near the lowest non-zero data point.
    3. Draw a linear regression between all points with 2 real, non zero x and y values (i.e. data points belonging to #1 and #2 are ignored)
    4. Calculate the pearson correlation for all points in #3 and report it in the plot title.
    
    This function is mainly intended for QC of experimental group replicates. As such it matches well with jwrangle.MQ_getSliceByPrefix()
    
    Returns: Nothing.
    '''                             
    
    # dfrep = df.copy()
    # dfrep.columns = pd.Series(df.columns).replace('.','-').tolist() 

    # Altair won't accept . characters in the coumn header so they must be replaced
    repcol = [i.replace('.','-') for i in list(df.columns)]
    df.columns = repcol
                             
    pairList = list(itertools.combinations(list(df), 2))

    chart = alt.hconcat().configure_mark(color = mark_color, strokeOpacity = 0.6
                                        ).configure_axis(labelFontSize=12, labelFontWeight='normal', titleFontSize = 13, titleFontWeight='normal'
                                        ).configure_title(fontWeight = 'normal', fontSize = 11
                                        )

    for pair in pairList:
        base_df = df[list(pair)].apply(pd.to_numeric)
        base_df = base_df.loc[(base_df != 0).any(1)]
        minLimit = min(base_df.replace(0,np.nan).min())-5
        maxLimit = max(base_df.max())+5

        base_df_replace = base_df.replace(0,minLimit+1.5)

        base = alt.Chart(base_df_replace).mark_point().encode(
        alt.X(pair[0]+':Q', scale=alt.Scale(domain=(minLimit, maxLimit))),
        alt.Y(pair[1]+':Q', scale=alt.Scale(domain=(minLimit, maxLimit))),
        ).properties(
        width=200,
        height=200
        )

        base_df_trim = base_df[(base_df != 0).all(1)]

        pearson = stats.pearsonr(x = base_df_trim[pair[0]], y = base_df_trim[pair[1]])
        pearson_round = [np.round(i, 4) for i in pearson]

        base_r = alt.Chart(base_df_trim).mark_point().encode(
        x=pair[0]+':Q',
        y=pair[1]+':Q',
        ).properties(
        width=200,
        height=200
        ).properties(title= 'pearsonr = ' + str(pearson_round[0]) + '; p-val = ' + str(pearson_round[1]))

        regression = base_r.transform_regression(on = pair[0],regression = pair[1], method = 'linear').mark_line(color = 'black')

        chart |= base + regression

    chart.display()

def BoxPlotByColumn_sbplot(df, prefix, x_lab = '', title = ''):
    '''
    Creates one boxplot per column and aggregates all boxplots into one figure based on the group's common prefix.
    The common prefix is removed from the boxplot label during plot generation.
    
    Input:
    df = the DataFrame to plot from. Often generated by jinspect.MQ_poolDataByCondition() 
    prefix = string. a common string element shared between column names
    x_lab = string. a label for the x-axis. default = ''
    title = string. a title. default = ''
    
    Return:
    No return. Draws a plot.
    '''
    
    cols = [i for i in df.columns if prefix in i]
    sub_df = df[cols]
    
    rename_cols = [i.replace(prefix,'') for i in sub_df]
    sub_df.columns = rename_cols

    plt.figure(title)
    ax = sns.boxplot(data=sub_df, orient="h", palette="Set2", linewidth = 0.5)
    ax.set_xlabel(x_lab)
    ax.set_title(title)

    plt.tight_layout()
    plt.close

def BarPlotByGroup_sbplot(df, x_col, y_col, title, pal = ['#1f78b4'], swarmplot=True, errorbars = 'SEM', 
                          xlabel = False, ylabel = False, xticklabels = False, xrange = None, yrange = None, font_scale = 1.3):
    
    '''Creates a barplot, swarmplot, SEM bar overlay
    
    Input:
    df = the DataFrame to plot from
    x_col = the dataframe column intended for the x axis
    y_col = the dataframe column intended for the y axis
    title = title
    pal = palette, must be submitted as a list or preset (i.e. 'colorblind', 'bright'). Default = blue
    swarmplot: bool. switches overlay data points on/off. Defualt = True
    errorbars: str/bool. switches errorbars on/off and specifies type. must be False, SD, SEM, or CI'. default = SEM
    xlabel: str/bool. Allows x labels to be manually specified
    ylabel: str/bool. Allows x labels to be manually specified
    xticklabels: tuple/bool. user can provide a strings and rotation in that order a tuple i.e. ([string labels], 30)
    Anymore customisation and we are nearing the useful limit of a premade function.
    
    Return:
    No return. Draws a plot.'''
    
    params = dict(data=df, x=df[x_col], y=y_col)
    
    # Preset Cosmetics
    bp_bar_edge = dict(linewidth=1)                                                     #applies to bar border
    bp_error_edge = dict(capsize=0.2, errcolor=".2", edgecolor=".2",errwidth=1.5)       #applies to error bar
    strip_swarm_circles_fill = dict(size=10, edgecolor='black', linewidth=1)            #applies to data points
    
    # sns.barplot, sns.stripplot, sns.pointplot #1    
    
    plt.figure(title)
    sns.set(font_scale=font_scale)
    fig, ax = plt.subplots(figsize=(len(list(set(df[x_col].tolist()))),5))
   
    ax = sns.barplot(palette=pal,                                                 #applies colour order
                            ci=None,                                                    #uses sem as error bar
                            zorder=0,
                            **bp_error_edge,
                            **bp_bar_edge,
                            **params)
    
    if swarmplot == True:
        ax = sns.swarmplot(**strip_swarm_circles_fill,                                      
                          palette=['#ffffff'],
                          zorder=1,                                                         #We must set the z-order. NB: this particular example only works without gridlines.
                          **params)
    
    if errorbars != False:
        if errorbars == 'SEM':
            ax = sns.pointplot(scale=0,                                                         #here, scale determines dot size instead of size. scale = 0, removes the dot
                               palette=['#000000'], 
                               ci=68,**bp_error_edge,                                           #hue shows in edge colour now. set palette to black if not desired.
                               **params)
        elif errorbars == 'SD':
            ax = sns.pointplot(scale=0,                                                         #here, scale determines dot size instead of size. scale = 0, removes the dot
                               palette=['#000000'], 
                               ci='sd',**bp_error_edge,                                           #hue shows in edge colour now. set palette to black if not desired.
                               **params)

        elif errorbars == 'CI':
            ax = sns.pointplot(scale=0,                                                         #here, scale determines dot size instead of size. scale = 0, removes the dot
                               palette=['#000000'], 
                               ci=95,**bp_error_edge,                                           #hue shows in edge colour now. set palette to black if not desired.
                               **params)
        
        else:
            print('wrong errorbar identifier, must be False, SD, SEM, or CI')

    if type(xlabel) == str:
        ax.set_xlabel(xlabel)
    
    if type(ylabel) == str:
        ax.set_ylabel(ylabel)
    
    if xticklabels != False:
        ax.set_xticklabels(xticklabels[0], rotation = xticklabels[1])
    elif xticklabels == False:
        plt.setp(ax.get_xticklabels(), rotation=90)
    
    if xrange != None:
        ax.set_xlim(bottom = xrange[0], top = xrange[1])
    
    if yrange != None:
        ax.set_ylim(bottom = yrange[0], top = yrange[1])
    
    ax.set_title(title, size = 20)
    sns.despine()
    plt.close
    
def LinePlotByGroup_sbplot(df, x_col, y_col, title, pal = ['#1f78b4'], swarmplot=True, errorbars = 'SEM', 
                          xlabel = False, ylabel = False, xticklabels = False):
    
    '''Creates a barplot, swarmplot, SEM bar overlay
    
    Input:
    df = the DataFrame to plot from
    x_col = the dataframe column intended for the x axis
    y_col = the dataframe column intended for the y axis
    title = title
    pal = palette, must be submitted as a list or preset (i.e. 'colorblind', 'bright'). Default = blue
    swarmplot: bool. switches overlay data points on/off. Defualt = True
    errorbars: str/bool. switches errorbars on/off and specifies type. must be False, SD, SEM, or CI'. default = SEM
    xlabel: str/bool. Allows x labels to be manually specified
    ylabel: str/bool. Allows x labels to be manually specified
    xticklabels: tuple/bool. user can provide a strings and rotation in that order a tuple i.e. ([string labels], 30)
    Anymore customisation and we are nearing the useful limit of a premade function.
    
    Return:
    No return. Draws a plot.'''
    
    params = dict(data=df, x=df[x_col], y=y_col)
    
    # Preset Cosmetics
    strip_swarm_circles_fill = dict(size=10, edgecolor='black', linewidth=1)            #applies to data points
    
    # sns.lineplot, sns.stripplot, sns.pointplot #1    
    
    plt.figure(title)
    
    fig, ax = plt.subplots(figsize=(len(list(set(df[x_col].tolist()))),5))
   
    ax = sns.lineplot(palette=pal,                                                 #applies colour order
                            ci=None,                                                    #uses sem as error bar
                            zorder=0,
                            **params)
    
    if swarmplot == True:
        ax = sns.swarmplot(**strip_swarm_circles_fill,                                      
                          palette=['#ffffff'],
                          zorder=1,                                                         #We must set the z-order. NB: this particular example only works without gridlines.
                          **params)
    
    if errorbars != False:
        if errorbars == 'SEM':
            ax = sns.pointplot(scale=0,                                                         #here, scale determines dot size instead of size. scale = 0, removes the dot
                               palette=['#000000'], 
                               ci=68,
                               **params)
        elif errorbars == 'SD':
            ax = sns.pointplot(scale=0,                                                         #here, scale determines dot size instead of size. scale = 0, removes the dot
                               palette=['#000000'], 
                               ci='sd',
                               **params)

        elif errorbars == 'CI':
            ax = sns.pointplot(scale=0,                                                         #here, scale determines dot size instead of size. scale = 0, removes the dot
                               palette=['#000000'], 
                               ci=95,
                               **params)
        
        else:
            print('wrong errorbar identifier, must be False, SD, SEM, or CI')

    if type(xlabel) == str:
        ax.set_xlabel(xlabel)
    
    if type(ylabel) == str:
        ax.set_ylabel(ylabel)
    
    if xticklabels != False:
        ax.set_xticklabels(xticklabels[0], rotation = xticklabels[1])
    elif xticklabels == False:
        plt.setp(ax.get_xticklabels(), rotation=90)
    
    ax.set_title(title, size = 16)
    sns.despine()
    plt.close
    
def ViolinStripPlot_sb(df, x, y, hue, title, hueNumber = 2, yrange = None):
    
    plt.figure()
    sns.set_style('whitegrid')
    ax = sns.violinplot(x=x, y=y, hue = hue, data=df, jitter = True,
                  dodge = True, palette = CommonPalettesAsHex['Set2_qual'], edgecolor = 'black', linewidth = 0.2)
    ax = sns.stripplot(x=x, y=y, hue = hue, data=df, jitter = True, 
                  dodge = True, palette = ['#ffffff', '#ffffff'], edgecolor = 'black', linewidth = 0.2)
    ax.set_title(title)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:hueNumber], labels[:hueNumber])
    
    if yrange != None:
        ax.set_ylim(bottom = yrange[0], top = yrange[1])

def ViolinCompare_sbplot(longForm, title = 'none', ylabel = 'Log2(Intensity)', hueNumber = 2, yrange = None, palette = ['yellow','aqua']):
    '''
    Will create a comparative violin plot with IQR markings. Useful for comparing normalisation effects during QC.
    
    Input:
    longForm = a longform dataframe such as one created by jwrangle.MQ_poolMulti()
    title
    ylabel
    
    Output
    As per description
    '''
    
    plt.figure(figsize = (15,6))
    # sns.set(style="darkgrid")

    ax = sns.violinplot(x="Sample", y="Measure", hue="Value", 
                         data=longForm, palette=palette, split=True,
                         scale="count", inner="quartile",
                         scale_hue=False, bw=.2, linewidth = 1)

    ax.set_title(title, fontsize = 16)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:hueNumber], labels[:hueNumber])
    
    plt.xticks(rotation=0, fontsize = 12)
    plt.ylabel(ylabel)

    if yrange != None:
        ax.set_ylim(bottom = yrange[0], top = yrange[1])


    plt.show()
    plt.close()


