# Extracts channels and fates
# To be called in as auxiliary library to produce the mtN summary file
import re
import os
import pandas as pd
import numpy as np
import lib_log_extract    # to use same conditions of log_extract.py


# get list of names of BBHs evolving or not evolving in the intermediate binary configuration
def GetNamesBBH(version,setname,fate_type,binary_type):
    path_to_set = f'./v_{version}'

    prog_BBHs =  pd.read_csv(f'{path_to_set}/progenitors/{fate_type}/{setname}_p_{fate_type}.csv')
    names_BBHs = prog_BBHs.name   # list of all names ended as fate_type

    prog_from_binary =  pd.read_csv(f'{path_to_set}/progenitors/{fate_type}_{binary_type}/{setname}_p_{fate_type}_{binary_type}.csv')
    names_BBH_from_binary = prog_from_binary.name   # list of names evolved as binary_type and ended as fate_type

    names_BBH_no_binary = names_BBHs[~names_BBHs.isin(names_BBH_from_binary)].to_numpy() 
    return names_BBHs, names_BBH_from_binary, names_BBH_no_binary

# get list from single type of binary
def GetNamesSingleBBH(version,setname,fate_type,binary_type, label_name):
    names_BBHs, names_BBH_from_binary, names_BBH_no_binary = GetNamesBBH(version,setname,fate_type,binary_type)
    if label_name == f'{fate_type}_from_all_prog':
        names_list = names_BBHs
    elif label_name == f'{fate_type}_from_{binary_type}':
        names_list = names_BBH_from_binary
    elif label_name == f'{fate_type}_not_from_{binary_type}':
        names_list = names_BBH_no_binary
    return names_list

    
# get list of names of MSBH evolving or not towards a BBH fate
def GetNamesMSBH(version,setname,fate_type,binary_type):
    path_to_set = f'./v_{version}'

    prog_MSBH =  pd.read_csv(f'{path_to_set}/progenitors/{binary_type}/{setname}_p_{binary_type}.csv')
    names_MSBHs = prog_MSBH.name   # list of all names of progenitors in the binary_type configuration

    prog_from_binary =  pd.read_csv(f'{path_to_set}/progenitors/{fate_type}_{binary_type}/{setname}_p_{fate_type}_{binary_type}.csv')
    names_BBH_from_binary = prog_from_binary.name   # list of names evolved as binary_type and ended as fate_type

    names_binary_no_BBHs = names_MSBHs[~names_MSBHs.isin(names_BBH_from_binary)].to_numpy()  # MSBHs not ending as BBH
    return names_MSBHs, names_BBH_from_binary, names_binary_no_BBHs

# get list from single type of binary
def GetNamesSingleMSBH(version,setname,fate_type,binary_type, label_name):
    names_MSBHs, names_BBH_from_binary, names_binary_no_BBHs = GetNamesMSBH(version,setname,fate_type,binary_type)
    if label_name == f'{binary_type}_to_all_fates':
        names_list = names_MSBHs
    elif label_name == f'{binary_type}_to_{fate_type}':
        names_list = names_BBH_from_binary
    elif label_name == f'{binary_type}_not_to_{fate_type}':
        names_list = names_binary_no_BBHs
    return names_list


# Create a column repeating the time of first SN 
# In each row of the correspondent binaries
#
def GetArrayWithTimes(Nrows,arrtimes):
    array = np.zeros(Nrows.sum())
    i=0
    for Nrow,arrtime in zip(Nrows,arrtimes):  
        array[i:i+Nrow] = np.repeat(arrtime,Nrow)
        i += Nrow
    return array

# Count how many time a binary satisfies a mass transfer condition
# e.g. how many times undergoes a successful stable RLO
#
def GetNcounts(history,condition,label):
    cond = history.loc[condition]
    cond[label] = cond.groupby(['name'])['event'].transform('count')
    cond = cond.drop_duplicates('name').set_index('name')[[label]]
    return cond

# From a complete history-like dataframe
# Characterize the number of BEevents
#
def GetNhistory(history,BElabels):
    history_singlebinaryline = history.drop_duplicates('name')[['name']].set_index('name')
    condmerge,condCEsurv,condRLOstable,condWR0,condWR1,condSN0,condSN1,condRem0,condRem1 = lib_log_extract.Condition(history)

    CEsurv = GetNcounts(history,condCEsurv,BElabels[0])
    RLOstable = GetNcounts(history,condRLOstable,BElabels[1])

    dflist = [history_singlebinaryline, CEsurv, RLOstable]
    for BElabel in BElabels[2:]:
        if BElabel == 'Nmerge':
            cond = condmerge
        dflist.append( GetNcounts(history,cond,BElabel) )    
    return pd.concat(dflist,axis=1).fillna(0)


# Extract mass transfer info with respect to the first CO formation
# of a list of names e.g. of names_BBHs as of output from
# names_BBHs, names_BBH_from_binary, names_BBH_no_binary = GetNamesBBH(version,setname,fate_type,binary_type)
#
def GetHistoryFirstSN(version,setname,fate_type,binary_type, label_name):
    # get the correct names list
    if label_name in {f'{fate_type}_from_all_prog',f'{fate_type}_from_{binary_type}',f'{fate_type}_not_from_{binary_type}'}:
        names_list = GetNamesSingleBBH(version,setname,fate_type,binary_type, label_name)
    elif label_name in {f'{binary_type}_to_all_fates',f'{binary_type}_to_{fate_type}',f'{binary_type}_not_to_{fate_type}'}:
        names_list = GetNamesSingleMSBH(version,setname,fate_type,binary_type, label_name)

    # read full history of mass transfer of catalogue
    history_full = pd.read_csv(f'./v_{version}/logs/{setname}_history_full.csv')
    history = history_full[history_full.name.isin(names_list)].sort_values(['name','time'])

    # history with only SN events to extract time of first SN
    historySNonly = history.loc[(history.event == f'SN0') | (history.event == f'SN1') ].sort_values(['name','time']) # select only rows with SN
    timefirstSN = historySNonly.groupby('name').apply(lambda x: x['time'].min()).to_numpy() #time of first SN event

    Nrows =  history.groupby('name').count()['time'].to_numpy()          # N of rows for each binary
    history['time_firstSN'] = GetArrayWithTimes(Nrows,timefirstSN)       # add column with time information

    # separate in two sub-dataframes pre or post the first SN
    historyprefirstSN = history.loc[ history.time <= history.time_firstSN ]
    historypostfirstSN = history.loc[ history.time > history.time_firstSN ]
    first_SNline = historySNonly.drop_duplicates(subset=['name'],keep='first',ignore_index=True)[['name','phase_comp','event']].set_index('name') # selects first remnant formation
    historyRemonly = history.loc[(history.event == 'Rem0') | (history.event == 'Rem1') ].sort_values(['name','time']) # select only rows with Rem
    last_Remline = historyRemonly.drop_duplicates(subset=['name'],keep='last',ignore_index=True)[['name','rem_type']].set_index('name') # selects second remnant formation
    historyWR1prefirstSN = historyprefirstSN.loc[historyprefirstSN.event == 'WR1'][['name','event']].sort_values(['name']).set_index('name') # rows where CO  companion is already a WR
    rows_WR1prefirstSN = historyWR1prefirstSN.rename(columns={'event':'WR1_prefirstSN'})    

    # get N of events with respect to first SN
    if label_name in {f'{fate_type}_from_all_prog',f'{fate_type}_from_{binary_type}',f'{fate_type}_not_from_{binary_type}'}:
        N = GetNhistory(history,['NCE','NRLOs'])
    elif label_name in {f'{binary_type}_to_all_fates',f'{binary_type}_to_{fate_type}',f'{binary_type}_not_to_{fate_type}'}:
        N = GetNhistory(history,['NCE','NRLOs','Nmerge'])
    NpreSN = GetNhistory(historyprefirstSN,['NCEpreSN','NRLOpreSN'])
    NpostSN = GetNhistory(historypostfirstSN,['NCEpostSN','NRLOpostSN'])
    NwithSN = pd.concat([N,NpreSN,NpostSN,first_SNline,last_Remline, rows_WR1prefirstSN],axis=1).fillna(0).rename(columns={'event':'firstSN'})
    # Note: 'phase_comp' is the companion phase at the first SN; 'rem_type' is the remnant_type of the second SN       

    return NwithSN

