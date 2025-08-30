# Extracts channels and fates
# To be called in as auxiliary library to produce the mt_history_full.csv summary file
import re
import os
import pandas as pd
import numpy as np


# Extract counts from the log file with regular expressions
# For instance, in the logfile, CE lines have the following structure
# B;name;ID;CE;time;ID1:M1:MHe1:MCO1:phase1:rem_type1:ID2:M2:MHe2:MCO2:phase2:rem_type2:a:afin:fate
# Where the labels '1' and '2' indicate respectively the values of the star
# That started (primary) or suffered(secondary) the CE
# Warning: the order of the output regex_list is important and will be then used to sort and classify the events
#
def LogRegex():
    exp = '\d+.\d+[e][+\-]?\d+'
    fl = '\d+.\d+'
    num_sign = '\-?\d'   # single digit with or without the minus sign before
    
    regex_CE= f'B;(\d+);\d+;(CE);({exp});\d+:{exp}:{exp}:{exp}:\d+:\d+:\d+:{exp}:{exp}:{exp}:\d+:\d+:{exp}:{exp}:(\d+)'  # add fateCE
    regex_RLOi = f'B;(\d+);\d+;(RLO_BEGIN);({exp});'
    regex_RLOf = f'B;(\d+);\d+;(RLO_END);({exp});'
    regex_COLL = f'B;(\d+);\d+;(COLLISION);({exp});'
    regex_MERGER = f'B;(\d+);\d+;(MERGER);({exp});'
    regex_SWALL = f'B;(\d+);\d+;(SWALLOWED);({exp});'
    regex_WR0 = f'S;(\d+);0;(HENAKED);({exp});'
    regex_WR1 = f'S;(\d+);1;(HENAKED);({exp});'
    regex_SN0 = f'B;(\d+);\d+;(BSN);({exp});0:{exp}:{exp}:{exp}:\d+:\d+:\d+:{exp}:{exp}:{exp}:(\d+)'                   # add phase of companion
    regex_SN1 = f'B;(\d+);\d+;(BSN);({exp});1:{exp}:{exp}:{exp}:\d+:\d+:\d+:{exp}:{exp}:{exp}:(\d+)'                   # add phase of companion
    regex_remtype0 = f'S;(\d+);0;(SN);({exp});{fl}:{fl}:{fl}:{fl}:({num_sign})'                   # add remtype formed after SN
    regex_remtype1 = f'S;(\d+);1;(SN);({exp});{fl}:{fl}:{fl}:{fl}:({num_sign})'                   # add remtype formed after SN
    return [regex_RLOi,regex_COLL,regex_CE,regex_MERGER,regex_SWALL,regex_RLOf,regex_WR0,regex_WR1,regex_SN0,regex_SN1,regex_remtype0,regex_remtype1]


# Extract single and binary events from logfile
# Events are previously identified through LogRegex
# To avoid possible errors in the workflow (especially for SWALLOWED events)
# We manually create empty dataframes following the same order of LogRegex
# And only for CE events (i==2) we use different column names
#
def LogToEvents(logname,path_to_sevn_output):
    logtext = open(path_to_sevn_output+logname,'r').read()
    regex_list = LogRegex()

    cols = ['name','event','time']
    df_event_list = []
    for i, regex_event in enumerate(regex_list):
        mask_event = re.findall(regex_event,logtext)
        if len(mask_event)!=0:
            if i==2:
                df_events = pd.DataFrame(np.asarray(mask_event), columns = cols + ['fateCE'])
            elif (i==8) or (i==9):
                df_events = pd.DataFrame(np.asarray(mask_event), columns = cols + ['phase_comp'])
            elif (i==10) or (i==11):
                df_events = pd.DataFrame(np.asarray(mask_event), columns = cols + ['rem_type'])
            else:
                df_events = pd.DataFrame(np.asarray(mask_event), columns = cols)
        else:
            if i==2:
                df_events = pd.DataFrame(columns = cols + ['fateCE'])
            elif (i==8) or (i==9):
                df_events = pd.DataFrame(columns = cols + ['phase_comp'])
            elif (i==10) or (i==11):
                df_events = pd.DataFrame(columns = cols + ['rem_type'])         
            else:
                df_events = pd.DataFrame(mask_event, columns = cols)         

        df_events['order'] = i  # add column to allow sorting event timeline for each single binary
        df_event_list.append(df_events)

    # List of dataframes, one dataframe corresponding to one mt event (in order): RLOi, COLL,CE,MERGER,SWALL,RLOf
    mtevents = pd.concat(df_event_list[:6]).astype({'name':'int','event':'str','time':'float','order':'int'}).sort_values(['name','time','order'])
    # List of dataframes, one dataframe corresponding to one singe event (in order): WR0,WR1,SN0,SN1, REM0,REM1
    singleevents = pd.concat(df_event_list[6:]).astype({'name':'int','event':'str','time':'float','order':'int'}).sort_values(['name','time','order'])   
    return mtevents, singleevents


# Classify binary events e.g. merging binaries or the posibility to survive CEs
#
def Condition(dfhistory):
    condmerge = dfhistory.BEvent.isin([17,10,12,15,13])   # correspond to possible SEVN values of 'BEvent'
    condCEsurv = dfhistory.BEvent.isin([11,14])           # correspond to possible SEVN values of 'BEvent'
    condRLOstable = dfhistory.BEvent.isin([5])            # correspond to possible SEVN values of 'BEvent'
    condWR0 = dfhistory.BEvent.isin([30])                 # manually added to follow the same formalism
    condWR1 = dfhistory.BEvent.isin([40])                 # manually added to follow the same formalism
    condSN0 = dfhistory.BEvent.isin([50])                 # manually added to follow the same formalism
    condSN1 = dfhistory.BEvent.isin([60])                 # manually added to follow the same formalism
    condrem0 = dfhistory.BEvent.isin([70])                 # manually added to follow the same formalism
    condrem1 = dfhistory.BEvent.isin([80])                 # manually added to follow the same formalism   
    return condmerge,condCEsurv,condRLOstable,condWR0,condWR1,condSN0,condSN1, condrem0, condrem1


# Create the 'history' dataframe where each binary or single star event is classified
# In particular, we distinguish the possibility that RLOs and CEs allow for survival or merging
# As well as the possibility that RLOs remain stable or become CEs
#
def GetHistory(logname,path_to_sevn_output):
    mtevents, singleevents = LogToEvents(logname,path_to_sevn_output)

    RLOiSWALLOWED = mtevents.loc[(mtevents['event'].eq('RLO_BEGIN')) & (mtevents['event'].shift(-1).eq('SWALLOWED'))][['name','time']]
    RLOiMERGE = mtevents.loc[(mtevents['event'].eq('RLO_BEGIN')) & (mtevents['event'].shift(-1).eq('MERGER'))][['name','time']]
    RLOiCEsurv = mtevents.loc[(mtevents['event'].eq('RLO_BEGIN')) & (mtevents['event'].shift(-1).eq('CE')) & ((mtevents['fateCE'].shift(-1).eq('0')) & (mtevents['event'].shift(-2).ne('MERGER')) & (mtevents['event'].shift(-2).ne('SWALLOWED')))][['name','time']] #RLO->CE surviving
    RLOiCEmerge = mtevents.loc[(mtevents['event'].eq('RLO_BEGIN')) & (mtevents['event'].shift(-1).eq('CE')) & ((mtevents['fateCE'].shift(-1).eq('1')) | (mtevents['event'].shift(-2).eq('MERGER')) | (mtevents['event'].shift(-2).eq('SWALLOWED')))][['name','time']] #RLO->CE merging
    RLOiRLOf = mtevents.loc[(mtevents['event'].eq('RLO_BEGIN')) & (mtevents['event'].shift(-1).eq('RLO_END'))][['name','time']]  # stable RLO

    COLLCEsurv = mtevents.loc[(mtevents['event'].eq('COLLISION')) & (mtevents['event'].shift(-1).eq('CE')) & ((mtevents['fateCE'].shift(-1).eq('0')) & (mtevents['event'].shift(-2).ne('MERGER')) & (mtevents['event'].shift(-2).ne('SWALLOWED')))][['name','time']]
    COLLCEmerge = mtevents.loc[(mtevents['event'].eq('COLLISION')) & (mtevents['event'].shift(-1).eq('CE')) & ((mtevents['fateCE'].shift(-1).eq('1')) | (mtevents['event'].shift(-2).eq('MERGER')) | (mtevents['event'].shift(-2).eq('SWALLOWED')))][['name','time']]
    COLLMERGE = mtevents.loc[(mtevents['event'].eq('COLLISION')) & (mtevents['event'].shift(-1).eq('MERGER'))][['name','time']]

    WR0 = singleevents.loc[singleevents['order'].eq(6)]   # use number assigned in LogToEvents function
    WR1 = singleevents.loc[singleevents['order'].eq(7)]   # use number assigned in LogToEvents function
    SN0 = singleevents.loc[singleevents['order'].eq(8)]   # use number assigned in LogToEvents function
    SN1 = singleevents.loc[singleevents['order'].eq(9)]   # use number assigned in LogToEvents function
    rem0 = singleevents.loc[singleevents['order'].eq(10)]   # use number assigned in LogToEvents function
    rem1 = singleevents.loc[singleevents['order'].eq(11)]   # use number assigned in LogToEvents function
    
    dflist = [RLOiSWALLOWED,RLOiMERGE,RLOiCEsurv,RLOiCEmerge,RLOiRLOf,COLLCEsurv,COLLCEmerge,COLLMERGE,WR0,WR1,SN0,SN1,rem0,rem1]
    BEventlist = [17,10,11,12,5,14,15,13,30,40,50,60,70,80]
    for df,BEvent in zip(dflist,BEventlist):
        df['BEvent'] = BEvent
    history = pd.concat(dflist)[['name','time','BEvent','phase_comp','rem_type']].sort_values(['name','time'])

    BEcond_list = Condition(history)
    BEfates = ['merging', 'CEsurv','RLOstable','WR0','WR1','SN0','SN1','Rem0','Rem1']
    for BEcond, BEfate in zip(BEcond_list,BEfates):
        history.loc[BEcond,'event'] = BEfate
    return history
