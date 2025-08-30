import numpy as np
import pandas as pd
import dask.dataframe as dd

# Functions to extract MS-BH and BBH dataframes from evolved and output files
# These functions need to be imported in reductionpipeline1.py


#######################################################
##          Auxiliary functions for dataframes       ##
#######################################################

# Selects first and last timestep of each binary, identified by its ID
# If in input there is a dask dataframe, it is necessary to set to True the optional argument [dask]
#
def first_last_ID(df,dask=False):   
    first = df.drop_duplicates(subset=['ID'], keep='first',ignore_index=True)
    last = df.drop_duplicates(subset=['ID'], keep='last',ignore_index=True)    
    if dask ==True:
        first = first.compute(num_workers = 4)
        last = last.compute(num_workers = 4)  
    df = pd.concat([first,last]).sort_values(['ID','BWorldtime'])    
    return df

# Concatenates dataframes in the same position in lists of lists
# Dataframes are iteratively concatenated to the ones in the 0th list
#
def ConcatDfLists(list_of_lists):
    for l in range(1,len(list_of_lists)):
        for d in range(0,len(list_of_lists[0])):
            list_of_lists[0][d] = pd.concat([list_of_lists[0][d], list_of_lists[l][d]],ignore_index=True)
    return list_of_lists[0]



##################################################
##        Binary configuration setup            ##
##################################################

# Classifiy BH if object satisfies
# - PhaseBSE = 14 	BH formation
# - Mass > MmaxNS        massive enough to be a BH in place of NS
# In SEVN defalut, a CO of 3 Msun is still a NS while a CO of 3.0001 is a BH
# [MmaxNS] sets max mass for NS classification: compact objects with M<=MmaxNS are NS
#
def BHcondition(df):
    MmaxNS= 3.   
    BH0 = (df['PhaseBSE_0'] == 14) & (df['Mass_0'] > MmaxNS)
    BH1 = (df['PhaseBSE_1'] == 14) & (df['Mass_1'] > MmaxNS)
    return BH0,BH1

# Extract BBH gravitationally bound and GW-merging
# The GW-merging condition considers time to merge + lifetime of binary
# [H0] is the Hubble time in Myr (by default H0=14000)
#
def boundGW(df):
    H0=14000  # Myr
    df_bound=df[~np.isnan(df['Semimajor'])].drop_duplicates(subset=['ID'], keep='last')
    df_GW=df_bound.loc[(df_bound['GWtime'] + df_bound['BWorldtime'])  <=H0]
    return df_bound, df_GW



# Preliminary set condition to identify a Main Sequence (MS) star 
# in SEVN corresponds to PhaseBSE =1
#
def MScondition(df):
    MS0 = (df['PhaseBSE_0'] == 1)
    MS1 = (df['PhaseBSE_1'] == 1)
    return MS0,MS1


# Extract infos from MSBH - a dask dataframe input
# keeps only first and last timestep for each ID to reduce ram occupation
# [n] sets the number of CPUs for dask dataframe computation
#
def MSBHextractif(MSBH):
    n=4                                                         
    MSBH = MSBH.compute(num_workers=n)
    MSBHif = first_last_ID(MSBH)                        
    return MSBHif



#######################################################
##                 Output analysis                   ##
#######################################################


# Read single evolved_N.csv file
# And use Kepler 3rd law to convert period into yrs
# [G4pi2] is G/(4 pi^2) in units of R_sun^3/(M_sun yr^2)
#
def ReadEvolved(evolvedname,path_to_sevn_output):
    G4pi2 =  9953108.1
    evo = pd.read_csv(path_to_sevn_output+evolvedname, sep='\s+').rename(columns = {'#ID':'ID'})
    evo['Period'] = np.sqrt(evo['a']**3 / (G4pi2*(evo['Mass_0']+evo['Mass_1'])))
    return evo



# Extract and classify binaries from a single raw output_N.csv file of SEVN
#
def OutputBBH(outputname,path_to_sevn_output):
    out=dd.read_csv(path_to_sevn_output+outputname, blocksize= 128* 1024 * 1024) 

    BH0,BH1 = BHcondition(out)
    BHBH = out.loc[BH0 & BH1]                                                  
    BHBH = first_last_ID(BHBH,dask=True)
    BHBH_bound, BHBH_GW = boundGW(BHBH)
    return [BHBH,BHBH_bound,BHBH_GW]


# Extract and classify MSBH binaries from single raw output_N.csv file of SEVN
# First compute dask dataframes to get all MSBHs and just initial and final timesteps
# Eventually delete MSBH and MSBH dataframes to free memory
# note: maintain the list output format so it easy expandable with other binary configurations and maintains consitency with the concatdflists
def OutputMSBH(outputname,path_to_sevn_output):
    out=dd.read_csv(path_to_sevn_output+outputname, blocksize= 128* 1024 * 1024) 

    MS0, MS1 = MScondition(out)                 
    BH0,BH1 = BHcondition(out)     
    MSBH = out.loc[ (BH0 & MS1) | (MS0 & BH1) ]                
    MSBHif = MSBHextractif(MSBH)
    del MSBH
    return [MSBHif]


