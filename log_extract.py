import warnings
warnings.simplefilter(action='ignore')
import os
import subprocess
import time
from multiprocessing import Pool
import dask.dataframe as dd
import pandas as pd
import numpy as np
import lib_log_extract    	# lib_log_extract.py, in the same folder, contains auxiliary functions to extract mass transfer channels




if __name__ == "__main__":
    # input parameters
    sevn_version = 'sevn_25ago'              # SEVN version adopted 
    angmoms = ['1.5']        # angular momentum model      # angmoms= ['-2', '-1', '1','1.5']
    Zs = ['0.00142']             # metallicity                 # Zs = ['0.00142','0.000142']
    tides = ['tides_simple','disabled']     # tides model	           # tides = ['tides_simple', 'disabled']

    
    # prepare lists
    name_list1 = ['BHBH','BHBH_bound','BHBH_GW']
    name_list2 = ['MSBH']
    word_list = ['BHBH', 'bound BHBH', 'GW-BHBH']
    
    # input    
    path_to_sevn = f'/tank1/korb/angmom/{sevn_version}/sevn'                   # original folder with SEVN output

    ##############################################################
    starttot = time.time()
    for angmom in angmoms:
        for Z in Zs:
            for tide in tides:
                #######################################################
                ##          Extract infos from raw outputs           ##
                #######################################################
                starttot=time.time()
                setname = f'Z{Z}_{angmom}_{tide}'
                print('################################ \n')
                print(f'     {setname} \n')
                print('################################')

                # input path
                path_to_sevn_output = f'{path_to_sevn}/run_scripts/{setname}/sevn_output/'  # original sevn_output folder with SEVN output

                # prepare new folders with results
                path_log_history = f'/tank1/korb/angmom/v_{sevn_version}/logs/'            # path to new folder with all useful results
                os.makedirs(path_log_history, exist_ok=True)

                # identify data of interest in SEVN log files
                file_names =os.listdir(path_to_sevn_output)
                logfile_list = [filename for filename in file_names if filename.split('_')[0]=='logfile']
                
                print('Parallel output analysis')
                startparallel=time.time()
     
                with Pool(processes = 6) as p:
                    print('Log files')
                    mt_history_list = p.starmap(lib_log_extract.GetHistory,zip(logfile_list,[path_to_sevn_output]*len(logfile_list) ) )

                # concatenate dataframes to a single one per type
                mt_history = pd.concat(mt_history_list)

                endparallel=time.time()
                print('Parallel output analysis concluded in [min]: ', (endparallel-startparallel)/60.)
                print('------------------------')


                # write logfiles analysis   
                startwrite=time.time()
                print('Writing mass transfer analysis')
                mt_history.to_csv(f'{path_log_history}/{setname}_history_full.csv', index=False)

                endwrite=time.time()
                print('Time to write the mass transfer dataframes [min]: ', (endwrite-startwrite)/60.)
                print('------------------------')

    endtot = time.time()
    print('************************************************************************************************')
    print('************************************************************************************************ \n')
    print('     All SEVN logfiles have been correctly analyzed in [min]: ', (endtot-starttot)/60., '\n')
    print('************************************************************************************************')
    print('************************************************************************************************')
    
    

