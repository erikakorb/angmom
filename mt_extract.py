import warnings
warnings.simplefilter(action='ignore')
import time
import os
import pandas as pd
import numpy as np
import lib_mt_extract    # import auxiliary library


if __name__ == "__main__":

    # input parameters
    sevn_version = 'sevn_25ago'              # SEVN version adopted 
    angmoms = ['1.5']        # angular momentum model      # angmoms= ['-2', '-1', '1','1.5']
    Zs = ['0.00142','0.000142']              # metallicity                 # Zs = ['0.00142','0.000142']
    tides = ['tides_simple','disabled']     # tides model	           # tides = ['tides_simple', 'disabled']

    fate_type = 'BHBH_GW'
    binary_type = 'MSBH'
    #label_name = f'{fate_type}_from_all_prog'       # if label_name in {f'{fate_type}_from_all_prog',f'{fate_type}_from_{binary_type}',f'{fate_type}_not_from_{binary_type}'}
    #label_name = f'{binary_type}_to_all_fates'       # elif label_name in {f'{binary_type}_to_all_fates',f'{binary_type}_to_{fate_type}',f'{binary_type}_not_to_{fate_type}'}
    label_names = [f'{fate_type}_from_all_prog',f'{binary_type}_to_all_fates']
    
    
    starttot = time.time()
    for label_name in label_names:
    	for angmom in angmoms:
        	for Z in Zs:
            		for tide in tides:
                		starttot=time.time()
                		setname = f'Z{Z}_{angmom}_{tide}'
                		print('################################ \n')
                		print(f'     {setname} - {label_name} \n')
                		print('################################')

                		# set paths
                		path_to_sevn = f'/tank1/korb/angmom/{sevn_version}/sevn/'  
                		path_mt_folder = f'/tank1/korb/angmom/v_{sevn_version}/mt_{label_name}' # new folder with mt summary

                		os.makedirs(path_mt_folder, exist_ok=True)

                		# get mt_summary file    
                		NwithSN = lib_mt_extract.GetHistoryFirstSN(sevn_version,setname,fate_type,binary_type, label_name)

                		# write output file
                		NwithSN.to_csv(f'{path_mt_folder}/{setname}_mt_{label_name}.csv')

