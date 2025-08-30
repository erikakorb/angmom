import warnings
warnings.simplefilter(action='ignore')
import os
import subprocess
import time
from multiprocessing import Pool
import dask.dataframe as dd
import pandas as pd
import numpy as np

# Read single evolved_N.csv file
# And get the simulated mass
#
def GetMass(evolvedname,path_to_sevn_output):
    evo = pd.read_csv(path_to_sevn_output+evolvedname, sep='\s+').rename(columns = {'#ID':'ID'})
    Mtot = evo['Mass_0'].sum() + evo['Mass_1'].sum() 
    return Mtot




if __name__ == "__main__":

    # input parameters
    sevn_version = 'sevn_25ago'              # SEVN version adopted 
    angmoms = ['1.5']        # angular momentum model      # angmoms= ['-2', '-1', '1','1.5']
    Zs = ['0.00142','0.000142']              # metallicity                 # Zs = ['0.00142','0.000142']
    tides = ['tides_simple','disabled']     # tides model             # tides = ['tides_simple', 'disabled']
    #sim_labels = ['IsoT','JeansT','LminT','LmaxT',
    #              'IsoNT','JeansNT','LminNT','LmaxNT']
    sim_labels = ['LmaxT','LmaxNT']

    # input path
    path_to_sevn = f'/tank1/korb/angmom/{sevn_version}/sevn'                   # original folder with SEVN output

    # output paths
    path_results = f'/tank1/korb/angmom/v_{sevn_version}'            # path to new folder with all useful results
    path_mass_sim = f'{path_results}/mass_simulated'                 

    os.makedirs(path_mass_sim, exist_ok=True)


    ##############################################################
    matrix = np.zeros((len(sim_labels),len(Zs))) # initialize matrix to store results
    
    starttot = time.time()
    icol = 0
    for Z in Zs:
        irow = 0
        for tide in tides:
            for angmom in angmoms:
                #######################################################
                ##          Extract infos from raw outputs           ##
                #######################################################
                start=time.time()
                setname = f'Z{Z}_{angmom}_{tide}'
                print('################################ \n')
                print(f'     {setname} \n')
                print('################################')

                # input path
                path_to_sevn_output = f'{path_to_sevn}/run_scripts/{setname}/sevn_output/'  # original sevn_output folder with SEVN output

                # identify data of interest in SEVN output and evolved files
                file_names =os.listdir(path_to_sevn_output)
                evolved_list = [filename for filename in file_names if filename.split('_')[0]=='evolved']

                print('Parallel output analysis')
                startparallel=time.time()

                with Pool(processes = 6) as p:
                    print('Evolved files')
                    Mtot_list = p.starmap(GetMass,zip(evolved_list,[path_to_sevn_output]*len(evolved_list)))                  

                # concatenate dataframes to a single one per type
                Mtot_sum = sum(Mtot_list)
                matrix[irow][icol] = Mtot_sum
                irow += 1

        icol += 1  

df = pd.DataFrame(matrix,columns = Zs,index=sim_labels)
df.to_csv(f'{path_mass_sim}/mass_simulated.csv')

