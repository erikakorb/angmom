import warnings
warnings.simplefilter(action='ignore')
import os
import subprocess
import time
from multiprocessing import Pool
import dask.dataframe as dd
import pandas as pd
import numpy as np
import lib_out_extract    	# lib_out_extract.py, in the same folder, contains auxiliary functions to extract MSBH and BBH dataframes




if __name__ == "__main__":
    # input parameters
    sevn_version = 'sevn_25ago'              # SEVN version adopted 
    angmoms = ['1.5']        # angular momentum model      # angmoms= ['-2', '-1', '1','1.5']
    Zs = ['0.00142']              # metallicity                 # Zs = ['0.00142','0.000142']
    tides = ['tides_simple','disabled']     # tides model	           # tides = ['tides_simple', 'disabled']

    
    # prepare lists
    name_list1 = ['BHBH','BHBH_bound','BHBH_GW']
    name_list2 = ['MSBH']
    word_list = ['BHBH', 'bound BHBH', 'GW-BHBH']

    # input path
    path_to_sevn = f'/tank1/korb/angmom/{sevn_version}/sevn'                   # original folder with SEVN output
    
    # output paths
    path_results = f'/tank1/korb/angmom/v_{sevn_version}'            # path to new folder with all useful results

    path_txt_results = f'{path_results}/results'                        # txt files with results in human-readable form
    path_df_prog = f'{path_results}/progenitors'                        # progenitors
    path_df_rem = f'{path_results}/remnants'                            # remnants
    path_df_in = f'{path_results}/initial'                              # initial timestamp for MSBH phase
    path_df_fin = f'{path_results}/final'                               # final timestamp for MSBH phase

    path_list = [path_results,path_txt_results,path_df_prog,path_df_rem,path_df_in,path_df_fin]
    for path in path_list:
    	os.makedirs(path, exist_ok=True)


    ##############################################################
    starttot = time.time()
    for angmom in angmoms:
        for Z in Zs:
            for tide in tides:
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
                output_list = [filename for filename in file_names if filename.split('_')[0]=='output']
                evolved_list = [filename for filename in file_names if filename.split('_')[0]=='evolved']
                
                print('Parallel output analysis')
                startparallel=time.time()
     
                with Pool(processes = 6) as p:
                    print('Evolved files')
                    df_evo_list = p.starmap(lib_out_extract.ReadEvolved,zip(evolved_list,[path_to_sevn_output]*len(evolved_list)))
                    print('Ouput files')
                    df1_list_of_lists = p.starmap(lib_out_extract.OutputBBH,zip(output_list,[path_to_sevn_output]*len(output_list)))
                    df2_list_of_lists = p.starmap(lib_out_extract.OutputMSBH,zip(output_list,[path_to_sevn_output]*len(output_list)))

                # concatenate dataframes to a single one per type
                evolved = pd.concat(df_evo_list, ignore_index=True)
                df_list1 = lib_out_extract.ConcatDfLists(df1_list_of_lists)
                df_list2 = lib_out_extract.ConcatDfLists(df2_list_of_lists)

                endparallel=time.time()
                print('Parallel output analysis concluded in [min]: ', (endparallel-startparallel)/60.)
                print('------------------------')



                #######################################################
                ##      Write dataframes and output analysis         ##
                #######################################################
                startwrite=time.time()

                with open(f'{path_txt_results}/{setname}_results.txt', 'w') as f:
                    f.write(f'Set with angmom: {angmom} \n')
                    f.write('Number of simulated binaries: ' + str(len(evolved.index)) + '\n')
                    for df2,name2 in zip(df_list2,name_list2):
                        f.write(f'Number of {name2} systems: ' + str(len(df2.drop_duplicates(subset='ID').index)) + '\n')

                    for df1,name1,word in zip(df_list1,name_list1,word_list):
                        print(f'Writing dataframes with {word}')
                        f.write('----------------------------------------\n')
                        f.write(f'Number of {word} systems: ' + str(len(df1.drop_duplicates(subset='ID').index)) + '\n')
                        
                        # select progenitors and remnants
                        remnant = df1.drop_duplicates(['ID'],keep='last')
                        IDs = set(remnant.ID.to_list())  # list of IDs
                        progenitor = evolved.query('ID in @IDs')

                        # ensure folders exist
                        os.makedirs(f'{path_df_prog}/{name1}', exist_ok=True)
                        os.makedirs(f'{path_df_rem}/{name1}', exist_ok=True)
                        
                        # write dataframes
                        progenitor.to_csv(f'{path_df_prog}/{name1}/{setname}_p_{name1}.csv')
                        remnant.to_csv(f'{path_df_rem}/{name1}/{setname}_r_{name1}.csv')

                        for df2,name2 in zip(df_list2,name_list2):
                            # ensure folders exist
                            os.makedirs(f'{path_df_in}/{name2}', exist_ok=True)
                            os.makedirs(f'{path_df_fin}/{name2}',exist_ok=True)
                            os.makedirs(f'{path_df_prog}/{name2}',exist_ok=True)

                            os.makedirs(f'{path_df_prog}/{name1}_{name2}', exist_ok=True)
                            os.makedirs(f'{path_df_in}/{name1}_{name2}', exist_ok=True)
                            os.makedirs(f'{path_df_fin}/{name1}_{name2}',exist_ok=True)
                            os.makedirs(f'{path_df_rem}/{name1}_{name2}',exist_ok=True)

                            # select and store only initial and final MSBH-like phase for all candidates
                            init = df2.drop_duplicates('ID', keep='first')
                            fin = df2.drop_duplicates('ID', keep='last')
                            init.to_csv(f'{path_df_in}/{name2}/{setname}_i_{name2}.csv')
                            fin.to_csv(f'{path_df_fin}/{name2}/{setname}_f_{name2}.csv')
                            
                            dfID = set(init.ID.to_list())
                            p = evolved.query('ID in @dfID')
                            p.to_csv(f'{path_df_prog}/{name2}/{setname}_p_{name2}.csv') 
                            
                            # ID of binaries that belong in both categories
                            binID = set(df1['ID']).intersection(df2['ID'])

                            # progenitors and remnants of such binaries
                            prog = evolved.query('ID in @binID')
                            remnant_timesteps = df1.query('ID in @binID')
                            rem = remnant_timesteps.drop_duplicates('ID', keep='last')
                            
                            # initial and final timestamps of WRBH phase evolution for each subgroup
                            timesteps = df2.query('ID in @binID')
                            initial = timesteps.drop_duplicates('ID', keep='first')
                            final = timesteps.drop_duplicates('ID', keep='last')

                            # write dataframes
                            prog.to_csv(f'{path_df_prog}/{name1}_{name2}/{setname}_p_{name1}_{name2}.csv')
                            rem.to_csv(f'{path_df_rem}/{name1}_{name2}/{setname}_r_{name1}_{name2}.csv')                                
                            initial.to_csv(f'{path_df_in}/{name1}_{name2}/{setname}_i_{name1}_{name2}.csv')
                            final.to_csv(f'{path_df_fin}/{name1}_{name2}/{setname}_f_{name1}_{name2}.csv')

                            # write useful results
                            f.write(f'Number of {word} systems formed after {name2}: ' + str(len(binID)) + '\n')

                endwrite=time.time()
                print('Time to write the useful dataframes [min]: ', (endwrite-startwrite)/60.)
                print('------------------------')

                ### output messages to check ###
                end=time.time()
                print(f'Time to elaborate all the {setname} output files [min]: ', (end-start)/60.)
                print(f'The files have been correctly stored in {path_results} \n')

    endtot = time.time()
    print('************************************************************************************************')
    print('************************************************************************************************ \n')
    print('     All SEVN output folders have been correctly analyzed in [min]: ', (endtot-starttot)/60., '\n')
    print('************************************************************************************************')
    print('************************************************************************************************')
    
    

