import re
import os
import pandas as pd
import numpy as np
import lib_mt_extract        # for auxiliary functions

sevn_version = 'sevn_25ago'              # SEVN version adopted 
angmoms = ['1.5']        # angular momentum model      # angmoms= ['-2', '-1', '1','1.5']
Zs = ['0.00142','0.000142']              # metallicity                 # Zs = ['0.00142','0.000142']
tides = ['tides_simple','disabled']     # tides model	           # tides = ['tides_simple', 'disabled']
#sim_labels = ['IsoT','JeansT','LminT','LmaxT',
#              'IsoNT','JeansNT','LminNT','LmaxNT']
sim_labels = ['LmaxT','LmaxNT']


fate_type = 'BHBH_GW'
binary_type = 'MSBH'
label_names = [f'{fate_type}_from_all_prog', f'{binary_type}_to_all_fates']
#label_name = f'{fate_type}_from_all_prog'       # if label_name in {f'{fate_type}_from_all_prog',f'{fate_type}_from_{binary_type}',f'{fate_type}_not_from_{binary_type}'}
#label_name = f'{binary_type}_to_all_fates'       # elif label_name in {f'{binary_type}_to_all_fates',f'{binary_type}_to_{fate_type}',f'{binary_type}_not_to_{fate_type}'}


for Z in Zs:
    for label_name in label_names:
        if label_name == f'{fate_type}_from_all_prog': # BBH
            name_list = ['N_MS_comp','N_noWR_comp', 'N_WR_comp',
                         'N_RLO_preSN','N_CE_preSN','N_noMTpreSN',
                         'N_RLO_postSN','N_CE_postSN','N_noMTpostSN']
        elif label_name == f'{binary_type}_to_all_fates':  # MSBH
            name_list = ['N_MSBH_to_BBH',
                         'N_RLO_postSN_to_BBH','N_CE_postSN_to_BBH','N_noMTpostSN_to_BBH',
                         'N_merge','N_WD','N_NS','N_BH',
                         'N_RLO_preSN','N_CE_preSN','N_noMTpreSN']
        
        matrix = np.zeros((len(sim_labels),len(name_list)+2)) # initialize matrix to store results
        i=0
        for tide in tides:
            for angmom in angmoms:
            
                setname = f'Z{Z}_{angmom}_{tide}'
                print('################################ \n')
                print(f'     {setname} \n')
                print('################################')

                # get number of all simulated binaries
                path_to_set = f'./v_{sevn_version}'
                path_res = f'{path_to_set}/results/{setname}_results.txt'
                row = []
                textfile = open(path_res,'r').read()
                Nsim = float(re.findall(f'Number of simulated binaries: (\d+)', textfile)[0])

                # read mass transfer summary
                df = pd.read_csv(f'./v_{sevn_version}/mt_{label_name}/{setname}_mt_{label_name}.csv')                
                N_all = len(df)
                
                # info for a BBH configuration
                if label_name == f'{fate_type}_from_all_prog':
                    N_prog_MS = len(df.loc[df.phase_comp == 1])       # progenitor in the MS-BH configuration
                    N_prog_HG_coreHe_AGB = len(df.loc[((df.phase_comp == 2) & (df.phase_comp == 3) & (df.phase_comp == 4) & (df.phase_comp == 5) & (df.phase_comp == 6) & (df.WR1_prefirstSN != 'WR1'))])   # progenitor in other configuration, but not a WR
                    N_prog_WR = len(df.loc[df.WR1_prefirstSN == 'WR1'])  # progenitor in the BH - WR configuration   

                    N_ch_RLOpreSN = len(df.loc[((df.NRLOpreSN > 0) & (df.NCEpreSN == 0))])
                    N_ch_CEpreSN = len(df.loc[df.NCEpreSN > 0])
                    N_ch_noMTpreSN = len(df.loc[((df.NRLOpreSN == 0) & (df.NCEpreSN == 0))])                                
                
                    N_ch_RLOpostSN = len(df.loc[((df.NRLOpostSN > 0) & (df.NCEpostSN == 0))])
                    N_ch_CEpostSN = len(df.loc[df.NCEpostSN > 0])
                    N_ch_noMTpostSN = len(df.loc[((df.NRLOpostSN == 0) & (df.NCEpostSN == 0))])

                    N_list = [N_prog_MS,N_prog_HG_coreHe_AGB,N_prog_WR,
                              N_ch_RLOpreSN,N_ch_CEpreSN,N_ch_noMTpreSN,
                              N_ch_RLOpostSN, N_ch_CEpostSN, N_ch_noMTpostSN]  
                                                                                     
                                             
                # info for MSBH configuration
                elif label_name == f'{binary_type}_to_all_fates':
                    names_MSBHs, names_BBH_from_binary, names_binary_no_BBHs = lib_mt_extract.GetNamesMSBH(sevn_version,setname,fate_type,binary_type)
                    df_MSBH_to_BBH = df.loc[df.name.isin(names_BBH_from_binary)]
                    N_MSBH_to_BBH = len(df_MSBH_to_BBH)
                    
                    # channels post SN only for MSBH becoming merging BBHs
                    N_ch_RLOpostSN = len(df_MSBH_to_BBH.loc[((df_MSBH_to_BBH.NRLOpostSN > 0) & (df_MSBH_to_BBH.NCEpostSN == 0))])
                    N_ch_CEpostSN = len(df_MSBH_to_BBH.loc[df_MSBH_to_BBH.NCEpostSN > 0])
                    N_ch_noMTpostSN = len(df_MSBH_to_BBH.loc[((df_MSBH_to_BBH.NRLOpostSN == 0) & (df_MSBH_to_BBH.NCEpostSN == 0))])
                    
                    merge_cond = ( (df.Nmerge > 0) | (df.rem_type == -1) )
                    surv_cond = ( (df.Nmerge == 0) | (df.rem_type != -1) )                        
                    N_ch_merge = len(df.loc[merge_cond])
                    N_ch_WD = len(df.loc[((df.rem_type == 1) | (df.rem_type == 2) | (df.rem_type == 3)) & surv_cond])
                    N_ch_NS = len(df.loc[((df.rem_type == 4) | (df.rem_type == 5)) & surv_cond])
                    N_ch_BH = N_all - N_MSBH_to_BBH - N_ch_merge - N_ch_WD - N_ch_NS    # other fates presume the formation of a BH in a wide or broken orbit, not merged
                    
                    # channels pre SN for all MSBH
                    N_ch_RLOpreSN = len(df.loc[((df.NRLOpreSN > 0) & (df.NCEpreSN == 0))])
                    N_ch_CEpreSN = len(df.loc[df.NCEpreSN > 0])
                    N_ch_noMTpreSN = len(df.loc[((df.NRLOpreSN == 0) & (df.NCEpreSN == 0))])   
                    
                    N_list = [N_MSBH_to_BBH,
                              N_ch_RLOpostSN, N_ch_CEpostSN, N_ch_noMTpostSN,
                              N_ch_merge, N_ch_WD, N_ch_NS, N_ch_BH,
                              N_ch_RLOpreSN,N_ch_CEpreSN,N_ch_noMTpreSN]
                matrix[i] = np.array([Nsim,N_all]+N_list,dtype=int)   # add row
                i+=1
                
        dfN = pd.DataFrame(matrix,columns = ['N_sim',f'N_{label_name}'] + name_list,index=sim_labels)
        
        os.makedirs(f'./v_{sevn_version}/summary/', exist_ok=True)
        dfN.to_csv(f'./v_{sevn_version}/summary/N_{label_name}_Z{Z}.csv')

