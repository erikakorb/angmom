# %%
"""
Simple script to produce a grid plot similar to Fig. 2 in Gallegos-Garcia+21 paper
(https://ui.adsabs.harvard.edu/abs/2021ApJ...922..110G/abstract)
"""
from sevnpy.sevn import SEVNmanager,sevnwrap,Star
from sevnpy.io import readlogstring
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import calclib as cl        # custom library with routines for tides, etc.


# initial period in days
def run_single(dict, m1,q,per,Zdonor):
    m2 = m1*q   # BH mass
    a = calc_a(per / 365.25,m1,m2)  # initial semimajor axis

    # initialize SEVN
    SEVNmanager.init( dict)

    # run
    dfP,lPdic=sevnwrap.evolve_binary(a,0.,m1,Zdonor,m2,Zdonor,
                                            spin_0 = 0,
                                            star_flag_1="BH", snmodel='delayed')

    SEVNmanager.close()

    # add effective radius i.e. stellar radius is limited to Roche lobe radius
    #  while RLO active
    dfP = pd.DataFrame(dfP)
    dfP['Radius_eff'] = np.select([(dfP.Radius_0 >= dfP.RL0), (dfP.Radius_0 < dfP.RL0)], 
                                [dfP['RL0'], dfP['Radius_0']], 
                                default=np.nan) 
    return dfP,lPdic




# %% [markdown]
# # Functions to plot grids

# %%
def find_channel_bugfix(df,log):
    """
    Find the formation channel based on the SEVN log
    """
    alog=readlogstring(log,events=["CE","RLO_BEGIN","MERGER","COLLISION","SWALLOWED"],capturing_cols=("time","event"))
    dfCE=alog["CE"]
    dfRL=alog["RLO_BEGIN"]
    dfM=alog["MERGER"]
    dfC=alog["COLLISION"]
    dfS=alog["SWALLOWED"]
    nosurv = (len(dfM)>0 or len(dfC)>0 or len(dfS)>0 or (find_status_bugfix(df)==-1))
    surv = (len(dfM)==0 and len(dfC)==0 and len(dfS)==0 and (find_status_bugfix(df) != -1))

    if len(dfRL)>0 and len(dfCE)==0 and surv:
        return 1
    elif len(dfCE)>0 and surv:
        return 2
    elif len(dfCE)>0 and nosurv:
        return 3
    elif len(dfRL)>0 and nosurv:
        return 4
    elif len(dfCE)==0 and len(dfRL)==0  and nosurv:
        return 5
    elif len(dfCE)==0 and len(dfRL)==0  and surv:
        return 6



def find_status_bugfix(df):
    """
    Find the final status of the  binary  based on the results of the SEVN evolution
    (-1 Destroyed binary, 1 Merge within an Hubble time, 0 Not merge within an Hubble time)
    """
    # ad hoc fix for systems where RLobe < He core radius, thus, should merge but for some bug are resetted as wide
    Mergers = np.where (df['RL0'] <= df['RHE_0'])[0]   # selects rows where RL < RHE. If at least one, system should have merged
    
    if (np.isnan(df["GWtime"][-1])) or (len(Mergers) > 0):
        return -1
    else:
        tdel= df["Worldtime"][-1] + df["GWtime"][-1]
        if tdel<=14000:return 1
        else: return 0

def find_status_single_bugfix(df):
    """
    Find the final status of the  binary  based on the results of the SEVN evolution
    (-1 Destroyed binary, 1 Merge within an Hubble time, 0 Not merge within an Hubble time)
    """
    # ad hoc fix for systems where RLobe < He core radius, thus, should merge but for some bug are resetted as wide
    Mergers = np.where (df['RL0'].values <= df['RHE_0'].values)[0]   # selects rows where RL < RHE. If at least one, system should have merged
    
    
    if np.isnan(df["GWtime"].values[-1]) or (len(Mergers) > 0):
        return -1
    else:
        tdel= df["Worldtime"].values[-1] + df["GWtime"].values[-1]
        if tdel<=14000:return 1
        else: return 0

def find_label_bugfix(df,log):
    """
    Find the label for the grid plot
    """
    ch = find_channel_bugfix(df,log)
    st = find_status_bugfix(df)

    if ch==3:
        return "MergerCE"
    elif ch==4:
        return "MergerRL"
    elif ch==2 and st==1:
        return "BBHmCE"
    elif ch==1 and st==1:
        return "BBHmSMT"
    #elif ch==5 or (ch==6 and st==1):
    #    return "other"
    else:
        return "wide"


################## end of bug fix functions ##########

def find_status(df):
    """
    Find the final status of the  binary  based on the results of the SEVN evolution
    (-1 Destroyed binary, 1 Merge within an Hubble time, 0 Not merge within an Hubble time)
    """
    
    if (np.isnan(df["GWtime"][-1])):
        return -1
    else:
        tdel= df["Worldtime"][-1] + df["GWtime"][-1]
        if tdel<=14000:return 1
        else: return 0

def find_channel(df,log):
    """
    Find the formation channel based on the SEVN log
    """
    alog=readlogstring(log,events=["CE","RLO_BEGIN","MERGER","COLLISION","SWALLOWED"],capturing_cols=("time","event"))
    dfCE=alog["CE"]
    dfRL=alog["RLO_BEGIN"]
    dfM=alog["MERGER"]
    dfC=alog["COLLISION"]
    dfS=alog["SWALLOWED"]
    nosurv = (len(dfM)>0 or len(dfC)>0 or len(dfS)>0 or (find_status(df)==-1))
    surv = (len(dfM)==0 and len(dfC)==0 and len(dfS)==0 and (find_status(df) != -1))

    if len(dfRL)>0 and len(dfCE)==0 and surv:
        return 1
    elif len(dfCE)>0 and surv:
        return 2
    elif len(dfCE)>0 and nosurv:
        return 3
    elif len(dfRL)>0 and nosurv:
        return 4
    elif len(dfCE)==0 and len(dfRL)==0  and nosurv:
        return 5
    elif len(dfCE)==0 and len(dfRL)==0  and surv:
        return 6

def find_status_single(df):
    """
    Find the final status of the  binary  based on the results of the SEVN evolution
    (-1 Destroyed binary, 1 Merge within an Hubble time, 0 Not merge within an Hubble time)
    """

    if np.isnan(df["GWtime"].values[-1]):
        return -1
    else:
        tdel= df["Worldtime"].values[-1] + df["GWtime"].values[-1]
        if tdel<=14000:return 1
        else: return 0


def find_label(df,log):
    """
    Find the label for the grid plot
    """
    ch = find_channel(df,log)
    st = find_status(df)

    if ch==3:
        return "MergerCE"
    elif ch==4:
        return "MergerRL"
    elif ch==2 and st==1:
        return "BBHmCE"
    elif ch==1 and st==1:
        return "BBHmSMT"
    #elif ch==5 or (ch==6 and st==1):
    #    return "other"
    else:
        return "wide"


def calc_a(per,m1,m2):
    """
    From period to semi-major axss
    :param per:  Periods in  year
    :param m1: Mass in Msun
    :param m2:  Mass in Msun
    :return: semi-major axis in Rsun
    """
    G = 3.925125598496094e8
    return (per*per*(G*(m1+m2))/(np.pi*np.pi*4))**(1./3.)

def plot_fig(ql,pl,ll,Mdonor,Zdonor):    
    plt.rcParams.update({
    "font.family": "cmr10",   # computer modern font, as in latex
    "font.size": 17,
    "axes.labelsize" : 17,
    "xtick.labelsize" : 17,
    "ytick.labelsize" : 17,
    "legend.fontsize": 15
    })

    fig,ax=plt.subplots(1,1,figsize=(5.4,5.8))
    idx=ll=="MergerRL"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="lightgrey",marker="s",label=r"Merger after SMT")
    idx=ll=="MergerCE"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="gold",marker="s",label=r"Merger after CE")
    idx=ll=="wide"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="lightskyblue",marker="s",label=r"Wide binary")
    idx=ll=="BBHmSMT"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="blue",marker="s",label=r"BBH merger after SMT")
    idx=ll=="BBHmCE"
    plt.scatter(ql[idx],np.log10(pl[idx]),c="orange",marker="s",label=r"BBH merger after CE")
    #plt.suptitle(r'$ M_{\rm donor}=$'+ str(Mdonor) +r'$ M_\odot$, Z='+ str(Zdonor))
    fig.legend(ncol=2, bbox_to_anchor=(0.98,1.12),frameon=False)
    plt.xlabel(r"$q_{\rm i} = M_{\rm BH}$/$M_{\rm donor}$")
    plt.ylabel(r"$\log_{10} (P_{\rm orb, i}) $ [days]")
    plt.xlim(xmin = 0.03, xmax = 1.02)
    plt.ylim(ymin = 0, ymax = 3.5)
    plt.subplots_adjust(bottom=0.02,left=0.01,right=0.995,top=0.90)
    plt.show()
    #fig.savefig(f'Mdonor{Mdonor}_Z{Zdonor}.jpg',dpi=200)


# %%
###### additional plots

def find_interactions(log):
    """
    Find the formation channel based on the SEVN log
    """
    CE = readlogstring(log,events=["CE"], capturing_cols = ['name','time', 'Phase_CEp', 'Semimajor_preCE','Semimajor_postCE'])   #
    RLOB= readlogstring(log,events=["RLO_BEGIN"], capturing_cols = ['name','time', 'Phase_RLOBd', 'Mass_RLOBd', 'MHE_RLOBd', 'Semimajor_RLOB'])
    RLOE= readlogstring(log,events=["RLO_END"], capturing_cols = ['name','time', 'Phase_RLOEa','Mass_RLOEa', 'MHE_RLOEa','Mlost_RLOE','Maccreted_RLOE', 'Semimajor_RLOB'])

    if len(CE)!=0 :
        perc_a_CE = 100* (CE['Semimajor_postCE'] - CE['Semimajor_preCE']) / CE['Semimajor_preCE']
        perc_a = perc_a_CE[0]
        last_semi = CE['Semimajor_postCE'][0]
        donor_phase = CE['Phase_CEp'][0]
    elif len(RLOE)!= 0:
        if len(RLOB) != len(RLOE): # not all RLOs that start actually end as stable
            RLOB = RLOB.iloc[:-1]  # remove last row, since is the last RLO that could drive e.g. a merger
        
        perc_a_RLO_all = 100* (RLOE['Semimajor_RLOB'] - RLOB['Semimajor_RLOB']) / RLOB['Semimajor_RLOB']
        perc_a_RLO = perc_a_RLO_all.to_numpy()
        idx_max = np.unravel_index(np.argmax(abs(perc_a_RLO), axis=None), perc_a_RLO.shape)  # account for possible multiple mass transfer, picking the one with largest variations
        perc_a = perc_a_RLO[idx_max]
        last_semi = RLOE['Semimajor_RLOB'].to_numpy()[idx_max]
        donor_phase = RLOB['Phase_RLOBd'].to_numpy()[idx_max]
    else:
        last_semi = 0
        perc_a = 0
        donor_phase = -1

    return last_semi, perc_a, donor_phase



def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_fig_semimajor(ql,pl,ll,ap):
    ap_GW = np.concatenate([ap[ll=='BBHmSMT'],ap[ll=='BBHmCE']])
    norm = colors.Normalize(vmin=min(ap_GW), vmax=max(ap_GW))
    cmap = truncate_colormap(plt.get_cmap('gnuplot2'), 0.13, 0.87)

    fig,ax=plt.subplots(1,1,figsize=(5.4,5.8))
    #idx=ll=="MergerRL"
    #plt.scatter(ql[idx],np.log10(pl[idx]),c=ap[idx], cmap=cmap, norm=norm, marker="x",label=r"Merger after SMT")
    #idx=ll=="MergerCE"
    #plt.scatter(ql[idx],np.log10(pl[idx]),c=ap[idx], cmap=cmap, norm=norm,marker="+",label=r"Merger after CE")
    #idx=ll=="wide"
    #plt.scatter(ql[idx],np.log10(pl[idx]),c=ap[idx], cmap=cmap, norm=norm, marker="*",label=r"Wide binary")
    idx=ll=="BBHmSMT"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=ap[idx], cmap=cmap, norm=norm,marker="s",label=r"BBH merger after SMT")
    idx=ll=="BBHmCE"
    ax = plt.scatter(ql[idx],np.log10(pl[idx]),c=ap[idx], cmap=cmap, norm=norm,marker="o",label=r"BBH merger after CE")
    #plt.suptitle(r'$ M_{\rm BH}=$'+ str(m2) +r'$ M_\odot$, Z='+ str(Zdonor))
    fig.legend(ncol=2, bbox_to_anchor=(0.98,1.12),frameon=False)
    plt.xlabel(r"$q_{\rm i} = M_{\rm BH}$/$M_{\rm donor}$")
    plt.ylabel(r"$\log_{10} (P_{\rm orb, i}) $ [days]")
    plt.xlim(xmin = 0.05, xmax = 1)
    plt.ylim(ymin = 0, ymax = 3.5)
    plt.subplots_adjust(bottom=0.02,left=0.01,right=0.995,top=0.90)
    cbar = plt.colorbar(ax, orientation='horizontal', label=r'Post-interaction semi-major axis [$R_\odot$]')
    cbar.set_ticks(np.arange(int(min(ap_GW)),int(max(ap_GW)),2))
    plt.show()
    #fig.savefig(f'Mdonor{Mdonor}_Z{Zdonor}.jpg',dpi=200)


def plot_fig_semimajor_var(ql,pl,ll,av):
    norm = colors.Normalize(vmin=min(av), vmax=100)
    cmap = truncate_colormap(plt.get_cmap('RdBu'), 0.0, 1)

    fig,ax=plt.subplots(1,1,figsize=(5.4,5.8))
    idx=ll=="MergerRL"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=av[idx], cmap=cmap, norm=norm, marker="x",label=r"Merger after SMT")
    idx=ll=="MergerCE"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=av[idx], cmap=cmap, norm=norm,marker="+",label=r"Merger after CE")
    idx=ll=="wide"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=av[idx], cmap=cmap, norm=norm, marker="*",label=r"Wide binary")
    idx=ll=="BBHmSMT"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=av[idx], cmap=cmap, norm=norm,marker="s",label=r"BBH merger after SMT")
    idx=ll=="BBHmCE"
    ax = plt.scatter(ql[idx],np.log10(pl[idx]),c=av[idx], cmap=cmap, norm=norm,marker="o",label=r"BBH merger after CE")
    #plt.suptitle(r'$ M_{\rm BH}=$'+ str(m2) +r'$ M_\odot$, Z='+ str(Zdonor))
    fig.legend(ncol=2, bbox_to_anchor=(0.98,1.12),frameon=False)
    plt.xlabel(r"$q_{\rm i} = M_{\rm BH}$/$M_{\rm donor}$")
    plt.ylabel(r"$\log_{10} (P_{\rm orb, i}) $ [days]")
    plt.xlim(xmin = 0.05, xmax = 1)
    plt.ylim(ymin = 0, ymax = 3.5)
    plt.subplots_adjust(bottom=0.02,left=0.01,right=0.995,top=0.90)
    plt.colorbar(ax, orientation='horizontal', label=r'($a_{\rm post} - a_{\rm pre}) / a_{\rm pre}$ (%)')
    plt.show()
    #fig.savefig(f'Mdonor{Mdonor}_Z{Zdonor}.jpg',dpi=200)



def plot_fig_phase_donor(ql,pl,ll,phd):
    norm = colors.Normalize(vmin=min(phd), vmax=max(phd))
    cmap = truncate_colormap(plt.get_cmap('gnuplot2'), 0.13, 0.87)

    fig,ax=plt.subplots(1,1,figsize=(5.4,5.8))
    idx=ll=="MergerRL"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=phd[idx], cmap=cmap, norm=norm, marker="x",label=r"Merger after SMT")
    idx=ll=="MergerCE"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=phd[idx], cmap=cmap, norm=norm,marker="+",label=r"Merger after CE")
    idx=ll=="wide"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=phd[idx], cmap=cmap, norm=norm, marker="*",label=r"Wide binary")
    idx=ll=="BBHmSMT"
    plt.scatter(ql[idx],np.log10(pl[idx]),c=phd[idx], cmap=cmap, norm=norm,marker="s",label=r"BBH merger after SMT")
    idx=ll=="BBHmCE"
    ax = plt.scatter(ql[idx],np.log10(pl[idx]),c=phd[idx], cmap=cmap, norm=norm,marker="o",label=r"BBH merger after CE")
    #plt.suptitle(r'$ M_{\rm BH}=$'+ str(m2) +r'$ M_\odot$, Z='+ str(Zdonor))
    fig.legend(ncol=2, bbox_to_anchor=(0.98,1.12),frameon=False)
    plt.xlabel(r"$q_{\rm i} = M_{\rm BH}$/$M_{\rm donor}$")
    plt.ylabel(r"$\log_{10} (P_{\rm orb, i}) $ [days]")
    plt.xlim(xmin = 0.05, xmax = 1)
    plt.ylim(ymin = 0, ymax = 3.5)
    plt.subplots_adjust(bottom=0.02,left=0.01,right=0.995,top=0.90)
    plt.colorbar(ax, orientation='horizontal', label='Phase donor at MT start')
    plt.show()
    #fig.savefig(f'Mdonor{Mdonor}_Z{Zdonor}.jpg',dpi=200)



def plot_grid(m1, Zdonor, qgrid,Pgrid, plot='all'):
    ql=[] # To store the mass-ratios
    pl=[] # To store the periods
    ll=[] # To store the label
    ap=[] # To store final semimajor
    av=[] # percentage of semimajor axis variation
    phd=[]  # phase of donor
    for p in Pgrid:
        for q in qgrid:
            m1=m1
            m2=m1*q
            a=calc_a(p/365.25,m1,m2)
            #Evolve a binary initialising the second star as a BH
            df,ldic=sevnwrap.evolve_binary(a,0.,m1,Zdonor,m2,Zdonor,
                                            star_flag_1="BH", snmodel='delayed')
            ll.append(find_label(df,ldic["Log"]))
            ql.append(q)
            pl.append(p)
            a,a_var, ph = find_interactions(ldic["Log"])
            ap.append( a )
            av.append( a_var )
            phd.append(ph)
    SEVNmanager.close()

    ql =  np.array(ql)
    pl =  np.array(pl)
    ll =  np.array(ll)
    ap =  np.array(ap)
    av = np.array(av)
    phd = np.array(phd)

    # plot
    if (plot == 'gg') or (plot=='all'):
        plot_fig(ql,pl,ll,m1,Zdonor)
    
    if (plot == 'mt') or (plot=='all'):
        plot_fig_semimajor(ql,pl,ll,ap)
        plot_fig_semimajor_var(ql,pl,ll,av)
        plot_fig_phase_donor(ql,pl,ll,phd)


# %% [markdown]
# ## Single tracks

# %%
def plot_single(dfP, xmin, xmax):
    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    plt.plot(dfP['Worldtime'],dfP['Semimajor'], color='dodgerblue', label='Semimajor orbit')
    plt.plot(dfP['Worldtime'],dfP['Radius_0'], color='k', linestyle='dashed', label='Radius')
    plt.plot(dfP['Worldtime'],dfP['RL0'], color='r', linestyle='dotted', label='Roche lobe')
    #plt.plot(dfP['Worldtime'],Darwin(dfP), label='Critical semimajor for Darwin instability')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start')
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=0.1)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel('Semimajor [Rsun]')
    plt.show()

    plt.plot(dfP['Worldtime'],dfP['Mass_0'], color='dodgerblue', label='Total')
    plt.plot(dfP['Worldtime'],dfP['MHE_0'], color='r', label='He core')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start')
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.legend()
    plt.xlabel('Time [Myr]')
    plt.ylabel('Stellar mass [Msun]')
    plt.show()

    plt.plot(dfP['Worldtime'],dfP['Inertia_0'], label='Inertia')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start')
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel('Inertia [Msun Rsun^2]')
    plt.show()

    plt.plot(dfP['Worldtime'],dfP['OmegaSpin_0'], label='Star')
    plt.plot(dfP['Worldtime'],2.0*np.pi/ dfP['Period'], label='Orbit')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start')
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel('Angular velocity [1/yr]')
    plt.show()

    plt.plot(dfP['Worldtime'],dfP['AngMomSpin_0'], label='Star')
    plt.plot(dfP['Worldtime'],dfP['AngMom'], label='Orbit')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start')
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel('Angular momentum [Msun Rsun^2 / yr]')
    plt.show()



def calc_a(per,m1,m2):
    """
    From period to semi-major axss
    :param per:  Periods in  year
    :param m1: Mass in Msun
    :param m2:  Mass in Msun
    :return: semi-major axis in Rsun
    """
    G = 3.925125598496094e8
    return (per*per*(G*(m1+m2))/(np.pi*np.pi*4))**(1./3.)


# calculate semimajor for onset of Darwin instability
# according to Hut 1981: instability if a<a_Darwin
# input: dictionary ; output: numpy array
def Darwin(df):
    mass_red = df['Mass_0'] * df['Mass_1'] / (df['Mass_0'] + df['Mass_1']) # reduced mass
    semi_darwin = ( 3*( df['Inertia_0'] + df['Inertia_1' ] ) / mass_red)**0.5
    return semi_darwin

########
###
#######




    #######################################
    
    ############## variations   ###########

    #######################################

def plot_dadt_scatter(dfP, xmin, xmax,ymin,ymax):
    plt.rcParams.update({
    "font.family": "cmr10",   # computer modern font, as in latex
    "font.size": 17,
    "axes.labelsize" : 17,
    "xtick.labelsize" : 17,
    "ytick.labelsize" : 17,
    "legend.fontsize": 15
    })

    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    # Semimajor variations da/dt
    fig = plt.figure(figsize=(10,8))
    plt.scatter(dfP['Worldtime'], dfP['dSemimajordt'], s=100, color='lightskyblue',  label=r'Total (+)')
    plt.scatter(dfP['Worldtime'],-dfP['dSemimajordt'], s=100, color='navy',  label=r'Total (-)')
    plt.scatter(dfP['Worldtime'], dfP['dadtTides'], s=50, color='gold',  label=r'Tides (+)')
    plt.scatter(dfP['Worldtime'],-dfP['dadtTides'], s=50, color='salmon',  label=r'Tides (-)')
    plt.scatter(dfP['Worldtime'],dfP['dadtRLO'], color='forestgreen',  label=r'RLO (+)')
    plt.scatter(dfP['Worldtime'],-dfP['dadtRLO'], color='limegreen',  label=r'RLO (-)')
    plt.scatter(dfP['Worldtime'],-dfP['dadtGW'], color='fuchsia',  label=r'GW (-)')
    plt.scatter(dfP['Worldtime'],dfP['dadtWinds'], color='indigo',  label=r'Winds (+)')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start')
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end')

    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.legend(loc='lower left', ncols=2)
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel(r'$|$da$|$/dt  [R$_\odot$ / Myr ]')
    plt.show()


def plot_dadt_lines(dfP, xmin, xmax,ymin, ymax):
    plt.rcParams.update({
    "font.family": "cmr10",   # computer modern font, as in latex
    "font.size": 17,
    "axes.labelsize" : 17,
    "xtick.labelsize" : 17,
    "ytick.labelsize" : 17,
    "legend.fontsize": 15
    })
        
    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    # Semimajor variations da/dt
    fig = plt.figure(figsize=(10,8))
    plt.plot(dfP['Worldtime'], dfP['dSemimajordt'].mask(dfP['dSemimajordt'] < 1e-6), lw=10, color='lightskyblue',  label=r'Total (+)')
    plt.plot(dfP['Worldtime'],-dfP['dSemimajordt'].mask(dfP['dSemimajordt'] > - 1e-6), lw=10, color='navy',  label=r'Total (-)')
    plt.plot(dfP['Worldtime'], dfP['dadtTides'].mask(dfP['dadtTides'] < 1e-6), lw=5, color='gold', linestyle='dashed', label=r'Tides (+)')
    plt.plot(dfP['Worldtime'],-dfP['dadtTides'].mask(dfP['dadtTides'] > - 1e-6), lw=5, color='salmon', linestyle='dashed', label=r'Tides (-)')
    plt.plot(dfP['Worldtime'],dfP['dadtRLO'].mask(dfP['dadtRLO'] < 1e-6), lw=3, linestyle='dashdot', color='forestgreen',  label=r'RLO (+)')
    plt.plot(dfP['Worldtime'],-dfP['dadtRLO'].mask(dfP['dadtRLO'] > - 1e-6), lw=3, linestyle='dashdot', color='limegreen',  label=r'RLO (-)')
    plt.plot(dfP['Worldtime'],dfP['dadtWinds'], color='indigo', lw=2,  label=r'Winds (+)')
    plt.plot(dfP['Worldtime'],-dfP['dadtGW'].mask(dfP['dadtGW'] > - 1e-6), color='fuchsia',  linestyle='dotted', label=r'GW (-)')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start',zorder=1)
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end',zorder=1)

    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.legend(loc='lower left', ncols=2)
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel(r'$|$da$|$/dt  [R$_\odot$ / Myr ]')
    plt.show()



def da_lines(dfP,xmin,xmax,ymin,ymax):
    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    # Semimajor variations da
    dt = dfP['BTimestep']    # use internal timestep in Myr

    fig = plt.figure(figsize=(10,8))
    plt.plot(dfP['Worldtime'], dfP['dSemimajordt'].mask(dfP['dSemimajordt'] < 1e-6)*dt, lw=10, color='lightskyblue',  label=r'Total (+)')
    plt.plot(dfP['Worldtime'],-dfP['dSemimajordt'].mask(dfP['dSemimajordt'] > - 1e-6)*dt, lw=10, color='navy',  label=r'Total (-)')
    plt.plot(dfP['Worldtime'], dfP['dadtTides'].mask(dfP['dadtTides'] < 1e-6)*dt, lw=5, color='gold', linestyle='dashed', label=r'Tides (+)')
    plt.plot(dfP['Worldtime'],-dfP['dadtTides'].mask(dfP['dadtTides'] > - 1e-6)*dt, lw=5, color='salmon', linestyle='dashed', label=r'Tides (-)')
    plt.plot(dfP['Worldtime'],dfP['dadtRLO'].mask(dfP['dadtRLO'] < 1e-6)*dt, lw=3, linestyle='dashdot', color='forestgreen',  label=r'RLO (+)')
    plt.plot(dfP['Worldtime'],-dfP['dadtRLO'].mask(dfP['dadtRLO'] > - 1e-6)*dt, lw=3, linestyle='dashdot', color='limegreen',  label=r'RLO (-)')
    plt.plot(dfP['Worldtime'],dfP['dadtWinds']*dt, color='indigo', lw=2,  label=r'Winds (+)')
    plt.plot(dfP['Worldtime'],-dfP['dadtGW'].mask(dfP['dadtGW'] > - 1e-6)*dt, color='fuchsia',  linestyle='dotted', label=r'GW (-)')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start',zorder=1)
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end',zorder=1)


    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.legend(loc='lower left', ncols=2)
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel(r'$|$da$|$  [R$_\odot$]')
    plt.show()


def plot_domega_lines(dfP, xmin, xmax, ymin,ymax):
    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    omega_crit = cl.OmegaCrit(dfP)
    omega_orb = cl.OmegaOrb(dfP)

    # Angular velocity
    fig = plt.figure(figsize=(10,8))
    plt.plot(dfP['Worldtime'], omega_orb, lw=2, color='navy', linestyle='dashed',  label=r'$\Omega_{\rm orb}$')
    plt.plot(dfP['Worldtime'], dfP['OmegaSpin_0'], lw=2, color='dodgerblue', linestyle='solid', label=r'$\Omega_{\rm star}$')
    #plt.plot(dfP['Worldtime'], omega_crit, lw=2, color='r', linestyle='dotted', label=r'$\Omega_{\rm star, crit}$')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start',zorder=1)
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end',zorder=1)

    #plt.scatter(dfP['Worldtime'],- (dadt_tides*1e6*dt_timestep + DA) , color='dodgerblue', s=10, label=r'- $\Delta a (RLO+tides)$')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.legend(loc='lower left', ncols=2)
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel(r' Angular velocity  [ 1/ yr ]')
    plt.show()

def plot_domega_scatter(dfP, xmin, xmax, ymin,ymax):
    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    omega_crit = cl.OmegaCrit(dfP)
    omega_orb = cl.OmegaOrb(dfP)

    # Angular velocity
    fig = plt.figure(figsize=(10,8))
    plt.scatter(dfP['Worldtime'], omega_orb, color='navy',  label=r'$\Omega_{\rm orb}$')
    plt.scatter(dfP['Worldtime'], dfP['OmegaSpin_0'], color='dodgerblue',  label=r'$\Omega_{\rm star}$')
    plt.scatter(dfP['Worldtime'], omega_crit,  color='r',  label=r'$\Omega_{\rm star, crit}$')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start',zorder=1)
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end',zorder=1)

    #plt.scatter(dfP['Worldtime'],- (dadt_tides*1e6*dt_timestep + DA) , color='dodgerblue', s=10, label=r'- $\Delta a (RLO+tides)$')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.legend(loc='lower left', ncols=2)
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel(r' Angular velocity  [ 1/ yr ]')
    plt.show()



def plot_time_sync_lines(dfP, xmin, xmax, ymin,ymax):
    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    # Time for tidal synchronization [yr]
    fig = plt.figure(figsize=(10,8))
    plt.plot(dfP['Worldtime'], cl.TimeSync(dfP), lw=2, color='navy',  label=r'Time sync')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start',zorder=1)
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end',zorder=1)

    #plt.scatter(dfP['Worldtime'],- (dadt_tides*1e6*dt_timestep + DA) , color='dodgerblue', s=10, label=r'- $\Delta a (RLO+tides)$')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.legend(loc='lower left', ncols=3)
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel(r'Time [yr]')
    plt.show()


def plot_time_sync_scatter(dfP, xmin, xmax, ymin,ymax):
    RLOstarttime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[0]
    RLOendtime = dfP.loc[dfP['Radius_0'] >= dfP['RL0']].Worldtime.values[-1]

    # Time for tidal synchronization [yr]
    fig = plt.figure(figsize=(10,8))
    plt.scatter(dfP['Worldtime'], cl.TimeSync(dfP), color='navy',  label=r'Time sync')
    plt.axvline(x=RLOstarttime, linestyle='dotted',color='k', label='RLO start',zorder=1)
    plt.axvline(x=RLOendtime, linestyle='dashdot',color='k', label='RLO end',zorder=1)

    #plt.scatter(dfP['Worldtime'],- (dadt_tides*1e6*dt_timestep + DA) , color='dodgerblue', s=10, label=r'- $\Delta a (RLO+tides)$')
    plt.xlim(xmin=xmin,xmax=xmax)
    plt.ylim(ymin=ymin, ymax=ymax)
    plt.legend(loc='lower left', ncols=3)
    plt.yscale('log')
    plt.xlabel('Time [Myr]')
    plt.ylabel(r'Time [yr]')
    plt.show()