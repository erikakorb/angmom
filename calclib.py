# %%
"""
Routines to calculate physical quantities for tidal effects and similar
"""

import numpy as np
import pandas as pd
import astropy.constants as const
import astropy.units as u

########## Angular velocity ###########


#### As in SEVN

# stellar critical angular velocity [1/yr]
# invert spin definition of Spin = omega / omega crit
# to omega crit = omega / spin
def OmegaCrit(df):
    return df['OmegaSpin_0'] / df['Spin_0']     

# orbital angular velocity [1/yr]
def OmegaOrb(df):
    return 2.0*np.pi / df['Period']

# eq 3 of Gagnier 2019
# https://www.aanda.org/articles/aa/full_html/2019/05/aa34599-18/aa34599-18.html#R48
# where
# SEVN uses stellar radius as proxy of the polar one, so the resulting Angular momentum is the minimum possible
# Note: the alternative would be to use it as proxy of the Equatorial radius (Req = R)
def OmegaCritFromRadius(M, R):
    G = const.G.to('R_sun3 / (M_sun yr2)').value  # [R_sun^3 / (M_sun * yr^2)], converted with astropy  
    Req = 1.5 *R
    return np.sqrt(G*M/(Req**3))



## alternative manual calculations

# angular velocity of the orbit [yr^-1]
#
#def omega_orbit(Mbin, a):
#    G = const.G.to('R_sun3 / (M_sun yr2)').value  # [R_sun^3 / (M_sun * yr^2)], converted with astropy  
#    return np.sqrt(G*Mbin) * a**(-1.5)    # fun fact: if -3.5 instead di -1.5, spin-orbit is synchronized
#
#def omega_orbit_from_SEVN(dfP):
#    mu = dfP['Mass_0'] * dfP['Mass_1'] / (dfP['Mass_0'] + dfP['Mass_1'])   # reduced mass
#    return dfP['AngMom'] / (mu * dfP['Semimajor']**2)






###############################

############ Tides ############

###############################



###########################################################################

## manual calculations of angular momentum and semimajor axis variations ##

############################################################################




def calculate_ang_mom_star(M,R, I, spin):
    # J = I * omega
    # [Msun Rsun^2 yr^-1] = [Msun Rsun^2] [yr^-1]
    omega_ang_crit = OmegaCritFromRadius(M, R)
    return I * spin * omega_ang_crit 



def calculate_ang_mom_BH(M, I, spin):
    c = const.c.to('R_sun / (yr)').value   # [Rsun / yr] speed of light
    G = const.G.to('R_sun3 / (M_sun yr2)').value    # gravitational constant
    R_Sch = 2*G*M *c**(-2)   # Schwarzschild radius
    omega_crit = OmegaCritFromRadius(M, R_Sch)
    return I * spin * omega_crit


def calculate_ang_mom_orb(M1, M2, a):
    G = const.G.to('R_sun3 / (M_sun yr2)').value  # [R_sun^3 / (M_sun * yr^2)], converted with astropy  
    Mbin = M1+M2
    return (M1*M2)*G**(0.5) * Mbin**(-0.5) * a**(0.5)   # as in Hut 1980

def calculate_ang_mom_tot(df):
    Jstar = calculate_ang_mom_star(df['Mass_0'],df['Radius_0'], df['Inertia_0'], df['Spin_0'])
    JBH = calculate_ang_mom_BH(df['Mass_1'], df['Inertia_1'], df['Spin_1'])
    Jorb = calculate_ang_mom_orb(df['Mass_0'], df['Mass_1'], df['Semimajor'])
    Jtot = Jorb + Jstar + JBH 
    Jcrit = 3* (Jstar + JBH)
    return Jorb, Jstar, JBH,  Jtot,  Jcrit

def ang_velocities(df):
    omega_star = df['Spin_0'] * OmegaCritFromRadius(df['Mass_0'], df['Radius_0']) 
    omega_BH = calculate_ang_mom_BH(df['Mass_1'], df['Inertia_1'], df['Spin_1']) / df['Inertia_1']
    return omega_star, omega_BH




#### Get quantities from SEVN


# effective inertia
# adopted only in SEVN
# and only after the sync timescale, thus the possible effective, timestep is calculated
# but it is the one used in the calculation of star angular momentum variation
# thus, in the orbital angular momentum variation
def Inertia_eff(df):
    Reff = df['Radius_eff']
    Rstar = df['Radius_0']
    return df['Inertia_0']* Reff*Reff/(Rstar*Rstar)




# calculate k_tides of Zahn+1977
# according to the routine "Tides_simple" in /sevn/src/binstar/orbit.cpp 
# Note: ONly SEVN calcluates it with the Effective Radius

def kt_rad(df):
    M = df['Mass_0']                     # assume primary as donor
    q = df['Mass_1'] / M
    R = df['Radius_eff']                # effective radius, is typical of SEVN!
    a = df['Semimajor']
    G = const.G.to('R_sun3 / (M_sun yr2)').value    # gravitational constant
    
    E2 = 1.58313e-9 * M**(2.84)   
    kt = np.sqrt(G * M * R*R / (a*a*a*a*a)) * (1+q)**(0.833333) * E2   # eq 42 Hurley+02 with typos fixed as in eq 37 of Iorio+23
    return kt

def kt_conv(df):
    M_env_conv = df['Qconv_0'] * df['Mass_0']   # convective envelope mass
    M = df['Mass_0']                     # assume primary as donor
    R = df['Radius_eff']                # effective radius, is typical of SEVN!

    tau_conv = df['Tconv_0']   # convective cells turnover time in yr
    omega_orb = 2*np.pi / df['Period']   # orbital angular velocity in 1/yr
    omega_spin = df['OmegaSpin_0']            # spin of the primary donor star
    P_tid = 2*np.pi/abs(omega_orb - omega_spin)  # equivalent to Eq. 33 of Hurley+2002
    f_conv = min(1, 0.25 * P_tid*P_tid / (tau_conv*tau_conv))
    kt = 0.095 * f_conv * M_env_conv / (tau_conv * M)    # eq 30 in Hurley+02 [1/yr]
    return kt

def k_tides(df):
    df['kt'] =  df.apply(lambda x: kt_rad(x) if (x['Qconv_0'] == 0) else kt_conv(x), axis=1)
    return df['kt']

# Hut's polynomials of eccentricity (Eq. 11 of Hut 1981)
def f_pol(df):
    e = df['Eccentricity']
    f1 = 1 + 14* e*e + 26.25 * e**4  + 8.75 * e**6 + 0.27 * e**8
    f2 = 1 + 7.5* e*e + 5.625 * e**4 + 0.31 * e**6
    f3 = 1 + 3.75* e*e + 1.875 * e**4 +  0.078 * e**6
    f4 = 1 + 1.5* e*e + 0.125 * e**4
    f5 = 1 + 3* e*e + 0.375 * e**4
    return f1,f2,f3,f4,f5 

# Equilibrium angular velocity 
# Where no angular momentum can be trasferred
# It is an estimate of the angular velocity where synchronization and circularization is reached
# Eq 34 of Hurley+2002
def OmegaEq(df):
    f1,f2,f3,f4,f5 = f_pol(df)           # Hut polynomials
    omega_orb = OmegaOrb(df)             # orbital angular velocity in 1/yr
    return f2*omega_orb / (f5*(1-df['Eccentricity']*df['Eccentricity'])**(1.5))

# time required to reach spin-orbit synchronization
# i.e. stellar angular velocity = OmegaEq calculated as in Eq.34 of Hurley+2002
# if time sync < timestep, time_sync is used as new effective timestep
# here it is implemented as in the SEVN routine inside Tides_simple
def TimeSync(df, Reff=False):
    return (OmegaEq(df) - df['OmegaSpin_0'])/ domegastar_tides(df, Reff)   # in yr





# variation of semi-major axis due to tides
# Eq. 9 in Hut 1981   and/or  eq. 34 in Iorio+2023  
# dadt_hut routine in /sevn/src/binstar/Orbit.cpp
def delta_a_tides(df):
    e = df['Eccentricity']
    M = df['Mass_0']                    # assume primary as donor and star affecged by tides
    q = df['Mass_1'] / M                # technically mass ratio of perturbing star / star affected by tides
    R = df['Radius_eff']                # effective radius, is typical of SEVN!
    a = df['Semimajor']
    R_over_a = R/a
    omega_orb = 2.0*np.pi / df['Period']   # orbital angular velocity in 1/yr
    omega_spin = df['OmegaSpin_0']          # spin of the primary donor star
    f1,f2,f3,f4,f5 = f_pol(df)
    kt = k_tides(df)
    dadt = -6.0 *kt*q*(q+1)*(R_over_a**8) * a *(1-e*e)**(-7.5) * (f1 - ((1-e*e)**(1.5)*f2* (omega_spin / omega_orb)))
    return dadt   # [Rsun / yr]

# variation of stellar angular velocity due to tides [1/yr^2]
# Eq. 11 in Hut 1981   and/or  eq. 36 in Iorio+2023  
# dspindt_hut routine in /sevn/src/binstar/Orbit.cpp
def domegastar_tides(df, Reff=False):
    e = df['Eccentricity']
    M = df['Mass_0']                    # assume primary as donor and star affected by tides
    q = df['Mass_1'] / M                # technically mass ratio of perturbing star / star affected by tides
    R = df['Radius_eff']                # effective radius, is typical of SEVN!
    a = df['Semimajor']
    R_over_a = R/a
    omega_orb = OmegaOrb(df)             # orbital angular velocity in 1/yr
    omega_spin = df['OmegaSpin_0']       # angular velocity of the primary donor star
    f1,f2,f3,f4,f5 = f_pol(df)
    kt = k_tides(df)
    if Reff==False:        # default case in SEVN i.e. true stellar radius is used (NB: rg was NOT corrected in the 7/7/25 commit)
        Rstar = df['Radius_0']              # real stellar radius of primary
    else:
        Rstar = df['Radius_eff']             # use effective radius if fixtides branch after 7th JUly 2025 commits(use_Reff=True)
    rg2 = df['Inertia_0']/ (M*Rstar*Rstar)  # I/(M*R**2) scaling for Inertia that uses real stellar radius because of inertia
    domegastar = 3.0 *kt*q*q/rg2 *(R_over_a**6) * omega_orb *(1-e*e)**(-6) * (f2 - ((1-e*e)**(1.5)*f5* (omega_spin / omega_orb)))   # /rg2 , as in the SEVN routine (possible bug though)
    return domegastar       






# variations of angular momentum during RLO

# Mass losses between two output timesteps [Msun in timestep of Myr]
# Due to stellar or binary processes
# Note: beware that SEVN applies separate calculations/outputs for separate stellar or binary processes
def MassVariations(df):
    # note on timestep:
    # 'Timestep' and 'BTtimestep' should be used, as in SEVN
    # but here, since I use outputs also for plotting, is better to keep consistency
    dt_timestep = df['Worldtime'].diff().shift(-1)  # timestep used, in Myr
    DMdwind = df['dMdt_0']*dt_timestep  # mass lost by donor by winds alone [Msun]
    DMawind = df['dMaccwinddt_1']*dt_timestep  # mass accreted by accretor by winds alone [Msun]
    DMdRLO = df['dMRLOdt_0']*dt_timestep  # mass lost by donor by RLO alone [Msun]
    DMaRLO = df['dMRLOdt_1']*dt_timestep  # mass accreted by accretor by RLO alone [Msun]
    DMlostbin = abs(abs(DMdRLO)-abs(DMaRLO))  # mass lost by the system as effect of RLO alone [Msun]
    return DMdwind, DMawind, DMdRLO, DMaRLO, DMlostbin  


# Jorb variation due to wind accretion [Msun Rsun^2 yr^-1] between two timesteps
# i.e. net effect of mass lost/accreted via winds only (Windaccretion.cpp module)
# it impacts semimajor as eq. 19 of Iorio+2023
# and derives from Hurley+02 eq. 17
def Jorb_variation_winds(df):
    DMdwind, DMawind, DMdRLO, DMaRLO, DMlostbin = MassVariations(df)
    Mdonor = df['Mass_0']
    Macc = df['Mass_1']
    Mbin = Mdonor + Macc
    mu_red = Mdonor*Macc/Mbin  # reduced mass
    omega_orb = OmegaOrb(df)
    a = df['Semimajor']
    Jorb = mu_red * a*a*omega_orb
    DeltaJorbWind = Jorb*(DMdwind*Macc - DMawind*Mdonor)/(Mdonor*Mbin)
    return DeltaJorbWind





# gamma factor for angular momentum lost by binary
# in SEVN there are 3 options
# as described in appendix A4.2
# following the definition of Soberman+1997
# 1) Jeans mode i.e. isotropic emission from primary
#       - rlo_gamma_angmom = -1 
# 2) Isotropic re-emission from secondary
#        - rlo_gamma_angomo = -2
# 3) circumbinary disk
#        - rlo_gamma_angomom = fraction oforbital angular momentum lost that deposits in disk
def gamma_RLO(df, rlo_gamma_angmom):
    Md = df['Mass_0']     # donor mass
    Ma = df['Mass_1']
    Mbin = Ma+Md
    if rlo_gamma_angmom == '-1':
        gammaRLO = Ma*Ma /(Mbin*Mbin)
    elif rlo_gamma_angmom == '-2':   # isotropic re-emission from secondary
        gammaRLO = Md*Md / (Mbin*Mbin)
    elif rlo_gamma_angmom >=0:
        gammaRLO = rlo_gamma_angmom
    return gammaRLO

# accretion radius
# following eq. 30 of Iorio+2023
def R_acc(df):
    q2 = 1/ df['Mass_0']/df['Mass_1']  # here it uses the donor-to-accretor ratio!!!
    Rmin = 0.0425 * (q2 + q2*q2)**0.25 * df['Semimajor']  # min radius of mass stream
    Ra = df['Radius_1']     # accretor radius
    Racc = np.select([(Rmin > Ra), (Rmin <= Ra)], 
                            [Ra, 1.7*Rmin],     # possible typo in Iorio+2023
                            default=np.nan) 
    #if Rmin > Ra:
    #    Racc = Ra     # mass accreted from accretor radius (possible typo in Iorio+2023)
    #else:
    #    Racc = 1.7 * Rmin # mass accreted from where the disk would have formed
    return Racc


# J orb variations during RLO
# variations in units of timestep in Myr
# Sez. 2.3.2.3 in Iorio+2023
# Eq. 28-31
def Jorb_variations_RLO(df, rlo_gamma_angmom):
    G = const.G.to('R_sun3 / (M_sun yr2)').value    # gravitational constant
    DMdwind, DMawind, DMdRLO, DMaRLO, DMlostbin = MassVariations(df)           # Msun in timestep in Myr
    gammaRLO = gamma_RLO(df, rlo_gamma_angmom)
    omega_orb = 2.0*np.pi/df['Period']              # 1/yr
    a = df['Semimajor']
    e = df['Eccentricity']
    RL = df['RL0']    # in SEVN, during RLO, stellar radius = Roche lobe radius (effective radius)

    ### note: 
    # In Iorio+2023, all contributions are calculated with their correct sign, then added
    # i.e. Jorblost (-), Jorb_donor (+), Jorb_acc (-)
    # as from the point of the gain/loss of the orbit
    # In routines, DM donor is absolute value, so contributions are again
    # intrinsic                 Jorblost (+)     Jorbdonor (+)  Jorbacc(+)
    # and increase/decrease     Djorblost (+)    Djorndonor (-) Jorbacc(+)
    # then everything is subtracted to assign the correct sign

    # Orbital angular momentum removed by the mass lost by the system, eq 27
    DJorb_lost = - DMlostbin*gammaRLO*a*a*np.sqrt(1-e*e)*omega_orb  
    
    # Orbital angular momentum transferred from donor to orbit, as donor loses mass, eq 28
    DJorb_donor = - DMdRLO*RL*RL*df['OmegaSpin_0']

    # Orbital angular momentum removed by mass accreted , eq 29
    DJorb_acc  = - DMaRLO*np.sqrt(G*df['Mass_1']*R_acc(df))

    # total J orb variation during RLO
    DJ_orb_RLO = DJorb_lost + DJorb_donor + DJorb_acc  

    # new quantities in next timestep because of RLO events
    Jorb_new = df['AngMom'] + DJ_orb_RLO      # Orbital Ang mom
    Md_new = df['Mass_0'] + DMdRLO               # Donor mass (DMd intrinsically negative)
    Ma_new = df['Mass_1'] + DMaRLO               # Accretor mass
    DA = Jorb_new*Jorb_new*(Md_new+Ma_new) / (G*(1-e*e)*Md_new*Md_new*Ma_new*Ma_new) - a   # semimajor variation
    return DJorb_lost, DJorb_donor, DJorb_acc,Jorb_new, DA   #[dquantity in timestep in Myr]


# Orbital variation due to GW emission
# As in Peters 64, implemented in SEVN according routine "Peters_gw" and eq.39 of Iorio+2023
# It could be relevant for tight systems
def DADT_GW(df):
    c = const.c.to('R_sun / (yr)').value   # [Rsun / yr] speed of light
    G = const.G.to('R_sun3 / (M_sun yr2)').value    # gravitational constant
    e = df['Eccentricity']
    a = df['Semimajor']
    M1 = df['Mass_0']
    M2 = df['Mass_1']
    k = 12.8* G*G*G/(c*c*c*c*c)
    fe = 1 + 3.04166 *e*e + 0.3854 * e*e*e*e
    ecc2 = (1-e*e)**(3.5)
    dadt = - k*fe* M1*M2*(M1+M2)/(a*a*a*ecc2)
    return dadt     # [Rsun/yr]
