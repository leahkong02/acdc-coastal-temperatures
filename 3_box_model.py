import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
from netCDF4 import Dataset
from scipy.stats import pearsonr, skew

global steps_per_day
steps_per_day = 20

def e_s(T):
    TC = T - 273.15
    # Function to calculate saturation mixing ratio. Returns pressure in hPa
    e_s = 6.11*10**(7.5*TC/(237.5+TC))  # TEMPERATURES IN CELSIUS!!!
    return(e_s)

def calc_q_s(T,P):
    # Function to calculate saturation mixing ratio as a function of temperature (in celsius) and pressure in hPa
    es = e_s(T)
    return(es*0.622/(P - 0.37*es))

def make_forcing(Nyears):

    #### RETURN F_forcing, T_forcing, Q_forcing

    days_per_year = 365

    Nmons = Nyears*12

    N = Nyears*days_per_year*steps_per_day

    Time = np.arange(Nmons)
    Time_model = np.linspace(0,Nmons-1,N)

    Fmean = 240
    Famp = 0
    #Famp = 80
    rand_F = 30

    F_cyc = -Famp*np.cos(2*np.pi*Time/12)
    #F_cyc = 0

    T_mean = 290
    T_amp = 0
    rand_T = 0
    T_cyc = -T_amp*np.cos(2*np.pi*Time/12 - np.pi/6)
    

    q_mean = 0.003
    q_amp  = 0
    rand_q = 0

    q_cyc = -q_amp*np.cos(2*np.pi*Time/12 - np.pi/6)

    i = 0
    qnoise = np.zeros(Nmons)
    Fnoise = np.zeros(Nmons)
    Tnoise = np.zeros(Nmons)

    qwalk = rand.normal(0,rand_T,size=Nmons)*rand_q
    Fwalk = rand.normal(0,rand_F,size=Nmons)
    Twalk = rand.normal(0,rand_T,size=Nmons)

    r = 0.15
    while i < Nmons-1:
        qnoise[i+1] = r*qnoise[i] + qwalk[i]
        Fnoise[i+1] = r*Fnoise[i] + Fwalk[i]
        Tnoise[i+1] = r*Tnoise[i] + Twalk[i]

        i+=1

    q_forcing = q_mean + q_cyc + qnoise
    F_forcing = Fmean + F_cyc + Fnoise
    T_forcing = T_mean + T_cyc + Tnoise

####################### Making a supplemental Figure

#    plt.figure(figsize=(15,6))
#    plt.subplot(3,1,1)
#    plt.plot(np.arange(120),F_forcing[:120],'k')
#    plt.subplot(3,1,2)
#    plt.plot(np.arange(120),T_forcing[:120],'k')
#    plt.subplot(3,1,3)
#    plt.plot(np.arange(120),RH_forcing[:120],'k')
#    plt.plot(np.arange(120),qs[:120])
#    plt.plot(np.arange(120),q_forcing[:120],'k')
#    plt.savefig('SI_Forcing.pdf')
#    plt.show()
#    f = breakhere

################################################################

    #### GETTING RARE CASES OF NEGATIVE NET DOWNWARD SOLAR
    np.putmask(F_forcing,F_forcing<0,0)

    F_forcing = np.interp(Time_model,Time,F_forcing)
    q_forcing = np.interp(Time_model,Time,q_forcing)
    T_forcing = np.interp(Time_model,Time,T_forcing)

    return(F_forcing,T_forcing,q_forcing)



def the_model(F, T_ocean):

    ##################################################################################
###################### SIMPLE MODEL FOR LAND-OCEAN #########################
##################################################################################

    #from scipy.stats import pearsonr, skew

    ################################ PHYSICAL CONSTANTS ###############################

    N = len(F)

    T_freeze = 273.15              # Kelvin
    P_s = 1000                     # hPa

    ######## PARAMETERS

    emis = 0.67                    # atmospheric emissivity
    r_ls = 250                    # land latent resistance [s/m]
    r_ss = 1000                    # land sensible resistance [s/m]
    r_lo = 250                    # ocean latent resistance [s/m]
    r_so = 1000                    # ocean sensible resistance [s/m]

    ######## BOUNDARY LAYER PARAMETERS

    h_o = 1000                     # ocean boundary layer depth [m]
    h_l = 1000                     # land boundary layer depth [m]

    #tau_o_ft = 86400               # ocean BL - FT mixing timescale [s]
    #tau_l_ft = 86400               # land BL - FT mixing timescale [s]
    tau_o_l = 86400*2               # ocean - land exchange timescale [s]
    w_e = 0.005                    # entrainment rate [m/s]
    q_ft = 0.003                   # free troposphere specific humidity

    ######## LAND SURFACE

    theta = 1                   # maximum volumetric soil moisture
    h_s = 0.1                      # land surface layer depth [m]

    ############### Physical Constants

    sig = 5.67e-08                 # Stefan-Boltzmann constant
    rho_a = 1.25                   # density of air [kg/m3]
    rho_l = 1000                   # density of water [kg/m3]
    rho_s = 1000                   # density of dry soil [kg/m3]
    c_ps = 1000                    # heat capacity dry soil [J/kg/K]
    c_pa = 1003                    # heat capacity air [J/kg/K]
    c_pl = 4182                    # heat capacity water [J/kg/K]
    L = 2257000                    # latent heat of vaporization [J/kg]

    ################### COMBINED PARAMETERS

    mu_s = rho_l*h_s             # water storage capacity [kg/m2]
    
    ############# INITIAL PARAMS
    
    theta_ocean_0 = 298
    theta_land_0 = 298
    theta_ft_0 = 290
    q_ocean_0 = 0.008
    q_land_0 = 0.005

    ################# STATE VARIABLES

    Ts = np.zeros(N)               # land surface temperature [K]
    m_s = np.zeros(N)              # land surface moisture
    theta_ocean = np.zeros(N)      # ocean BL potential temperature [K]
    theta_land = np.zeros(N)       # land BL potential temperature [K]
    theta_ft = np.zeros(N)         # free troposphere potential temperature [K]
    q_ocean = np.zeros(N)          # ocean atm specific humidity
    q_land = np.zeros(N)           # land atm specific humidity

    ################# FLUX VARIABLES

    H_land = np.zeros(N)           # land sensible heat flux [W/m2]
    H_ocean = np.zeros(N)          # ocean sensible heat flux [W/m2]

    LE_land = np.zeros(N)          # land latent heat flux [W/m2]
    LE_ocean = np.zeros(N)         # ocean latent heat flux [W/m2]

    E_land = np.zeros(N)           # land evaporation [kg/m2/s]
    E_ocean = np.zeros(N)          # ocean evaporation [kg/m2/s]

    LHF = np.zeros(N)              # total latent heat flux [W/m2]

    P = np.zeros(N)                # precipitation [m]

    ################# INITIAL CONDITIONS

    Ts[0] = theta_land_0

    m_s[0] = 0.01

    theta_ocean[0] = theta_ocean_0
    theta_land[0] = theta_land_0
    theta_ft[0] = theta_ft_0

    q_ocean[0] = q_ocean_0
    q_land[0] = q_land_0
    ############## PRECIPITATION

    P_avg = 8                       # average precipitation intensity [mm]
    a_1 = 0.2                      # precip frequency [1/days]
    omega = a_1 / steps_per_day

    ###########  TIME STEP

    sec_per_day = 86400         # seconds per day
    dt = 86400. / steps_per_day         # time increment (10 chunks per day)
    i = 0

    while i < N-1:

        ########### LONGWAVE RADIATION FOR LAND SURFACE
        OLR = sig*(Ts[i]**4)
        DLR = sig*emis*(theta_land[i]**4)
        
        ########### LONGWAVE RADIATION FOR ATM OVER LAND
        LR_l = sig*emis*(Ts[i]**4) + sig*emis*(theta_ft[i]**4) - 2*sig*emis*(theta_land[i]**4)
        
        ########### LONGWAVE RADIATION FOR ATM OVER OCEAN
        LR_o = sig*emis*(T_ocean[i]**4) + sig*emis*(theta_ft[i]**4) - 2*sig*emis*(theta_ocean[i]**4)

        ########## LAND SENSIBLE HEAT FLUX

        H_land[i] = rho_a*c_pa*(Ts[i] - theta_land[i])/r_ss

        ########### OCEAN SENSIBLE HEAT FLUX

        H_ocean[i] = rho_a*c_pa*(T_ocean[i] - theta_ocean[i])/r_so

        ########### LAND SURFACE HUMIDITY

        es_land = e_s(Ts[i])
        qs_land = (es_land*0.622/(P_s - 0.37*es_land))
        q_diff_land = qs_land - q_land[i]

        ###########  LAND EVAPORATION

        if q_diff_land > 0:
            E_land[i] = rho_a*m_s[i]*q_diff_land/(theta*r_ls)
        else:
            E_land[i] = 0

        ########### LAND LATENT HEAT FLUX

        LE_land[i] = L*E_land[i]

        ###########  OCEAN SURFACE HUMIDITY

        es_ocean = e_s(T_ocean[i])
        qs_ocean = (es_ocean*0.622/(P_s - 0.37*es_ocean))
        q_diff_ocean = qs_ocean - q_ocean[i]

        ###########  OCEAN EVAPORATION

        if q_diff_ocean > 0:
            E_ocean[i] = rho_a*q_diff_ocean/r_lo
        else:
            E_ocean[i] = 0

        ############  OCEAN LATENT HEAT FLUX

        LE_ocean[i] = L*E_ocean[i]

        ############  PRECIPITATION

        does_rain = rand.rand()
        if does_rain < omega:
            P[i] = rand.gamma(P_avg,scale=1)/1000
        else:
            P[i] = 0

        #############  TOTAL LATENT HEAT FLUX

        LHF[i] = LE_land[i] + LE_ocean[i]

        #############  LAND SURFACE HEAT CAPACITY

        C_eff = h_s*(c_ps*rho_s + c_pl*m_s[i]*rho_l)

        ############# LAND SURFACE TEMPERATURE

        dTs_dt = (F[i] - OLR + DLR - LE_land[i] - H_land[i])/C_eff

        #############  LAND SOIL MOISTURE

        dm_s_dt = P[i]/(h_s*theta) - E_land[i]/mu_s

        #############  OCEAN BL TEMPERATURE

        dtheta_ocean_dt = ((

            # Ocean sensible heating
            H_ocean[i]
            
            # Longwave radiation
            + LR_o
            
            ) / (rho_a*c_pa*h_o)
            # Mixing with free troposphere
            - (theta_ocean[i] - theta_ft[i])*w_e/h_o

            # Ocean-land exchange
            - (theta_ocean[i] - theta_land[i])/tau_o_l
        )

        ##############  LAND BL TEMPERATURE

        dtheta_land_dt = ((

            # Sensible heat from land
            H_land[i]
            
            # Longwave radiation
            + LR_l
            
            ) / (rho_a*c_pa*h_l)
            # Mixing with free troposphere
            - (theta_land[i] - theta_ft[i])*w_e/h_l

            # Ocean-land exchange
            + (theta_ocean[i] - theta_land[i])/tau_o_l
        )

        ############## FREE TROPOSPHERE TEMPERATURE

        dtheta_ft_dt = 0
        
        ########### OCEAN SPECIFIC HUMIDITY DYNAMICS
        
        dqo_dt = (q_land[i] - q_ocean[i])/tau_o_l + (E_ocean[i]/(h_o*rho_a)) + (q_ft - q_ocean[i])*w_e/h_o
        
        ########### LAND SPECIFIC HUMIDITY DYNAMICS
        
        dql_dt = (q_ocean[i] - q_land[i])/tau_o_l + (E_land[i]/(h_l*rho_a)) + (q_ft - q_land[i])*w_e/h_l

        ############  UPDATE STATE VARIABLES

        Ts[i+1] = Ts[i] + dTs_dt*dt

        m_s[i+1] = m_s[i] + dm_s_dt*dt

        theta_ocean[i+1] = theta_ocean[i] + dtheta_ocean_dt*dt

        theta_land[i+1] = theta_land[i] + dtheta_land_dt*dt

        theta_ft[i+1] = theta_ft[i] + dtheta_ft_dt*dt
        
        q_ocean[i+1] = q_ocean[i] + dqo_dt*dt
        
        q_land[i+1] = q_land[i] + dql_dt*dt

        ############  MOISTURE LIMITS

        if m_s[i+1] > theta:
            m_s[i+1] = theta
        if m_s[i+1] < 0:
            m_s[i+1] = 0

        i += 1

    return (
        Ts,
        m_s,
        theta_ocean,
        theta_land,
        P,
        q_ocean,
        q_land
    )


F,T,q = make_forcing(100)
Ts,ms,thetao,thetal,P,qo, ql = the_model(F, T)

#plt.hist(Ts)
#plt.savefig('ACDC_hist.png')
#plt.close()


Ts_C = [x-273 for x in Ts]
# plt.plot(np.arange(len(Ts)),Ts)
plt.hist(Ts_C, bins = 40)
#plt.hist(ms, bins=40)
plt.xlabel('land surface temp in C')
plt.suptitle(f'mean = {np.mean(Ts_C)},skew = {skew(Ts_C)}, oceantemp = {T[0]}')
plt.savefig('ACDC_hist.png')
plt.close()
