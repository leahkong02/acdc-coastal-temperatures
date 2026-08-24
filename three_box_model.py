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

def make_forcing(Nyears, T_mean=290):

	#### RETURN F_forcing, T_forcing, Q_forcing

	days_per_year = 365

	Nmons = Nyears*12

	N = Nyears*days_per_year*steps_per_day

	Time = np.arange(Nmons)
	Time_model = np.linspace(0,Nmons-1,N)

	Fmean = 240
	Famp = 0
	rand_F = 30

	F_cyc = -Famp*np.cos(2*np.pi*Time/12)

	T_amp = 0
	rand_T = 0
	T_cyc = -T_amp*np.cos(2*np.pi*Time/12 - np.pi/6) # keep fixed and just see what happens for different forcings
	# do perpetual summer for forcing instead (no winter ! )

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

#	plt.figure(figsize=(15,6))
#	plt.subplot(3,1,1)
#	plt.plot(np.arange(120),F_forcing[:120],'k')
#	plt.subplot(3,1,2)
#	plt.plot(np.arange(120),T_forcing[:120],'k')
#	plt.subplot(3,1,3)
#	plt.plot(np.arange(120),RH_forcing[:120],'k')
#	plt.plot(np.arange(120),qs[:120])
#	plt.plot(np.arange(120),q_forcing[:120],'k')
#	plt.savefig('SI_Forcing.pdf')
#	plt.show()
#	f = breakhere

################################################################

	#### GETTING RARE CASES OF NEGATIVE NET DOWNWARD SOLAR
	np.putmask(F_forcing,F_forcing<0,0)

	F_forcing = np.interp(Time_model,Time,F_forcing)
	q_forcing = np.interp(Time_model,Time,q_forcing)
	T_forcing = np.interp(Time_model,Time,T_forcing)

	return(F_forcing,T_forcing,q_forcing)

def the_model(F,To,q_R):

##################################################################################
###################### SIMPLE MODEL FOR FIGURE 3 #########################
##################################################################################


################################ PHYSICAL CONSTANTS ###############################

	N = len(F)	 		# Number of soil moisture values we simulate
	T_freeze = 273.15	# Kelvin
	P_s	 = 1000		# hPa

######## PARAMETERS

	emis = 0.67 	# atmospheric emissivity
	r_lo = 250
	r_so = 1000
	r_ls = 250
	r_ss = 1000
	we_l = 0.005
	we_o = 0.005
	th_ft = 290
	q_ft = 0.003 #calc_q_s(th_ft, 900)*0.5
	tau = 86400*2 #stochasticity would be very fun or calculate statistically 
	h_s = 0.1
	h_bl = 1000
	theta = 0.45

############### Physical Constants
	
	sig = 5.67e-08			# SB
	rho_a = 1.25			# density of air [kg/m^3]
	rho_l = 1000			# denisty of water [kg/m^3]
	rho_s = 1000			# density of dry soil [kg/m^3]
	c_pl  = 1000			# heat capacity of dry soil [J/kg/K]
	c_po  =  4182			# heat capacity of water [J/kg/K]
	c_pa = 1003
	L = 2257000				# Latent enthalpy of vaporization [J/kg]
	mu_s = rho_l*h_s

################# State Variables

	Ts = np.zeros(N)
	m_s = np.zeros(N)
	th_o = np.zeros(N)
	th_l = np.zeros(N)
	q_o = np.zeros(N)
	q_l = np.zeros(N)

############# Functions of state variables

	LHF_o = np.zeros(N)
	LHF_l = np.zeros(N)
	P = np.zeros(N)

	Ts[0] = 298			# Initial condition for surface temperature
	m_s[0] = 0.01			# ""			surface moisture
	th_l[0] = 298
	th_o[0] = 298
	q_l[0] = 0.005
	q_o[0] = 0.008
############## PRECIP STUFF ###################################################

	P_avg = 8				# average precipitation intensity [mm] precip very different across west and east coasts
	a_1 = .2				# precip frequency [1/days]
	omega = a_1/steps_per_day

####################################################################################

	sec_per_day = 86400		# seconds per day
	dt = 86400./steps_per_day		# time increment
	i = 0

	while i < N-1:

		I_l = sig*(Ts[i]**4) - sig*emis*(th_l[i]**4)
		I_o = sig*(To[i]**4) - sig*emis*(th_o[i]**4)
		SHF_l = c_pa*rho_a*(Ts[i] - th_l[i])/r_ss
		SHF_o = c_pa*rho_a*(To[i] - th_o[i])/r_so

		es_sat_l = e_s(Ts[i])
		qs_sat_l = es_sat_l*0.622/(P_s - 0.37*es_sat_l)
		q_diff_l = qs_sat_l - q_l[i]			# Humidity difference near surface

		es_sat_o = e_s(To[i])
		qs_sat_o = es_sat_o*0.622/(P_s - 0.37*es_sat_o)
		q_diff_o = qs_sat_o - q_o[i]
		
		if q_diff_l > 0:
			E_s_l = rho_a*m_s[i]*q_diff_l/(r_ls)		# Surface Transpiration - no theta?
		else:
			E_s_l = 0

		if q_diff_o > 0:
			E_s_o = rho_a*q_diff_o/(r_lo)
		else:
			E_s_o = 0

		does_rain = rand.rand()
		if does_rain < omega:
			P[i] = rand.gamma(P_avg,scale=1)/1000 # m rainfall (intensity)
		else:
			P[i] = 0
	
		LHF_o[i] = L*E_s_o		# Latent heat flux over ocean
		LHF_l[i] = L*E_s_l
		C_eff = h_s*(c_pl*rho_s + c_po*m_s[i]*rho_l)	# Effective heat capacity of storage [J/m^2/K]

		dTs_dt = (F[i] - I_l - LHF_l[i] - SHF_l)/C_eff
		dms_dt = P[i]/(h_s) - E_s_l/mu_s
		dthl_dt = (th_o[i] - th_l[i])/tau + (SHF_l + I_l - sig*emis*(th_l[i]**4) + sig*emis*emis*th_ft**4)/(c_pa*rho_a*h_bl) + (th_ft - th_l[i])*we_l/h_bl
		dtho_dt = -(th_o[i] - th_l[i])/tau + (SHF_o + I_o - sig*emis*(th_o[i]**4) + sig*emis*emis*th_ft**4)/(c_pa*rho_a*h_bl) + (th_ft - th_o[i])*we_o/h_bl
		dql_dt = (q_o[i] - q_l[i])/tau + E_s_l/(h_bl*rho_a) + (q_ft - q_l[i])*we_l/h_bl
		dqo_dt = -(q_o[i] - q_l[i])/tau + E_s_o/(h_bl*rho_a) + (q_ft - q_o[i])*we_o/h_bl

		Ts[i+1] 	= Ts[i] + dTs_dt*dt
		m_s[i+1]	= min(m_s[i] + dms_dt*dt, 1) # cap soil moisture at 1?? 
		th_l[i+1]	= th_l[i] + dthl_dt*dt
		th_o[i+1]	= th_o[i] + dtho_dt*dt
		q_l[i+1]	= q_l[i] + dql_dt*dt
		q_o[i+1]	= q_o[i] + dqo_dt*dt
	
		i+=1

	return(Ts,m_s,th_l,th_o,q_l,q_o)
	
F,T,q = make_forcing(100, T_mean=290)
Ts,m_s,th_l,th_o,q_l,q_o = the_model(F,T,q)

Ts_C = [x-273 for x in Ts]
# plt.plot(np.arange(len(Ts)),m_s)
plt.hist(Ts_C, bins = 40)
plt.xlabel('land surface temp in C')
plt.suptitle(f'mean = {np.mean(Ts_C)},skew = {skew(Ts_C)}, oceantemp = {T[0]-273}')
plt.show()
# plt.savefig('lucas_code.png')
