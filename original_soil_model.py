import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
from netCDF4 import Dataset

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

	Fmean = 160
	Famp = 80
	rand_F = 30

	F_cyc = -Famp*np.cos(2*np.pi*Time/12)

	T_mean = 290
	T_amp = 4
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

def the_model(F,T_R,q_R):

##################################################################################
###################### SIMPLE MODEL FOR FIGURE 3 #########################
##################################################################################

	from scipy.stats import pearsonr

################################ PHYSICAL CONSTANTS ###############################

	N = len(F)	 		# Number of soil moisture values we simulate
	T_freeze = 273.15	# Kelvin
	P_s	 = 1000		# hPa

######## PARAMETERS

	emis = 0.67 	# atmospheric emissivity
	r_a = 1000
	r_s = 1000
	r_v = 250
	v_c = 1e-06	

############ Geometry

	theta = 0.45
	h_s = 0.1			# surface layer depth [m]
	h_D = 0.9			# deep layer
	Jackson = 0.93			# From Jackson et al. 
	r = 1 - Jackson**(h_s*100) 	# Roots by Jackson

############### Physical Constants
	
	sig = 5.67e-08			# SB
	rho_a = 1.25			# density of air [kg/m^3]
	rho_l = 1000			# denisty of water [kg/m^3]
	rho_s = 1000			# density of dry soil [kg/m^3]
	c_ps  = 1000			# heat capacity of dry soil [J/kg/K]
	c_pa = 1003			# heat capacity of dry air [J/kg/K]
	c_pl  =  4182			# heat capacity of water [J/kg/K]
	L = 2257000			# Latent enthalpy of vaporization [J/kg]

################### Combined parameters

	mu_s = rho_l*h_s
	mu_d = rho_l*h_D		# storage capacity of root layer

################# State Variables

	Ts = np.zeros(N)
	m_s = np.zeros(N)
	m_d = np.zeros(N)
	LHF = np.zeros(N)
	TRANSP = np.zeros(N)
	P = np.zeros(N)

	Ts[0] = T_R[0]			# Initial condition for surface temperature
	m_s[0] = 0			# ""			surface moisture
	m_d[0] = 0			# ""			rood mosisture

############## PRECIP STUFF ###################################################

	P_avg = 8				# average precipitation intensity [mm]
	a_1 = .2				# precip frequency [1/days]
	omega = a_1/steps_per_day

####################################################################################

	sec_per_day = 86400		# seconds per day
	dt = 86400./steps_per_day		# time increment (10 chunks per day)
	i = 0

	while i < N-1:

		OLR = sig*(Ts[i]**4)
		DLR = sig*emis*(T_R[i]**4)
		H = c_pa*rho_a*(Ts[i] - T_R[i])/r_a

		es_sat = e_s(Ts[i])
		qs_sat = es_sat*0.622/(P_s - 0.37*es_sat)
		q_diff = qs_sat - q_R[i]			# Humidity difference near surface
		
		if q_diff > 0:
			E_s 	= rho_a*m_s[i]*q_diff/(theta*r_s)		# Surface Transpiration
		else:
			E_s = 0

		if q_diff > 0:

			beta_s = m_s[i]/theta
			beta_r = m_d[i]/theta

			T_s 	= rho_a*r*beta_s*q_diff/(r_v)			# Surface Evaporation [kg/m^2/s]
			T_r     = rho_a*(1-r)*beta_r*q_diff/(r_v)		# Root Transpiration  [kg/m^2/s]
		else:
			T_r = 0
			T_s = 0

		does_rain = rand.rand()
		if does_rain < omega:
			P[i] = rand.gamma(P_avg,scale=1)/1000 # m rainfall (intensity)
		else:
			P[i] = 0
	
		TRANSP[i] = T_r + T_s				# Total transpiration
		LHF[i] = L*(E_s + T_r + T_s)			# Latent heat flux
		Cap = rho_l*v_c*(m_d[i] - m_s[i])		# Capilary exchange
		C_eff = h_s*(c_ps*rho_s + c_pl*m_s[i]*rho_l)	# Effective heat capacity of storage [J/m^2/K]

		dT_dt = (F[i] - OLR + DLR - LHF[i] - H)/C_eff
		dms_dt = P[i]/(h_s*theta) + (Cap - E_s - T_s)/mu_s # precip/(m)
		dmd_dt = -(T_r + Cap)/mu_d

		Ts[i+1] 	= Ts[i] + dT_dt*dt
		m_s[i+1]	= m_s[i] + dms_dt*dt 
		m_d[i+1]	= m_d[i] + dmd_dt*dt


		if m_s[i+1] > m_d[i+1]:
			surface_ro = (m_s[i+1] - m_d[i+1])*h_s/h_D
			m_s[i+1] = m_d[i+1]
			m_d[i+1] = m_d[i+1] + surface_ro

		if m_d[i+1] > theta:
			m_d[i+1] = theta	
		if m_s[i+1] < 0:
			m_s[i+1] = 0
		if m_d[i+1] < 0:
			m_d[i+1] = 0
	
		i+=1

	return(Ts,m_s*theta,m_d*theta,P)
	

F,T,q = make_forcing(10)
Ts,ms,md,P = the_model(F,T,q)

plt.plot(np.arange(len(Ts)),Ts)
plt.show()
plt.savefig('lucas_code.png')
