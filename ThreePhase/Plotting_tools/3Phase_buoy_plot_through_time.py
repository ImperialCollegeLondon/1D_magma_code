### import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from celluloid import Camera
import imageio
import ffmpeg
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



FPSN=5


prompt1="Input the number of the initial input: "
kSTART = int(input(prompt1))

prompt2="Input the number of the final input: "
kEND = int(input(prompt2))

#prompt3="Input the time interval between inputs: "
#t_int = float(input(prompt3))

prompt4="What the critical strength of the rock? "
crit_rock_input = float(input(prompt4))


#Initialise vectors
time=[None]*(kEND-kSTART)

time_1=[]
phi_top_1=[]
phi_bot_1=[]
buoy_top_1=[]
buoy_bot_1=[]
buoy_hphi_top_1=[]
buoy_hphi_bot_1=[]
h_buoy_1=[]
h_rti_1=[]
t2_1=[]
over_hrti_1=[]
over_layer_1=[]
over_total_1=[]
crit_therm_1=[]
crit_strength_1=[]
crit_rock_1=[]

time_2=[]
phi_top_2=[]
phi_bot_2=[]
buoy_top_2=[]
buoy_bot_2=[]
buoy_hphi_top_2=[]
buoy_hphi_bot_2=[]
h_buoy_2=[]
h_rti_2=[]
t2_2=[]
over_hrti_2=[]
over_layer_2=[]
over_total_2=[]
crit_therm_2=[]
crit_strength_2=[]
crit_rock_2=[]

time_3=[]
phi_top_3=[]
phi_bot_3=[]
buoy_top_3=[]
buoy_bot_3=[]
buoy_hphi_top_3=[]
buoy_hphi_bot_3=[]
h_buoy_3=[]
h_rti_3=[]
t2_3=[]
over_hrti_3=[]
over_layer_3=[]
over_total_3=[]
crit_therm_3=[]
crit_strength_3=[]
crit_rock_3=[]

time_4=[]
phi_top_4=[]
phi_bot_4=[]
buoy_top_4=[]
buoy_bot_4=[]
buoy_hphi_top_4=[]
buoy_hphi_bot_4=[]
h_buoy_4=[]
h_rti_4=[]
t2_4=[]
over_hrti_4=[]
over_layer_4=[]
over_total_4=[]
crit_therm_4=[]
crit_strength_4=[]
crit_rock_4=[]

time_5=[]
phi_top_5=[]
phi_bot_5=[]
buoy_top_5=[]
buoy_bot_5=[]
buoy_hphi_top_5=[]
buoy_hphi_bot_5=[]
h_buoy_5=[]
h_rti_5=[]
t2_5=[]
over_hrti_5=[]
over_layer_5=[]
over_total_5=[]
crit_therm_5=[]
crit_strength_5=[]
crit_rock_5=[]

time_6=[]
phi_top_6=[]
phi_bot_6=[]
buoy_top_6=[]
buoy_bot_6=[]
buoy_hphi_top_6=[]
buoy_hphi_bot_6=[]
h_buoy_6=[]
h_rti_6=[]
t2_6=[]
over_hrti_6=[]
over_layer_6=[]
over_total_6=[]
crit_therm_6=[]
crit_strength_6=[]
crit_rock_6=[]

# initialise figure

fig, (ax1) = plt.subplots(1,1)
fig.set_size_inches(30,10)

# define plot
for k in range(kSTART,kEND,1):
	
	print(k)

			
	dtype1 = np.dtype ([('time','f8'),('no','f8'),('prev_no','f8'),('phi_top','f8'),('phi_base','f8'),
				('buoy_top','f8'),('buoy_base','f8'),('buoyCMFtop','f8'),('buoyCMFbase','f8'),
				('buoy_phi','f8'),('buoy_hphi','f8'),('phi_top_km','f8'),('phi_base_km','f8')
				,('buoy_top_km','f8'),('buoy_base_km','f8'),('buoy_CMF_top_km','f8'),('buoy_CMF_base_km','f8')
				,('buoy_phi_m','f8'),('buoy_hphi_m','f8'),('crust_dens','f8'),('layer_dens','f8')
				,('hrti_flag','f8'),('t2','f8'),('hrti_0','f8'),('hrti','f8')
				,('over_rti','f8'),('over_layer','f8'),('over_total','f8'),('crit_therm','f8')
				,('crit_over','f8')])
				

			


	desired_width=4
	if k==0:
		filename="BUOYphi_output_0.txt" 
	else:
		numb=k#"{:0>{}}".format(k)#format(k,desired_width)
		filename="BUOYphi_output_%s.txt" %numb

	data = np.loadtxt(filename, dtype=dtype1, skiprows=1)




	phi_top = np.array([data['phi_top_km']])
	phi_bot = np.array([data['phi_base_km']])
	buoy_top = np.array([data['buoy_top_km']])
	buoy_bot = np.array([data['buoy_base_km']])
	buoy_hphi_top = np.array([data['buoy_CMF_top_km']])
	buoy_hphi_bot = np.array([data['buoy_CMF_base_km']])

	h_buoy = np.array([data['buoy_phi_m']])
	h_rti = np.array([data['hrti']])

	t2 = np.array([data['t2']])
	t2 = t2/(60*60*24*365.25*1000)

	over_hrti = np.array([data['over_rti']])
	over_layer = np.array([data['over_layer']])
	over_total = np.array([data['over_total']])
	crit_therm = np.array([data['crit_therm']])
	crit_strength = np.array([data['crit_over']])
	crit_rock = crit_rock_input*(over_hrti/over_hrti)

	numb_layer = data['no']

	
		



	# calculating time
	time = np.array([data['time']])
	
	
	print(numb_layer.size)

	if numb_layer.size==0:
		print('ok 0')
	elif numb_layer.size==1:
		print('ok 1')
		
		time_1.append(time[0])
		phi_top_1.append(phi_top[0])
		phi_bot_1.append(phi_bot[0])
		buoy_top_1.append(buoy_top[0])
		buoy_bot_1.append(buoy_bot[0])
		buoy_hphi_top_1.append(buoy_hphi_top[0])
		buoy_hphi_bot_1.append(buoy_hphi_bot[0])
		h_buoy_1.append(h_buoy[0])
		h_rti_1.append(h_rti[0])
		t2_1.append(t2[0])
		over_hrti_1.append(over_hrti[0])
		over_layer_1.append(over_layer[0])
		over_total_1.append(over_total[0])
		crit_therm_1.append(crit_therm[0])
		crit_strength_1.append(crit_strength[0])
		crit_rock_1.append(crit_rock[0])	
		
	elif numb_layer.size==2:
		print('ok 2')
		
		
		time_1.append(time[0,0])
		phi_top_1.append(phi_top[0,0])
		phi_bot_1.append(phi_bot[0,0])
		buoy_top_1.append(buoy_top[0,0])
		buoy_bot_1.append(buoy_bot[0,0])
		buoy_hphi_top_1.append(buoy_hphi_top[0,0])
		buoy_hphi_bot_1.append(buoy_hphi_bot[0,0])
		h_buoy_1.append(h_buoy[0,0])
		h_rti_1.append(h_rti[0,0])
		t2_1.append(t2[0,0])
		over_hrti_1.append(over_hrti[0,0])
		over_layer_1.append(over_layer[0,0])
		over_total_1.append(over_total[0,0])
		crit_therm_1.append(crit_therm[0,0])
		crit_strength_1.append(crit_strength[0,0])
		crit_rock_1.append(crit_rock[0,0])
		
		time_2.append(time[0,1])
		phi_top_2.append(phi_top[0,1])
		phi_bot_2.append(phi_bot[0,1])
		buoy_top_2.append(buoy_top[0,1])
		buoy_bot_2.append(buoy_bot[0,1])
		buoy_hphi_top_2.append(buoy_hphi_top[0,1])
		buoy_hphi_bot_2.append(buoy_hphi_bot[0,1])
		h_buoy_2.append(h_buoy[0,1])
		h_rti_2.append(h_rti[0,1])
		t2_2.append(t2[0,1])
		over_hrti_2.append(over_hrti[0,1])
		over_layer_2.append(over_layer[0,1])
		over_total_2.append(over_total[0,1])
		crit_therm_2.append(crit_therm[0,1])
		crit_strength_2.append(crit_strength[0,1])
		crit_rock_2.append(crit_rock[0,1])
		
		
	elif numb_layer.size==3:
		print('ok 3')
		
		time_1.append(time[0,0])
		phi_top_1.append(phi_top[0,0])
		phi_bot_1.append(phi_bot[0,0])
		buoy_top_1.append(buoy_top[0,0])
		buoy_bot_1.append(buoy_bot[0,0])
		buoy_hphi_top_1.append(buoy_hphi_top[0,0])
		buoy_hphi_bot_1.append(buoy_hphi_bot[0,0])
		h_buoy_1.append(h_buoy[0,0])
		h_rti_1.append(h_rti[0,0])
		t2_1.append(t2[0,0])
		over_hrti_1.append(over_hrti[0,0])
		over_layer_1.append(over_layer[0,0])
		over_total_1.append(over_total[0,0])
		crit_therm_1.append(crit_therm[0,0])
		crit_strength_1.append(crit_strength[0,0])
		crit_rock_1.append(crit_rock[0,0])
		
		time_2.append(time[0,1])
		phi_top_2.append(phi_top[0,1])
		phi_bot_2.append(phi_bot[0,1])
		buoy_top_2.append(buoy_top[0,1])
		buoy_bot_2.append(buoy_bot[0,1])
		buoy_hphi_top_2.append(buoy_hphi_top[0,1])
		buoy_hphi_bot_2.append(buoy_hphi_bot[0,1])
		h_buoy_2.append(h_buoy[0,1])
		h_rti_2.append(h_rti[0,1])
		t2_2.append(t2[0,1])
		over_hrti_2.append(over_hrti[0,1])
		over_layer_2.append(over_layer[0,1])
		over_total_2.append(over_total[0,1])
		crit_therm_2.append(crit_therm[0,1])
		crit_strength_2.append(crit_strength[0,1])
		crit_rock_2.append(crit_rock[0,1])
		
		time_3.append(time[0,2])
		phi_top_3.append(phi_top[0,2])
		phi_bot_3.append(phi_bot[0,2])
		buoy_top_3.append(buoy_top[0,2])
		buoy_bot_3.append(buoy_bot[0,2])
		buoy_hphi_top_3.append(buoy_hphi_top[0,2])
		buoy_hphi_bot_3.append(buoy_hphi_bot[0,2])
		h_buoy_3.append(h_buoy[0,2])
		h_rti_3.append(h_rti[0,2])
		t2_3.append(t2[0,2])
		over_hrti_3.append(over_hrti[0,2])
		over_layer_3.append(over_layer[0,2])
		over_total_3.append(over_total[0,2])
		crit_therm_3.append(crit_therm[0,2])
		crit_strength_3.append(crit_strength[0,2])
		crit_rock_3.append(crit_rock[0,2])
		
		
	elif numb_layer.size==4:
		print('ok 4')
		
		time_1.append(time[0,0])
		phi_top_1.append(phi_top[0,0])
		phi_bot_1.append(phi_bot[0,0])
		buoy_top_1.append(buoy_top[0,0])
		buoy_bot_1.append(buoy_bot[0,0])
		buoy_hphi_top_1.append(buoy_hphi_top[0,0])
		buoy_hphi_bot_1.append(buoy_hphi_bot[0,0])
		h_buoy_1.append(h_buoy[0,0])
		h_rti_1.append(h_rti[0,0])
		t2_1.append(t2[0,0])
		over_hrti_1.append(over_hrti[0,0])
		over_layer_1.append(over_layer[0,0])
		over_total_1.append(over_total[0,0])
		crit_therm_1.append(crit_therm[0,0])
		crit_strength_1.append(crit_strength[0,0])
		crit_rock_1.append(crit_rock[0,0])
		
		time_2.append(time[0,1])
		phi_top_2.append(phi_top[0,1])
		phi_bot_2.append(phi_bot[0,1])
		buoy_top_2.append(buoy_top[0,1])
		buoy_bot_2.append(buoy_bot[0,1])
		buoy_hphi_top_2.append(buoy_hphi_top[0,1])
		buoy_hphi_bot_2.append(buoy_hphi_bot[0,1])
		h_buoy_2.append(h_buoy[0,1])
		h_rti_2.append(h_rti[0,1])
		t2_2.append(t2[0,1])
		over_hrti_2.append(over_hrti[0,1])
		over_layer_2.append(over_layer[0,1])
		over_total_2.append(over_total[0,1])
		crit_therm_2.append(crit_therm[0,1])
		crit_strength_2.append(crit_strength[0,1])
		crit_rock_2.append(crit_rock[0,1])
		
		time_3.append(time[0,2])
		phi_top_3.append(phi_top[0,2])
		phi_bot_3.append(phi_bot[0,2])
		buoy_top_3.append(buoy_top[0,2])
		buoy_bot_3.append(buoy_bot[0,2])
		buoy_hphi_top_3.append(buoy_hphi_top[0,2])
		buoy_hphi_bot_3.append(buoy_hphi_bot[0,2])
		h_buoy_3.append(h_buoy[0,2])
		h_rti_3.append(h_rti[0,2])
		t2_3.append(t2[0,2])
		over_hrti_3.append(over_hrti[0,2])
		over_layer_3.append(over_layer[0,2])
		over_total_3.append(over_total[0,2])
		crit_therm_3.append(crit_therm[0,2])
		crit_strength_3.append(crit_strength[0,2])
		crit_rock_3.append(crit_rock[0,2])
		
		time_4.append(time[0,3])
		phi_top_4.append(phi_top[0,3])
		phi_bot_4.append(phi_bot[0,3])
		buoy_top_4.append(buoy_top[0,3])
		buoy_bot_4.append(buoy_bot[0,3])
		buoy_hphi_top_4.append(buoy_hphi_top[0,3])
		buoy_hphi_bot_4.append(buoy_hphi_bot[0,3])
		h_buoy_4.append(h_buoy[0,3])
		h_rti_4.append(h_rti[0,3])
		t2_4.append(t2[0,3])
		over_hrti_4.append(over_hrti[0,3])
		over_layer_4.append(over_layer[0,3])
		over_total_4.append(over_total[0,3])
		crit_therm_4.append(crit_therm[0,3])
		crit_strength_4.append(crit_strength[0,3])
		crit_rock_4.append(crit_rock[0,3])
		
	elif numb_layer.size==5:
		print('ok 5')
		
		time_1.append(time[0,0])
		phi_top_1.append(phi_top[0,0])
		phi_bot_1.append(phi_bot[0,0])
		buoy_top_1.append(buoy_top[0,0])
		buoy_bot_1.append(buoy_bot[0,0])
		buoy_hphi_top_1.append(buoy_hphi_top[0,0])
		buoy_hphi_bot_1.append(buoy_hphi_bot[0,0])
		h_buoy_1.append(h_buoy[0,0])
		h_rti_1.append(h_rti[0,0])
		t2_1.append(t2[0,0])
		over_hrti_1.append(over_hrti[0,0])
		over_layer_1.append(over_layer[0,0])
		over_total_1.append(over_total[0,0])
		crit_therm_1.append(crit_therm[0,0])
		crit_strength_1.append(crit_strength[0,0])
		crit_rock_1.append(crit_rock[0,0])
		
		time_2.append(time[0,1])
		phi_top_2.append(phi_top[0,1])
		phi_bot_2.append(phi_bot[0,1])
		buoy_top_2.append(buoy_top[0,1])
		buoy_bot_2.append(buoy_bot[0,1])
		buoy_hphi_top_2.append(buoy_hphi_top[0,1])
		buoy_hphi_bot_2.append(buoy_hphi_bot[0,1])
		h_buoy_2.append(h_buoy[0,1])
		h_rti_2.append(h_rti[0,1])
		t2_2.append(t2[0,1])
		over_hrti_2.append(over_hrti[0,1])
		over_layer_2.append(over_layer[0,1])
		over_total_2.append(over_total[0,1])
		crit_therm_2.append(crit_therm[0,1])
		crit_strength_2.append(crit_strength[0,1])
		crit_rock_2.append(crit_rock[0,1])
		
		time_3.append(time[0,2])
		phi_top_3.append(phi_top[0,2])
		phi_bot_3.append(phi_bot[0,2])
		buoy_top_3.append(buoy_top[0,2])
		buoy_bot_3.append(buoy_bot[0,2])
		buoy_hphi_top_3.append(buoy_hphi_top[0,2])
		buoy_hphi_bot_3.append(buoy_hphi_bot[0,2])
		h_buoy_3.append(h_buoy[0,2])
		h_rti_3.append(h_rti[0,2])
		t2_3.append(t2[0,2])
		over_hrti_3.append(over_hrti[0,2])
		over_layer_3.append(over_layer[0,2])
		over_total_3.append(over_total[0,2])
		crit_therm_3.append(crit_therm[0,2])
		crit_strength_3.append(crit_strength[0,2])
		crit_rock_3.append(crit_rock[0,2])
		
		time_4.append(time[0,3])
		phi_top_4.append(phi_top[0,3])
		phi_bot_4.append(phi_bot[0,3])
		buoy_top_4.append(buoy_top[0,3])
		buoy_bot_4.append(buoy_bot[0,3])
		buoy_hphi_top_4.append(buoy_hphi_top[0,3])
		buoy_hphi_bot_4.append(buoy_hphi_bot[0,3])
		h_buoy_4.append(h_buoy[0,3])
		h_rti_4.append(h_rti[0,3])
		t2_4.append(t2[0,3])
		over_hrti_4.append(over_hrti[0,3])
		over_layer_4.append(over_layer[0,3])
		over_total_4.append(over_total[0,3])
		crit_therm_4.append(crit_therm[0,3])
		crit_strength_4.append(crit_strength[0,3])
		crit_rock_4.append(crit_rock[0,3])
		
		time_5.append(time[0,4])
		phi_top_5.append(phi_top[0,4])
		phi_bot_5.append(phi_bot[0,4])
		buoy_top_5.append(buoy_top[0,4])
		buoy_bot_5.append(buoy_bot[0,4])
		buoy_hphi_top_5.append(buoy_hphi_top[0,4])
		buoy_hphi_bot_5.append(buoy_hphi_bot[0,4])
		h_buoy_5.append(h_buoy[0,4])
		h_rti_5.append(h_rti[0,4])
		t2_5.append(t2[0,4])
		over_hrti_5.append(over_hrti[0,4])
		over_layer_5.append(over_layer[0,4])
		over_total_5.append(over_total[0,4])
		crit_therm_5.append(crit_therm[0,4])
		crit_strength_5.append(crit_strength[0,4])
		crit_rock_5.append(crit_rock[0,4])
		
	elif numb_layer.size==6:
		print('ok 6')
		
		time_1.append(time[0,0])
		phi_top_1.append(phi_top[0,0])
		phi_bot_1.append(phi_bot[0,0])
		buoy_top_1.append(buoy_top[0,0])
		buoy_bot_1.append(buoy_bot[0,0])
		buoy_hphi_top_1.append(buoy_hphi_top[0,0])
		buoy_hphi_bot_1.append(buoy_hphi_bot[0,0])
		h_buoy_1.append(h_buoy[0,0])
		h_rti_1.append(h_rti[0,0])
		t2_1.append(t2[0,0])
		over_hrti_1.append(over_hrti[0,0])
		over_layer_1.append(over_layer[0,0])
		over_total_1.append(over_total[0,0])
		crit_therm_1.append(crit_therm[0,0])
		crit_strength_1.append(crit_strength[0,0])
		crit_rock_1.append(crit_rock[0,0])
		
		time_2.append(time[0,1])
		phi_top_2.append(phi_top[0,1])
		phi_bot_2.append(phi_bot[0,1])
		buoy_top_2.append(buoy_top[0,1])
		buoy_bot_2.append(buoy_bot[0,1])
		buoy_hphi_top_2.append(buoy_hphi_top[0,1])
		buoy_hphi_bot_2.append(buoy_hphi_bot[0,1])
		h_buoy_2.append(h_buoy[0,1])
		h_rti_2.append(h_rti[0,1])
		t2_2.append(t2[0,1])
		over_hrti_2.append(over_hrti[0,1])
		over_layer_2.append(over_layer[0,1])
		over_total_2.append(over_total[0,1])
		crit_therm_2.append(crit_therm[0,1])
		crit_strength_2.append(crit_strength[0,1])
		crit_rock_2.append(crit_rock[0,1])
		
		time_3.append(time[0,2])
		phi_top_3.append(phi_top[0,2])
		phi_bot_3.append(phi_bot[0,2])
		buoy_top_3.append(buoy_top[0,2])
		buoy_bot_3.append(buoy_bot[0,2])
		buoy_hphi_top_3.append(buoy_hphi_top[0,2])
		buoy_hphi_bot_3.append(buoy_hphi_bot[0,2])
		h_buoy_3.append(h_buoy[0,2])
		h_rti_3.append(h_rti[0,2])
		t2_3.append(t2[0,2])
		over_hrti_3.append(over_hrti[0,2])
		over_layer_3.append(over_layer[0,2])
		over_total_3.append(over_total[0,2])
		crit_therm_3.append(crit_therm[0,2])
		crit_strength_3.append(crit_strength[0,2])
		crit_rock_3.append(crit_rock[0,2])
		
		time_4.append(time[0,3])
		phi_top_4.append(phi_top[0,3])
		phi_bot_4.append(phi_bot[0,3])
		buoy_top_4.append(buoy_top[0,3])
		buoy_bot_4.append(buoy_bot[0,3])
		buoy_hphi_top_4.append(buoy_hphi_top[0,3])
		buoy_hphi_bot_4.append(buoy_hphi_bot[0,3])
		h_buoy_4.append(h_buoy[0,3])
		h_rti_4.append(h_rti[0,3])
		t2_4.append(t2[0,3])
		over_hrti_4.append(over_hrti[0,3])
		over_layer_4.append(over_layer[0,3])
		over_total_4.append(over_total[0,3])
		crit_therm_4.append(crit_therm[0,3])
		crit_strength_4.append(crit_strength[0,3])
		crit_rock_4.append(crit_rock[0,3])
		
		time_5.append(time[0,4])
		phi_top_5.append(phi_top[0,4])
		phi_bot_5.append(phi_bot[0,4])
		buoy_top_5.append(buoy_top[0,4])
		buoy_bot_5.append(buoy_bot[0,4])
		buoy_hphi_top_5.append(buoy_hphi_top[0,4])
		buoy_hphi_bot_5.append(buoy_hphi_bot[0,4])
		h_buoy_5.append(h_buoy[0,4])
		h_rti_5.append(h_rti[0,4])
		t2_5.append(t2[0,4])
		over_hrti_5.append(over_hrti[0,4])
		over_layer_5.append(over_layer[0,4])
		over_total_5.append(over_total[0,4])
		crit_therm_5.append(crit_therm[0,4])
		crit_strength_5.append(crit_strength[0,4])
		crit_rock_5.append(crit_rock[0,4])
		
		time_6.append(time[0,5])
		phi_top_6.append(phi_top[0,5])
		phi_bot_6.append(phi_bot[0,5])
		buoy_top_6.append(buoy_top[0,5])
		buoy_bot_6.append(buoy_bot[0,5])
		buoy_hphi_top_6.append(buoy_hphi_top[0,5])
		buoy_hphi_bot_6.append(buoy_hphi_bot[0,5])
		h_buoy_6.append(h_buoy[0,5])
		h_rti_6.append(h_rti[0,5])
		t2_6.append(t2[0,5])
		over_hrti_6.append(over_hrti[0,5])
		over_layer_6.append(over_layer[0,5])
		over_total_6.append(over_total[0,5])
		crit_therm_6.append(crit_therm[0,5])
		crit_strength_6.append(crit_strength[0,5])
		crit_rock_6.append(crit_rock[0,5])
		
	else:
		raise Exception("More melt layers then is accounted for in code! Edit code")
	
	
####################################################here#

#########################################OVERVIEW MELT LAYERS########################

ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Overview of melt layers')

# axis 1
ax1.plot(time_1,phi_top_1, 'ro', label="Top of melt layer 1") 
ax1.plot(time_1,phi_bot_1, 'bo', label="Base of melt layer 1") 

ax1.plot(time_2,phi_top_2, 'rx', label="Top of melt layer 2") 
ax1.plot(time_2,phi_bot_2, 'bx', label="Base of melt layer 2") 

ax1.plot(time_3,phi_top_3, 'r^', label="Top of melt layer 3") 
ax1.plot(time_3,phi_bot_3, 'b^', label="Base of melt layer 3")  

ax1.plot(time_4,phi_top_4, 'r*', label="Top of melt layer 4") 
ax1.plot(time_4,phi_bot_4, 'b*', label="Base of melt layer 4")

ax1.plot(time_5,phi_top_5, 'r+', label="Top of melt layer 5") 
ax1.plot(time_5,phi_bot_5, 'b+', label="Base of melt layer 5")

ax1.plot(time_6,phi_top_6, 'rs', label="Top of melt layer 6") 
ax1.plot(time_6,phi_bot_6, 'bs', label="Base of melt layer 6")


ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
ax1.set_ylim([-35, -10]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_phidepths_overview.svg") 

#########################################LAYER 1 MELT LAYERS########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 1')

# axis 1
ax1.plot(time_1,phi_top_1, 'ro', label="Top of melt") 
ax1.plot(time_1,phi_bot_1, 'bo', label="Base of melt") 

ax1.plot(time_1,buoy_top_1, 'rx', label="Top of buoy melt") 
ax1.plot(time_1,buoy_bot_1, 'bx', label="Base of buoy melt") 

ax1.plot(time_1,buoy_hphi_top_1, 'r^', label="Top of buoy magma") 
ax1.plot(time_1,buoy_hphi_bot_1, 'b^', label="Base of buoy magma") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
ax1.set_xlim([0, 1000])
ax1.set_ylim([-35, -20]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_phidepths_Layer1_zoom.svg") 

#########################################LAYER 2 MELT LAYERS########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 2')

# axis 1
ax1.plot(time_2,phi_top_2, 'ro', label="Top of melt") 
ax1.plot(time_2,phi_bot_2, 'bo', label="Base of melt") 

ax1.plot(time_2,buoy_top_2, 'rx', label="Top of buoy melt") 
ax1.plot(time_2,buoy_bot_2, 'bx', label="Base of buoy melt") 

ax1.plot(time_2,buoy_hphi_top_2, 'r^', label="Top of buoy magma") 
ax1.plot(time_2,buoy_hphi_bot_2, 'b^', label="Base of buoy magma") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
ax1.set_ylim([-35, -10]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_phidepths_Layer2.svg") 

#########################################LAYER 3 MELT LAYERS########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 3')

# axis 1
ax1.plot(time_3,phi_top_3, 'ro', label="Top of melt") 
ax1.plot(time_3,phi_bot_3, 'bo', label="Base of melt") 

ax1.plot(time_3,buoy_top_3, 'rx', label="Top of buoy melt") 
ax1.plot(time_3,buoy_bot_3, 'bx', label="Base of buoy melt") 

ax1.plot(time_3,buoy_hphi_top_3, 'r^', label="Top of buoy magma") 
ax1.plot(time_3,buoy_hphi_bot_3, 'b^', label="Base of buoy magma") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
ax1.set_ylim([-35, -10]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_phidepths_Layer3.svg") 

#########################################LAYER 4 MELT LAYERS########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 4')

# axis 1
ax1.plot(time_4,phi_top_4, 'ro', label="Top of melt") 
ax1.plot(time_4,phi_bot_4, 'bo', label="Base of melt") 

ax1.plot(time_4,buoy_top_4, 'rx', label="Top of buoy melt") 
ax1.plot(time_4,buoy_bot_4, 'bx', label="Base of buoy melt") 

ax1.plot(time_4,buoy_hphi_top_4, 'r^', label="Top of buoy magma") 
ax1.plot(time_4,buoy_hphi_bot_4, 'b^', label="Base of buoy magma") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
ax1.set_ylim([-35, -10]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_phidepths_Layer4.svg") 


#########################################LAYER 5 MELT LAYERS########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 5')

# axis 1
ax1.plot(time_5,phi_top_5, 'ro', label="Top of melt") 
ax1.plot(time_5,phi_bot_5, 'bo', label="Base of melt") 

ax1.plot(time_5,buoy_top_5, 'rx', label="Top of buoy melt") 
ax1.plot(time_5,buoy_bot_5, 'bx', label="Base of buoy melt") 

ax1.plot(time_5,buoy_hphi_top_5, 'r^', label="Top of buoy magma") 
ax1.plot(time_5,buoy_hphi_bot_5, 'b^', label="Base of buoy magma") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
ax1.set_ylim([-35, -10]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_phidepths_Layer5.svg") 


#########################################LAYER 6 MELT LAYERS########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 6')

# axis 1
ax1.plot(time_6,phi_top_6, 'ro', label="Top of melt") 
ax1.plot(time_6,phi_bot_6, 'bo', label="Base of melt") 

ax1.plot(time_6,buoy_top_6, 'rx', label="Top of buoy melt") 
ax1.plot(time_6,buoy_bot_6, 'bx', label="Base of buoy melt") 

ax1.plot(time_6,buoy_hphi_top_6, 'r^', label="Top of buoy magma") 
ax1.plot(time_6,buoy_hphi_bot_6, 'b^', label="Base of buoy magma") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
ax1.set_ylim([-35, -10]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_phidepths_Layer6.svg") 

#########################################LAYER 1 hb hrti########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 1')

# axis 1
ax1.plot(time_1,h_buoy_1, 'ro', label="hb") 
ax1.plot(time_1,h_rti_1, 'b^', label="hrti") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Thickness (m)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_hbhrti_Layer1.svg") 

#########################################LAYER 2 hb hrti########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 2')

# axis 1
ax1.plot(time_2,h_buoy_2, 'ro', label="hb") 
ax1.plot(time_2,h_rti_2, 'b^', label="hrti") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Thickness (m)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_hbhrti_Layer2.svg") 

#########################################LAYER 3 hb hrti########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 3')

# axis 1
ax1.plot(time_3,h_buoy_3, 'ro', label="hb") 
ax1.plot(time_3,h_rti_3, 'b^', label="hrti") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Thickness (m)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_hbhrti_Layer3.svg") 

#########################################LAYER 4 hb hrti########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 4')

# axis 1
ax1.plot(time_4,h_buoy_4, 'ro', label="hb") 
ax1.plot(time_4,h_rti_4, 'b^', label="hrti") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Thickness (m)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_hbhrti_Layer4.svg") 

#########################################LAYER 5 hb hrti########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 5')

# axis 1
ax1.plot(time_5,h_buoy_5, 'ro', label="hb") 
ax1.plot(time_5,h_rti_5, 'b^', label="hrti") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Thickness (m)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_hbhrti_Layer5.svg") 

#########################################LAYER 6 hb hrti########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 6')

# axis 1
ax1.plot(time_6,h_buoy_6, 'ro', label="hb") 
ax1.plot(time_6,h_rti_6, 'b^', label="hrti") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Thickness (m)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_hbhrti_Layer6.svg") 

#########################################LAYER 1 t2########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 1')

# axis 1
ax1.plot(time_1,t2_1, 'ko') 


ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("t2 (ka)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_t2_Layer1.svg") 

#########################################LAYER 2 t2########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 2')

# axis 1
ax1.plot(time_2,t2_2, 'ko') 


ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("t2 (ka)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_t2_Layer2.svg") 

#########################################LAYER 3 t2########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 3')

# axis 1
ax1.plot(time_3,t2_3, 'ko') 


ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("t2 (ka)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_t2_Layer3.svg") 

#########################################LAYER 4 t2########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 4')

# axis 1
ax1.plot(time_4,t2_4, 'ko') 


ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("t2 (ka)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_t2_Layer4.svg") 

#########################################LAYER 1 ovepressure########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 1')

# axis 1
ax1.plot(time_1,over_hrti_1, 'bx', label="overpressure hrti") 
ax1.plot(time_1,over_layer_1, 'b+', label="overpressure hb") 
ax1.plot(time_1,over_total_1, 'b^', label="total overpressure") 

ax1.plot(time_1,crit_therm_1, 'ro', label="Thermal death criteria") 
ax1.plot(time_1,crit_rock_1, 'rs', label="Rock strength") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Overpressure (Pa)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_overp_Layer1.svg") 

#########################################LAYER 2 ovepressure########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 2')

# axis 1
ax1.plot(time_2,over_hrti_2, 'bx', label="overpressure hrti") 
ax1.plot(time_2,over_layer_2, 'b+', label="overpressure hb") 
ax1.plot(time_2,over_total_2, 'b^', label="total overpressure") 

ax1.plot(time_2,crit_therm_2, 'ro', label="Thermal death criteria") 
ax1.plot(time_2,crit_rock_2, 'rs', label="Rock strength") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Overpressure (Pa)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_overp_Layer2.svg") 

#########################################LAYER 3 ovepressure########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 3')

# axis 1
ax1.plot(time_3,over_hrti_3, 'bx', label="overpressure hrti") 
ax1.plot(time_3,over_layer_3, 'b+', label="overpressure hb") 
ax1.plot(time_3,over_total_3, 'b^', label="total overpressure") 

ax1.plot(time_3,crit_therm_3, 'ro', label="Thermal death criteria") 
ax1.plot(time_3,crit_rock_3, 'rs', label="Rock strength") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Overpressure (Pa)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_overp_Layer3.svg") 

#########################################LAYER 4 ovepressure########################
ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

plt.title('Layer 4')

# axis 1
ax1.plot(time_4,over_hrti_4, 'bx', label="overpressure hrti") 
ax1.plot(time_4,over_layer_4, 'b+', label="overpressure hb") 
ax1.plot(time_4,over_total_4, 'b^', label="total overpressure") 

ax1.plot(time_4,crit_therm_4, 'ro', label="Thermal death criteria") 
ax1.plot(time_4,crit_rock_4, 'rs', label="Rock strength") 

ax1.legend()
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Overpressure (Pa)",fontsize=16)  

fig.tight_layout()

fig.savefig("A_fig_buoy_output_time_3Phase_overp_Layer4.svg") 




