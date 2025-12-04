### import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from celluloid import Camera
import imageio
import ffmpeg
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

### import input file
dtype2 = np.dtype ([('name','str'),('value','f8')])
input1 = np.loadtxt('1AA_Plot_video_output.txt', dtype=dtype2, delimiter=',')#, skiprows=1)
input_value=input1['value']
# start file for output
kSTART=int(input_value[0])
# end file for output
kEND=int(input("The last output no.: "))#int(input_value[1])
# time interval between outputs
t_int=input_value[2]
# top depth of the plots
top_plot=input_value[3]
# bottom of the plots
bot_plot=input_value[4]
#min mgo
min_mgo=input_value[5]
#max mgo
max_mgo=input_value[6]
#min sio2
min_sio2=input_value[7]
#max sio2
max_sio2=input_value[8]
# frames per second
FPSN = input_value[9]
time=[None]*(kEND-kSTART)

### plot initialising
#fig, (ax1, ax2, ax3) = plt.subplots(1,3) #camera method
# camera=Camera(fig) #camera method


prompt='What depth would you like to take data from? '
depth_interest = float(input(prompt))



#frames1=[]
fig, (ax2, ax3) = plt.subplots(2,1)
fig.set_size_inches(8,12)


t_length=len(time)
phi_i = [0]*t_length
h_mf_i = [0]*t_length
c_mf_i = [0]*t_length
r_mf_i = [0]*t_length
c_cb_i = [0]*t_length
r_cb_i = [0]*t_length

h_mf_t_i = [0]*t_length
c_mf_t_i = [0]*t_length
r_mf_t_i = [0]*t_length
c_cb_t_i = [0]*t_length
r_cb_t_i = [0]*t_length

cb_i = [0]*t_length

mf_check = [0]*t_length
mf_real = [0]*t_length

cb_check = [0]*t_length
cb_real = [0]*t_length

# define plot
#def plot(k):
t_length=len(time)
for k in range(len(time)):
	print(k)
	if k==0:
	# DEPTH, PHI, H, T, TS, TL, CB, CB SIO2, CB MGO, CS MGO, CL MGO, OL, PX
		dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
				('Tl','f8'),('cb','f8'),('cb_mg','f8'),
				('cs_mg','f8'),('cl_mg','f8'),('ol','f8'),('opx','f8'),
				('cpx','f8'),('feld','f8'),('h_mf','f8'),('c_mf','f8')
						,('r_mf','f8'),('h_mf_t','f8'),('c_mf_t','f8')
						,('r_mf_t','f8'),('c_cb','f8'),('r_cb','f8'),('c_cb_t','f8')
						,('r_cb_t','f8')])


	else:

	# DEPTH, PHI, H, T, TS, TL, CB, CB SIO2, CB MGO, CS MGO, CL MGO, OL, PX
		dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
				('Tl','f8'),('cb','f8'),('cb_mg','f8'),
				('cs_mg','f8'),('cl_mg','f8'),('ol','f8'),('opx','f8'),
				('cpx','f8'),('feld','f8'),('dens','f8'),('h_mf','f8'),('c_mf','f8')
						,('r_mf','f8'),('r_mf_dum','f8'),('h_mf_t','f8'),('c_mf_t','f8')
						,('r_mf_t','f8'),('r_mf_t_dum','f8'),('c_cb','f8'),('r_cb','f8'),('r_cb_dum','f8'),('c_cb_t','f8')
						,('r_cb_t','f8'),('r_cb_t_dum','f8')])
        

	desired_width=4
	if k==0:
		filename="output_0_CELLS.txt" 
	else:
		numb=k#"{:0>{}}".format(k)#format(k,desired_width)
		filename="output_%s_CELLS.txt" %numb
		#filepart={'output','.txt'}
		#filename=numb.join(filepart)
		# print(filename)
	
	data = np.loadtxt(filename, dtype=dtype1, skiprows=2)

	# converting solidus
	#solidus=data['ts']-273.15



	# calculating time
	time[k] = t_int*k/1000
	depth = data['z']
	phi = data['phi']
	
	
	#heating/cooling mf - cumu
	h_mf = data['h_mf']
	#compaction mf - cumu
	c_mf = data['c_mf']
	#reactive flow mf - cumu
	if k>0:
		r_mf = data['r_mf_dum']
	else:
		r_mf = data['r_mf']	
	#compaction cb - cumu
	c_cb = data['c_cb']
	#reactive flow cb - cumu
	if k>0:
		r_cb = data['r_cb_dum']
	else:
		r_cb = data['r_cb']
	#heating/cooling mf - step
	h_mf_t = data['h_mf_t']
	#compaction mf - step
	c_mf_t = data['c_mf_t']
	#reactive flow mf - step
	if k>0:
		r_mf_t = data['r_mf_t_dum']
	else:
		r_mf_t = data['r_mf_t']	
	#compaction cb - step
	c_cb_t = data['c_cb_t']
	#reactive flow cb - step
	if k>0:
		r_cb_t = data['r_cb_t_dum']
	else:
		r_cb_t = data['r_cb_t']
		
		
		
	
	cb = data['cb_mg']
	
	#t_length=len(time)-1
	
	

	
	#print(len(time))

	for i in range(len(depth)):
		#print(depth[i], depth_interest)
		if depth[i] == depth_interest:
			
			phi_i[k] = phi[i]
			h_mf_i[k] = h_mf[i]
			c_mf_i[k] = c_mf[i]
			r_mf_i[k] = r_mf[i]
			c_cb_i[k] = c_cb[i]
			r_cb_i[k] = r_cb[i]
			
			h_mf_t_i[k] = h_mf_t[i]
			c_mf_t_i[k] = c_mf_t[i]
			r_mf_t_i[k] = r_mf_t[i]
			c_cb_t_i[k] = c_cb_t[i]
			r_cb_t_i[k] = r_cb_t[i]
			
			cb_i[k] = cb[i]
			break
			
			
		
	mf_check[k] = 1 + h_mf_i[k] + c_mf_i[k] + r_mf_i[k]
	mf_real[k] = phi_i[k]
	
	cb_check[k] = 18.9 + c_cb_i[k] + r_cb_i[k]
	
	cb_real[k] = cb_i[k]





#cumu_total_end_mf = abs(h_mf_i[-1]) + abs(c_mf_i[-1]) + abs(r_mf_i[-1])
#cumu_total_end_cb = abs(c_cb_i[-1]) + abs(r_cb_i[-1])
new_h_mf_t = [0]*(t_length-1)
new_c_mf_t = [0]*(t_length-1)
new_r_mf_t = [0]*(t_length-1)
new_c_cb_t = [0]*(t_length-1)
new_r_cb_t = [0]*(t_length-1)
new_time = [0]*(t_length-1)
cumu_total_end_mf = [0]*(t_length-1)
cumu_total_end_cb = [0]*(t_length-1)	
	
for k in range(len(time)-1):

	
	new_time[k] = time[k+1]
	
	new_h_mf_t[k] = h_mf_i[k+1] - h_mf_i[k]
	
	new_c_mf_t[k] = c_mf_i[k+1] - c_mf_i[k]
	
	new_r_mf_t[k] = r_mf_i[k+1] - r_mf_i[k]
	
	new_c_cb_t[k] = c_cb_i[k+1] - c_cb_i[k]
	
	new_r_cb_t[k] = r_cb_i[k+1] - r_cb_i[k]
	
	cumu_total_end_mf[k] = abs(new_h_mf_t[k]) + abs(new_c_mf_t[k]) + abs(new_r_mf_t[k])
	cumu_total_end_cb[k] = abs(new_c_cb_t[k]) + abs(new_r_cb_t[k])
	


scale_h_mf = [0]*(t_length-1)	
scale_c_mf = [0]*(t_length-1)
scale_r_mf = [0]*(t_length-1)
	
scale_c_cb = [0]*(t_length-1)
scale_r_cb = [0]*(t_length-1)
	

	
for k in range(len(time)-1):
	
	if cumu_total_end_mf[k]>0:
		scale_h_mf[k] = new_h_mf_t[k]/cumu_total_end_mf[k]	
		scale_c_mf[k] = new_c_mf_t[k]/cumu_total_end_mf[k]
		scale_r_mf[k] = new_r_mf_t[k]/cumu_total_end_mf[k]
	else:
		scale_h_mf[k] = 0	
		scale_c_mf[k] = 0
		scale_r_mf[k] = 0
	
	
	if cumu_total_end_cb[k]>0:		
		scale_c_cb[k] = new_c_cb_t[k]/cumu_total_end_cb[k]
		scale_r_cb[k] = new_r_cb_t[k]/cumu_total_end_cb[k]
	else:
		scale_c_cb[k] = 0
		scale_r_cb[k] = 0
	
	with open("inc_total_change.txt", "w") as file:
		file.write("time\th_mf\tc_mf\tr_mf\tc_cb\tr_cb\n")
		for t, h_mf, c_mf, r_mf, c_cb, r_cb in zip(new_time, scale_h_mf, scale_c_mf, scale_r_mf, scale_c_cb, scale_r_cb):
			file.write(f"{t:<10.10f}\t{h_mf:<10.10f}\t{c_mf:<10.10f}\t{r_mf:<10.10f}\t{c_cb:<10.10f}\t{r_cb}\n")
	
	  

    ### FIGURE ###

#t_length=len(time)-1
#for k in range(t_length):
	#plot(k)
        

ax2.clear()
ax3.clear()

# title of figure

fig.tight_layout()
#ax2.set_title("Cumulative contribution" ) 
   

# axis 2

#ax2.plot(new_time,new_h_mf_t, 'b:') 
#ax2.plot(new_time,new_c_mf_t, 'r:') 
#ax2.plot(new_time,new_r_mf_t, 'k:') 
ax2.plot(new_time,scale_h_mf, 'b', label = 'Heating/cooling') 
ax2.plot(new_time,scale_c_mf, 'r', label = 'Compaction')  
ax2.plot(new_time,scale_r_mf, 'k', label = 'Reactive flow')  

ax2.set_xlim([0, time[-1]])
ax2.set_ylim([-1, 1])
#ax2.xaxis.set_minor_locator(MultipleLocator(0.1))    
#ax2.set_ylim([bot_plot, top_plot])
ax2.set_frame_on(True)
ax2.tick_params(direction='in', top=True, right=True, which='major')
ax2.tick_params(direction='in', top=True, right=True, which='minor')
ax2.minorticks_on()
ax2.legend(loc="best")
ax2.set_xlabel("Time (ka)")
ax2.set_ylabel("Scaled stepwise contribution to melt fraction change")   


#ax22 = ax2.twinx()
#ax22.plot(time,phi_i, color='k',marker='+',linestyle=':',linewidth=2)
#ax22.set_ylabel("Melt Fraction (-)") 
#ax22.set_xlim([0, time[-1]])  
#ax22.set_ylim([0, 1.01])

# axis 3

ax3.plot(new_time,scale_c_cb, 'r', label = 'Compaction')  
ax3.plot(new_time,scale_r_cb, 'k', label = 'Reactive flow')
ax3.legend(loc="best")  
ax3.set_xlim([0, time[-1]])
ax3.set_ylim([-1, 1])
ax3.set_frame_on(True)
ax3.tick_params(direction='in', top=True, right=True, which='major')
ax3.tick_params(direction='in', top=True, right=True, which='minor')
ax3.minorticks_on()
#ax3.set_ylim([bot_plot, top_plot])
ax3.set_xlabel("Time (ka)")
ax3.set_ylabel("Scaled stepwise contribution to bulk MgO change (-)")  

#ax33 = ax3.twinx()
#ax33.plot(time,phi_i, color='k',marker='+',linestyle=':',linewidth=2)
#ax33.set_ylabel("Melt Fraction (-)") 
#ax33.set_xlim([0, time[-1]])  
#ax33.set_ylim([0, 1.01])



fig.tight_layout()
numb=depth_interest
fig.savefig("AA_scaled_stepwise_cont_at_%skm.svg" %numb)    


