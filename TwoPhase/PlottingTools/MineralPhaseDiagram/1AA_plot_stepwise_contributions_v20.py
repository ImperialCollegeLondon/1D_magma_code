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
r_cb_i_abs = [0]*t_length

r_mf_i_dum = [0]*t_length
r_cb_i_dum = [0]*t_length
r_mf_t_i_dum = [0]*t_length
r_cb_t_i_dum = [0]*t_length


cb_i = [0]*t_length

mf_check = [0]*t_length
mf_real = [0]*t_length
mf_res = [0]*t_length

cb_check = [0]*t_length
cb_real = [0]*t_length
cb_res = [0]*t_length

phi_compact =[0]*(t_length-1)
phi_real = [0]*(t_length-1)

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
	mf_res[k] =  abs(mf_real[k] - mf_check[k])
	
	cb_check[k] = 18.9 + c_cb_i[k] + r_cb_i[k]
	
	cb_real[k] = cb_i[k]
	
	cb_res[k] = abs(cb_real[k] - cb_check[k]) 
	cb_res[0]=0
	mf_res[0]=0
	r_cb_i_abs[k] = abs(r_cb_i[k])




new_h_mf_t = [0]*(t_length-1)
new_c_mf_t = [0]*(t_length-1)
new_r_mf_t = [0]*(t_length-1)
new_c_cb_t = [0]*(t_length-1)
new_r_cb_t = [0]*(t_length-1)


new_r_mf_t_dum = [0]*(t_length-1)
new_r_cb_t_dum = [0]*(t_length-1)

new_time = [0]*(t_length-1)

phi_compact[0] = 1	
	
for k in range(len(time)-1):

	
	new_time[k] = time[k+1]
	
	new_h_mf_t[k] = h_mf_i[k+1] - h_mf_i[k]
	
	new_c_mf_t[k] = c_mf_i[k+1] - c_mf_i[k]
	
	new_r_mf_t[k] = r_mf_i[k+1] - r_mf_i[k]
	
	new_c_cb_t[k] = c_cb_i[k+1] - c_cb_i[k]
	
	new_r_cb_t[k] = r_cb_i[k+1] - r_cb_i[k]
	
	if k>0:
	
		phi_compact[k] = phi_compact[k-1] + new_c_mf_t[k] #- new_r_mf_t[k]
	
	phi_real[k] = phi_i[k+1]
	
	
	#if k>0:
	#
	#	new_r_mf_t_dum[k] = r_mf_i_dum[k+1] - r_mf_i_dum[k]
	#	new_r_cb_t_dum[k] = r_cb_i_dum[k+1] - r_cb_i_dum[k]
		
	#else:
	#	new_r_mf_t_dum[k] = r_mf_i_dum[k+1] - 0
	#	new_r_cb_t_dum[k] = r_cb_i_dum[k+1] - 0
	
	
	with open("new_step_calc.txt", "w") as file:
		file.write("time\tnew_h_mf_t\tnew_c_mf_t\tnew_r_mf_t\tnew_c_cb_t\tnew_r_cb_t\n")
		for t, new_h_mf, new_c_mf, new_r_mf, new_c_cb, new_r_cb in zip(new_time, new_h_mf_t, new_c_mf_t, new_r_mf_t, new_c_cb_t, new_r_cb_t):
			file.write(f"{t:<10.10f}\t{new_h_mf:<10.10f}\t{new_c_mf:<10.10f}\t{new_r_mf:<10.10f}\t{new_c_cb:<10.10f}\t{new_r_cb}\n")
		
		
	
	
print(phi_real)	
	  

    ### FIGURE ###

#t_length=len(time)-1
#for k in range(t_length):
	#plot(k)
        


################################################################################################
fig.clf()
fig, (ax2,ax3) = plt.subplots(2,1)
fig.set_size_inches(8,12)


# title of figure
ax2.clear()
ax3.clear()

fig.tight_layout()


ax2.plot(new_time,new_h_mf_t, color='b',label='Heating/cooling') 
ax2.plot(new_time,new_r_mf_t, color='k',label='Reactive flow') 
ax2.plot(new_time,new_c_mf_t, color='r',label='Compaction') 
ax2.set_xlim([0, time[-1]])
#ax2.set_ylim([-1, 1])
ax2.set_frame_on(True)
ax2.tick_params(direction='in', top=True, right=True, which='major')
ax2.tick_params(direction='in', top=True, right=True, which='minor')
ax2.minorticks_on()
ax2.legend(loc="best")
ax2.set_xlabel("Time (ka)")
ax2.set_ylabel("Stepwise contribution to melt fraction change (-)") 

ax4 = ax2.twinx()
ax4.plot(new_time,phi_compact, color="#808080", linestyle=':')
ax4.plot(new_time, phi_real, color="#808080")#, linestyle='-')
ax4.set_ylim([0, 1.01])
ax4.set_ylabel("Phi (-)")
ax4.set_frame_on(True)
ax4.tick_params(direction='in', top=True, right=True, which='major')
ax4.tick_params(direction='in', top=True, right=True, which='minor')
ax4.minorticks_on()

ax3.plot(new_time,new_r_cb_t, color='k',label='Reactive flow') 
ax3.plot(new_time,new_c_cb_t, color='r',label='Compaction') 
ax3.set_xlim([0, time[-1]])
#ax3.set_ylim([-1, 1])
ax3.set_frame_on(True)
ax3.tick_params(direction='in', top=True, right=True, which='major')
ax3.tick_params(direction='in', top=True, right=True, which='minor')
ax3.minorticks_on()
ax3.legend(loc="best")
ax3.set_xlabel("Time (ka)")
ax3.set_ylabel("Stepwise contribution to bulk MgO change (-)")

 

fig.tight_layout()
numb=depth_interest
fig.savefig("AA_stepwise_contr_at_%skm.svg" %numb) 











