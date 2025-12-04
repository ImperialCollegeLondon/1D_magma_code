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
kSTART=int(input("The first output no.: "))
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
fig, (ax2,ax3,ax4,ax5) = plt.subplots(4,1)
fig.set_size_inches(8,24)


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
cl_i = []*t_length
time_cl = []*t_length
cs_i = []*t_length
time_cs = []*t_length
cb_initial = [0]*t_length

ol_i = [0]*t_length
opx_i = [0]*t_length
cpx_i = [0]*t_length
feld_i = [0]*t_length

T_i = [0]*t_length
Ts_i = [0]*t_length
Tl_i = [0]*t_length

h_i = [0]*t_length
hs_i = [0]*t_length
hl_i = [0]*t_length



liq_i=0
sol_i=0

temp_cs=[]*t_length
temp_cl=[]*t_length

# define plot
#def plot(k):

t_length=len(time)

j=0


for k in range(kSTART, kEND,1):
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
	time[j] = t_int*k/1000
	depth = data['z']
	phi = data['phi']
		
	
	
	cb = data['cb_mg']
	cl = data['cl_mg']
	cs = data['cs_mg']
	
	ol = data['ol']
	opx = data['opx']
	cpx = data['cpx']
	feld = data['feld']
	
	T = data['T']
	Ts = data['Ts'] 
	Tl = data['Tl']
	
	h = data['h']
	
	#t_length=len(time)-1

	
	#print(len(time))
	
	

	for i in range(len(depth)):
		#print(depth[i], depth_interest)
		if depth[i] == depth_interest:
			
			phi_i[j] = phi[i]
			
			cb_i[j] = cb[i]
			
			cb_initial[j] = 18.9
			
			if phi[i]>0:
				cl_i.append(cl[i])
				time_cl.append(t_int*k/1000)
				temp_cl.append(T[i])
				#cl_i[liq_i] = cl[i]
				#time_cl[liq_i] = t_int*k/1000
				#liq_i=liq_i+1
				
			
			if phi[i]<1:
				cs_i.append(cs[i])
				time_cs.append(t_int*k/1000)
				temp_cs.append(T[i])
				#cs_i[sol_i] = cs[i]
				#time_cs[sol_i] = t_int*k/1000
				#sol_i=sol_i+1
			#cs_i[k] = cs[i]
			
			ol_i[j] = ol[i]
			opx_i[j] = opx[i]
			cpx_i[j] = cpx[i]
			feld_i[j] = feld[i]
			
			T_i[j] = T[i]
			Ts_i[j] = Ts[i]
			Tl_i[j] = Tl[i]
			
			h_i[j] = h[i]
			hs_i[j] = Ts[i]*1100
			hl_i[j] = Tl[i]*1100+550e3
			break

	j=j+1
    ### FIGURE ###



for k in range(len(phi_i)):
	with open("overview_outputs.txt", "w") as file:
		file.write("time\tTl\tTs\tT\tCb\tphi\tol\topx\tcpx\tfeld\n")
		for time_i, tl,ts,t,cb,phi,ol,opx,cpx,feld in zip(time, Tl_i, Ts_i, T_i, cb_i, phi_i, ol_i, opx_i, cpx_i, feld_i):
			file.write(f"{time_i:<10.10f}\t{tl:<10.10f}\t{ts:<10.10f}\t{t:<10.10f}\t{cb:<10.10f}\t{phi:<10.10f}\t{ol:<10.10f}\t{opx:<10.10f}\t{cpx:<10.10f}\t{feld}\n")

for k in range(len(cl_i)):
	with open("cl_outputs.txt", "w") as file:
		file.write("time\ttemp_cl\tcl\n")
		for time_i, T_cl,cl in zip(time_cl, temp_cl, cl_i):
			file.write(f"{time_i:<10.10f}\t{T_cl:<10.10f}\t{cl:<10.10f}\n")


for k in range(len(cs_i)):
	with open("cs_outputs.txt", "w") as file:
		file.write("time\ttemp_cs\tcs\n")
		for time_i, T_cs,cs in zip(time_cs, temp_cs, cs_i):
			file.write(f"{time_i:<10.10f}\t{T_cs:<10.10f}\t{cs:<10.10f}\n")



#t_length=len(time)-1
#for k in range(t_length):
	#plot(k)
        

ax2.clear()

# title of figure

fig.tight_layout()
#ax2.set_title("Temperature" ) 
   

# axis 2
ax2.plot(time,T_i, color='b', label = 'Temperature')
ax2.plot(time,Tl_i, 'k', label = 'Liquidus') 
ax2.plot(time,Ts_i, 'k', linestyle='dashed', label = 'Solidus') 
ax2.set_xlim([time[0], time[-1]])
ax2.set_ylim([600, 1600])
ax2.set_frame_on(True)
ax2.tick_params(direction='in', top=True, right=True, which='major')
ax2.tick_params(direction='in', top=True, right=True, which='minor')
ax2.minorticks_on()
ax2.legend(loc="best")
ax2.set_xlabel("Time (ka)")
ax2.set_ylabel("$Temperature (^oC)")   



#fig.tight_layout()
#numb=depth_interest
#fig.savefig("AA_MF_Temp_depth%skm.svg" %numb)    

ax5.clear()
ax5.plot( cb_initial,T_i, color = '#808080', label = 'Initial bulk MgO') 
ax5.plot( cb_i,T_i, color = 'b', marker='.', linestyle = 'none', label = 'Bulk MgO')  
ax5.set_xlim([0, 57])
ax5.set_ylim([700, 1740])
ax5.set_frame_on(True)
ax5.tick_params(direction='in', top=True, right=True, which='major')
ax5.tick_params(direction='in', top=True, right=True, which='minor')
ax5.minorticks_on()
ax5.legend(loc="best")
ax5.set_ylabel("Temperature (^oC)")
ax5.set_xlabel("MgO %")   



# axis 2
ax3.plot(time,cb_i, color='b', label = 'Bulk MgO')
ax3.plot(time_cs,cs_i, color='k', label = 'Solid MgO') 
ax3.plot(time_cl,cl_i, color='r', label = 'Melt MgO') 
ax3.set_xlim([time[0], time[-1]])  
ax3.set_ylim([0, 60])
ax3.set_frame_on(True)
ax3.tick_params(direction='in', top=True, right=True, which='major')
ax3.tick_params(direction='in', top=True, right=True, which='minor')
ax3.minorticks_on()
ax3.legend(loc="best")
ax3.set_xlabel("Time (ka)")
ax3.set_ylabel("MgO (%)")   

  

# axis 2
ax4.plot(time,ol_i, color='#337538', label = 'ol') 
ax4.plot(time,opx_i, color='#2e2585', label = 'opx') 
ax4.plot(time,cpx_i, color='#c26a77', label = 'cpx') 
ax4.plot(time,feld_i, color='#dcd57d', label = 'fsp') 
ax4.plot(time,phi_i, color='k')
ax4.set_xlim([time[0], time[-1]])
ax4.set_frame_on(True)
ax4.tick_params(direction='in', top=True, right=True, which='major')
ax4.tick_params(direction='in', top=True, right=True, which='minor')
ax4.minorticks_on()
ax4.legend(loc="best")   
ax4.set_ylim([0, 1.01])
ax4.set_xlabel("Time (ka)")
ax4.set_ylabel("%")   




fig.tight_layout()
numb=depth_interest
fig.savefig("AA_depth_data_overTime_@%skm.svg" %numb)    


 
