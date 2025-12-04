### import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from celluloid import Camera
import imageio
import ffmpeg
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


prompt1="What phase diagram was used? Enter: 1 = Hu et al. 22, JPet. 2 = Solid phase diagram. 3 = Layered intrusion phase diagram. "
PD = int(input(prompt1))

prompt2="Would you like to plot by depth (answer 0) or by thickness (1)? "
depth_flag = int(input(prompt2))

if depth_flag == 0:
	depth_corrector = 0
else:
	prompt3 = "What was the intrusion depth? "
	depth_corrector = float(input(prompt3))

### import input file
dtype2 = np.dtype ([('name','str'),('value','f8')])
input = np.loadtxt('1AA_Plot_video_output_1.txt', dtype=dtype2, delimiter=',')#, skiprows=1)
input_value=input['value']
# start file for output
kSTART=int(input_value[0])
# end file for output
kEND=int(input_value[1])
# time interval between outputs
t_int=input_value[2]
# top depth of the plots
top_plot=input_value[3]
# bottom of the plots
bot_plot=input_value[4]

if depth_flag==1:
	top_plot=abs(top_plot+depth_corrector)*1000
	bot_plot=abs(bot_plot+depth_corrector)*1000

if PD==3:
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
	# contribution outputted?
	cont_out = input_value[10]
else:
	#min sio2
	min_sio2=input_value[5]
	#max sio2
	max_sio2=input_value[6]
	# frames per second
	FPSN = input_value[7]
	# contribution outputted?
	cont_out = input_value[8]

time=[None]*(kEND-kSTART)
### plot initialising
#fig, (ax1, ax2, ax3) = plt.subplots(1,3) #camera method
# camera=Camera(fig) #camera method

#frames1=[]
if PD==3:
	if cont_out==1:
		fig, (ax1,ax4, ax5) = plt.subplots(1,3)
	else:
		fig, (ax1, ax3, ax4) = plt.subplots(1,3)
else:
	if cont_out==1:
		fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1,5)
	else:
		fig, (ax1, ax2, ax3) = plt.subplots(1,3)
		
fig.set_size_inches(20,10.5)

# define plot
def plot(k):
	
	print(k)
	if PD==1 or PD==2:

		if cont_out==1:
			if k==0:
				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
							('Tl','f8'),('cb','f8'),('cb_mg','f8'),
							('cs_mg','f8'),('cl_mg','f8'),('h_mf','f8'),('c_mf','f8')
							,('r_mf','f8'),('h_mf_t','f8'),('c_mf_t','f8')
							,('r_mf_t','f8'),('c_cb','f8'),('r_cb','f8'),('c_cb_t','f8')
							,('r_cb_t','f8')])
				

			else:

				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
							('Tl','f8'),('cb','f8'),('cb_mg','f8'),
							('cs_mg','f8'),('cl_mg','f8'),('dens','f8'),('h_mf','f8'),('c_mf','f8')
							,('r_mf','f8'),('h_mf_t','f8'),('c_mf_t','f8')
							,('r_mf_t','f8'),('c_cb','f8'),('r_cb','f8'),('c_cb_t','f8')
							,('r_cb_t','f8')])
							
		else:
			if k==0:
				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
							('Tl','f8'),('cb','f8'),('cb_mg','f8'),
							('cs_mg','f8'),('cl_mg','f8')])
				

			else:

				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
							('Tl','f8'),('cb','f8'),('cb_mg','f8'),
							('cs_mg','f8'),('cl_mg','f8'),('dens','f8')])



	elif PD==3:

		if cont_out==1:
			if k==0:
		
				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
						('Tl','f8'),('cb','f8'),('cb_mg','f8'),
						('cs_mg','f8'),('cl_mg','f8'),('ol','f8'),('opx','f8'),
						('cpx','f8'),('feld','f8'),('h_mf','f8'),('c_mf','f8')
						,('r_mf','f8'),('h_mf_t','f8'),('c_mf_t','f8')
						,('r_mf_t','f8'),('c_cb','f8'),('r_cb','f8'),('c_cb_t','f8')
						,('r_cb_t','f8')])
			

			else:

	
				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
					('Tl','f8'),('cb','f8'),('cb_mg','f8'),
					('cs_mg','f8'),('cl_mg','f8'),('ol','f8'),('opx','f8'),
					('cpx','f8'),('feld','f8'),('dens','f8'),('h_mf','f8'),('c_mf','f8')
					,('r_mf','f8'),('r_mf_dum','f8'),('h_mf_t','f8'),('c_mf_t','f8')
					,('r_mf_t','f8'),('r_mf_t_dum','f8'),('c_cb','f8'),
					('r_cb','f8'),('r_cb_dum','f8'),('c_cb_t','f8'),('r_cb_t','f8')
					,('r_cb_t_dum','f8')])
						
		else:
			if k==0:
		
				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
						('Tl','f8'),('cb','f8'),('cb_mg','f8'),
						('cs_mg','f8'),('cl_mg','f8'),('ol','f8'),('opx','f8'),
						('cpx','f8'),('feld','f8')])
			

			else:

		
				dtype1 = np.dtype ([('z','f8'),('phi','f8'),('h','f8'),('T','f8'),('Ts','f8'),
							('Tl','f8'),('cb','f8'),('cb_mg','f8'),
						('cs_mg','f8'),('cl_mg','f8'),('ol','f8'),('opx','f8'),
						('cpx','f8'),('feld','f8'),('dens','f8')])
						
        

	desired_width=4
	if k==0:
		filename="output_0_CELLS.txt" 
	else:
		numb=k#"{:0>{}}".format(k)#format(k,desired_width)
		filename="output_%s_CELLS.txt" %numb

	data = np.loadtxt(filename, dtype=dtype1, skiprows=2)

	if PD==3:
	
		Ts = data['Ts']
		Tl = data['Tl']
		
		depth_initial = data['z']
	
		depth = [0]*len(depth_initial)
		
		hs = Ts*1100
		hl = Tl*1100+550e3
		hin_cl = [0]*len(Tl)
		hin_cs = [0]*len(Tl)
		og_cb = [0]*len(Tl)
		
		for i in range(len(Tl)):
			hin_cl[i]=0.6*100
			hin_cs[i]=100-hin_cl[i]
			og_cb[i] = 18.9
		
			if depth_flag == 0:
				depth[i] = depth_initial[i]+depth_corrector
			else:
				depth[i] = -(depth_initial[i]+depth_corrector)*1000

		if cont_out==1:
			comp_cb_dt = data['c_cb_t']
			react_cb_dt = data['r_cb_t']
			phi = data['phi']
			cl = data['cl_mg']
			cs = data['cs_mg']
			
			comp_phi = [0]*len(cl)
			react_phi = [0]*len(cl)
			clcs_phi = [0]*len(cl)
			
			for i in range(len(cl)):
				if phi[i] == 0:
					comp_phi[i] = 0
					react_phi[i] = 0
					clcs_phi[i] = 0
				else:
					comp_phi[i] = comp_cb_dt[i]/phi[i]
					react_phi[i] = react_cb_dt[i]/phi[i]
					clcs_phi[i] = (cs[i] - cl[i])/phi[i]  
					
		
		
		# finding depths of magma and mush
	
		phi_plot = [0, 100]
		
		top_mush_depth=[]
		top_magma_depth=[]
		bot_mush_depth=[]
		bot_magma_depth=[]
		
		for i in range(len(depth)):
			if phi[i]>=0.6 and phi[i-1]<0.6:
				top_magma_depth.append(depth[i])
			
			if phi[i]>=0.6 and phi[i+1]<0.6:
				bot_magma_depth.append(depth[i])
				
			if phi[i]>0 and phi[i]<=0.6 and (phi[i-1]>0.6):
				top_mush_depth.append(depth[i-1])
				
			if phi[i]>0 and phi[i]<=0.6 and (phi[i+1]>0.6):
				bot_mush_depth.append(depth[i+1])
				
			if phi[i]>0 and phi[i]<=0.6 and (phi[i-1]==0):# or phi[i-1]>0.6):
				top_mush_depth.append(depth[i])
				
			if phi[i]>0 and phi[i]<=0.6 and (phi[i+1]==0 ):#or phi[i+1]>0.6):
				bot_mush_depth.append(depth[i])

		magma_plot = len(top_magma_depth)
		mush_plot = len(top_mush_depth)
    
    

    # calculating time
	time[k] = t_int*k/1000


    ### FIGURE ###

        # initialising figure, fig


	ax1.clear()
	ax4.clear()
	ax5.clear()
# title of figure

	fig.tight_layout()
	plt.rcParams['font.size']='16'
	if PD==3:
		# axis 1
		ax1.plot(data['T'],depth, color='b', label="Temperature") # temperature
		ax1.plot(data['Ts'],depth, color='k',linestyle='dashed',label="Solidus") # solidus
		ax1.plot(data['Tl'],depth, color='k',label="Liquidus") # liquidus
		#ax1.legend(loc="lower left")
		ax1.set_xlim([0, 2000])
		ax1.xaxis.set_minor_locator(MultipleLocator(100))
		ax1.set_ylim([bot_plot, top_plot])
		ax1.yaxis.set_minor_locator(MultipleLocator(10)) 
		ax1.set_xlabel("Temperature ($^o$C)",fontsize=16)
		ax1.set_ylabel("Depth (km)",fontsize=16)
		
		    

		 

		# axis 3
		ax4.set_title(f'Time = {time[k]} ka' )#% time[k] ) 
		#ax4.plot((data['ol']+data['opx']+data['cpx']+data['feld'])*100,data['z'],color='k', label='Total')
		ax4.plot(100-data['phi']*100,depth, color='k',linewidth=2) # melt fraction
	
		ax4.fill_between((data['ol']+data['opx']+data['cpx']+data['feld'])*100,depth,color='#dcd57d',  label='feld')
		ax4.fill_between((data['ol']+data['opx']+data['cpx'])*100,depth,color='#c26a77',  label='cpx')
		ax4.fill_between((data['ol']*100+data['opx']*100),depth,color='#2e2585',  label='opx')
		ax4.fill_between((data['ol']*100),depth, color='#337538', label='ol')
		ax4.fill_between(100-data['phi']*100,depth, color='#DDDDDD', label='melt')
		#ax4.plot(hin_cl,depth, color='k', linestyle='dashed')
		
		for i in range(magma_plot):
			ax4.plot(phi_plot, [top_magma_depth[i], top_magma_depth[i]], color='#ff4500', linestyle='dashed')
			ax4.plot(phi_plot, [bot_magma_depth[i], bot_magma_depth[i]], color='#ff4500', linestyle='dashed')
			ax4.text(30, (top_magma_depth[i]+bot_magma_depth[i])/2, 'magma', color='#ff4500')
			
		for i in range(mush_plot):
			ax4.plot(phi_plot, [top_mush_depth[i], top_mush_depth[i]], color='#ff4500', linestyle='dashed')
			ax4.plot(phi_plot, [bot_mush_depth[i], bot_mush_depth[i]], color='#ff4500', linestyle='dashed')
			
			if abs((top_mush_depth[i]-bot_mush_depth[i]))>20:
				ax4.text(30, (top_mush_depth[i]+bot_mush_depth[i])/2, 'mush', color='#ff4500')
			else:
				ax4.text(30, top_mush_depth[i], 'mush', color='#ff4500')
				
				
		#ax4.legend(loc="lower left")
		ax4.set_xlim([0, 100])
		if k==0:
			ax44 = ax4.twiny()
			ax44.set_xticks([0, 20, 40, 60, 80, 100])
			ax44.set_xticklabels([100, 80, 60, 40, 20, 0])
			ax44.set_xlabel("Melt fraction (%)",fontsize=16)
			
			
			
		#ax4.xaxis.set_minor_locator(MultipleLocator(0.1))    
		ax4.set_ylim([bot_plot, top_plot])
		ax4.yaxis.set_minor_locator(MultipleLocator(10))
		ax4.set_xlabel("Solid fraction (%)",fontsize=16)
		ax4.set_ylabel("Depth (km)",fontsize=16)  
		
		
		

		# axis 4
		ax5.plot(og_cb,depth, color='#808080',label='Initial bulk') # initial bulk mgo
		ax5.plot(data['cs_mg'],depth, color='k', label='Solid') #solid
		ax5.plot(data['cl_mg'],depth, color='r', label='Melt') #melt mgo
		ax5.plot(data['cb_mg'],depth, color='b',label='Bulk') #bulk mgo
		
		ax5.set_xlim([min_mgo, max_mgo])
		ax5.set_ylim([bot_plot, top_plot])
		#ax5.set_xticks([45,55,65,75])
		#ax5.set_xticklabels(['45','55','65','75'])
		ax5.xaxis.set_minor_locator(MultipleLocator(5))
		ax5.yaxis.set_minor_locator(MultipleLocator(10))
		#ax5.legend(loc="lower left")
		ax5.set_xlim([0, 60])
		ax5.set_xlabel("MgO (%)",fontsize=16)
		ax5.set_ylabel("Depth (km)",fontsize=16) 


	
			

	else:
		# axis 1
		ax1.plot(data['T'],data['z'], color='b', label="Temperature") # temperature
		ax1.plot(data['Ts'],data['z'], color='r',linestyle='dashed',label="Solidus") # solidus
		ax1.plot(data['Tl'],data['z'], color='r',label="Liquidus") # liquidus
		#ax1.legend(loc="lower left")
		ax1.set_xlim([0, 2000])
		ax1.xaxis.set_minor_locator(MultipleLocator(100))
		ax1.set_ylim([bot_plot, top_plot])
		ax1.set_xlabel("Temperature ($^o$C)")
		ax1.set_ylabel("Depth (km)")    

		# axis 2
		#ax2.set_title(f'Time = {time[k]} ka' )#% time[k] ) 
		ax2.plot(data['phi']*100,data['z'], color='b',label='Melt') # melt fraction
		ax2.set_xlim([0, 100])
		#ax2.xaxis.set_minor_locator(MultipleLocator(0.1))    
		ax2.set_ylim([bot_plot, top_plot])
		ax2.set_xlabel("Melt fraction (%)")
		ax2.set_ylabel("Depth (km)")   

		# axis 3
		ax3.plot(data['cs_mg'],data['z'], color='g', label='Solid') #solid
		ax3.plot(data['cl_mg'],data['z'], color='r', label='Melt') #melt mgo
		ax3.plot(data['cb_mg'],data['z'], color='b',label='Bulk') #bulk mgo
		ax3.set_xlim([min_sio2, max_sio2])
		ax3.set_ylim([bot_plot, top_plot])
		#ax3.set_xticks([45,55,65,75])
		#ax3.set_xticklabels(['45','55','65','75'])
		ax3.xaxis.set_minor_locator(MultipleLocator(5))
		ax3.legend(loc="lower left")
		ax3.set_xlabel("Sio2 (%)")
		ax3.set_ylabel("Depth (km)") 

		if cont_out==1:
			# axis 4
			ax4.plot(data['h_mf'],data['z'], color='b', label='Heating/Cooling') 
			ax4.plot(data['r_mf'],data['z'], color='g', label='Reactive flow') 
			ax4.plot(data['c_mf'],data['z'], color='r', label='Compaction') 
			#ax4.set_xlim([0, 60])
			ax4.set_ylim([bot_plot, top_plot])
			#ax4.xaxis.set_minor_locator(MultipleLocator(5))
			ax4.legend(loc="lower left")
			ax4.set_xlabel("$\sum \Delta \phi$ (-)")
			ax4.set_ylabel("Depth (km)") 


			# axis 5
			ax5.plot(data['r_cb'],data['z'], color='g', label='Reactive flow') 
			ax5.plot(data['c_cb'],data['z'], color='r', label='Compaction') 
			#ax5.set_xlim([0, 60])
			ax5.set_ylim([bot_plot, top_plot])
			#ax5.xaxis.set_minor_locator(MultipleLocator(5))
			ax5.legend(loc="lower left")
			ax5.set_xlabel("Sum of bulk SiO2 change (%)")
			ax5.set_ylabel("Depth (km)")

		ax44.clear()
	fig.tight_layout()
    	


ani = manimation.FuncAnimation(fig,plot,frames=kEND-kSTART)
FFwrite=manimation.FFMpegWriter(fps=FPSN)   
if cont_out==1:
	ani.save("A_video_output_fill_cont3p_opp_NAME.mp4", writer=FFwrite) 
else:
	ani.save("A_video_output_fill6p_NAME.mp4", writer=FFwrite) 
#animation = camera.animate() #camera method
#animation.save('animation.gif', writer='PillowWriter', fps=2) #camera method

