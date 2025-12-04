### import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from celluloid import Camera
import imageio
import ffmpeg
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)









# initialise figure

fig, (ax1,ax2,ax3) = plt.subplots(3,1)
fig.set_size_inches(30,20)

# define plot

	

			
dtype1 = np.dtype ([('time','f8'),('mass_cons','f8'),('mass_res','f8')])
filename="mass_cons.txt" 
data = np.loadtxt(filename, dtype=dtype1, skiprows=1)

dtype1 = np.dtype ([('time','f8'),('precision','f8'),('res','f8')])
filename="In_Depth_break_data.txt" 
data1 = np.loadtxt(filename, dtype=dtype1, skiprows=1)


dtype1 = np.dtype ([('time','f8'),('percent','f8'),('breaks','f8'),('timesteps','f8')])
filename="Break_of_tolerance.txt" 
data2 = np.loadtxt(filename, dtype=dtype1, skiprows=1)



#phi_top = np.array([data['phi_top_km']])


time_mass = data['time']
mass_percent = data['mass_cons']

time_break = data1['time']
precision = data1['precision']
residual = data1['res']

time_tol = data2['time']
percent = data2['percent']

	
		



	
	
#################################################### Figures #################

ax1.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'


# axis 1
ax1.plot(time_mass,mass_percent, 'kx-') 

#ax1.legend()
#ax1.set_xlim([0, 4])
#ax1.set_ylim([-50, -25]) 
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Percentage of mass conserved (%)",fontsize=16)  

# axis 2
ax2.plot(time_tol,percent, 'kx-') 
#ax2.legend()
#ax2.set_xlim([0, 4])
#ax2.set_ylim([-50, -25]) 
ax2.set_frame_on(True)
ax2.tick_params(direction='in', top=True, right=True, which='major')
ax2.tick_params(direction='in', top=True, right=True, which='minor')
ax2.minorticks_on()
ax2.set_xlabel("Time (ka)",fontsize=16)
ax2.set_ylabel("Percentage of timesteps failed to conserve (%)",fontsize=16) 

# axis 3
ax3.plot(time_break,precision, 'rx-', label='Precision') 
ax3.plot(time_break,residual, 'kx-', label='Residual') 

ax3.legend(loc="best")
ax3.set_frame_on(True)
ax3.tick_params(direction='in', top=True, right=True, which='major')
ax3.tick_params(direction='in', top=True, right=True, which='minor')
ax3.minorticks_on()
#ax3.set_xlim([0, 4])
#ax3.set_ylim([-50, -25]) 
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("(-)",fontsize=16) 

fig.tight_layout()

fig.savefig("A_3Phase_errors.svg") 





