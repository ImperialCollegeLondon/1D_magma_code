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
fig.set_size_inches(20,15)

# define plot

	

			
dtype1 = np.dtype ([('time','f8'),('evac_from','f8'), ('evac_to', 'f8'),('thick','f8'),('vol','f8'),('comp','f8'), ('nodes_evac','f8'), ('nodes_to','f8')])
filename="Averages_evacuated.txt" 
data = np.loadtxt(filename, dtype=dtype1, skiprows=1)

time = data['time']
thick = data['thick']
comp = data['comp']

evac_to = data['evac_to']
evac_from = data['evac_from']

	
		



	
	
#################################################### Figures #################

ax1.clear()
ax2.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'


# axis 1
ax1.plot(time,thick, 'ko') 

#ax1.legend()
ax1.set_xlim([0, 1800])
#ax1.set_ylim([-50, -25]) 
ax1.set_frame_on(True)
ax1.tick_params(direction='in', top=True, right=True, which='major')
ax1.tick_params(direction='in', top=True, right=True, which='minor')
ax1.minorticks_on()
ax1.set_xlabel("Time (ka)",fontsize=16)
ax1.set_ylabel("Thickness evacuated (km)",fontsize=16)  

# axis 2
ax2.plot(time,comp, 'ko') 
#ax2.legend()
ax2.set_xlim([0, 1800])
#ax2.set_ylim([-50, -25]) 
ax2.set_frame_on(True)
ax2.tick_params(direction='in', top=True, right=True, which='major')
ax2.tick_params(direction='in', top=True, right=True, which='minor')
ax2.minorticks_on()
ax2.set_xlabel("Time (ka)",fontsize=16)
ax2.set_ylabel("Composition evacuated (SiO2 %)",fontsize=16) 


# axis 3
ax3.plot(time,evac_from, 'ko', label='Evacuation depth') 
ax3.plot(time,evac_to, 'ro', label='Intrusion depth') 
ax3.legend(loc="best")
ax3.set_xlim([0, 1800])
#ax3.set_ylim([-50, -25]) 
ax3.set_frame_on(True)
ax3.tick_params(direction='in', top=True, right=True, which='major')
ax3.tick_params(direction='in', top=True, right=True, which='minor')
ax3.minorticks_on()
ax3.set_xlabel("Time (ka)",fontsize=16)
ax3.set_ylabel("Depth (km)",fontsize=16) 
 

fig.tight_layout()

fig.savefig("A_3Phase_evacuations.svg") 





