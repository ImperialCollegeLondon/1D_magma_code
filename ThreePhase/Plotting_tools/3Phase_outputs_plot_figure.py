### import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from celluloid import Camera
import imageio
import ffmpeg
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



FPSN=5


prompt1="Input the number of the output to plot: "
k = int(input(prompt1))



#prompt3="Input the time interval between inputs: "
#t_int = float(input(prompt3))

prompt4="What is the time interval (ka) between outputs? "
t_intka = float(input(prompt4))



time=[None]

numb=k


# initialise figure

fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4)
fig.set_size_inches(20,15.5)

# define plot

	
print(k)

if k==0:

		
	dtype1 = np.dtype ([('depth','f8'),('phi','f8'),('frac_v','f8'),('H','f8'),('T','f8'),
				('Ts','f8'),('Tl','f8'),('cb','f8'),('cl','f8'),
				('cs','f8'),('vb','f8'),('vl','f8'),('vs','f8')
				,('vg','f8')])
				
				
				
else:
				
	dtype1 = np.dtype ([('depth','f8'),('phi','f8'),('frac_v','f8'),('H','f8'),('T','f8'),
				('Ts','f8'),('Tl','f8'),('cb','f8'),('cl','f8'),
				('cs','f8'),('vb','f8'),('vl','f8'),('vs','f8')
				,('vg','f8'),('bulk_rho','f8')])		

		


desired_width=4
if k==0:
	filename="output_0_CELLS.txt" 
else:
	numb=k#"{:0>{}}".format(k)#format(k,desired_width)
	filename="output_%s_CELLS.txt" %numb

data = np.loadtxt(filename, dtype=dtype1, skiprows=2)


# data that will be plotted

depth = data['depth']

phi = data['phi']
frac_v = data['frac_v']

cb = data['cb']
cl = data['cl']
cs = data['cs']

vb = data['vb']
vl = data['vl']
vs = data['vs']
vg = data['vg']

T = data['T']
Tl = data['Tl']
Ts = data['Ts']


	



# calculating time
time = t_intka*k


####################################################here#

### FIGURE ###

# initialising figure, fig


ax1.clear()
ax2.clear()
ax3.clear()
ax4.clear()

# title of figure

fig.tight_layout()
plt.rcParams['font.size']='16'

# axis 1
ax1.plot(frac_v,depth, 'r', label="Volatile fraction") 
ax1.plot(phi,depth, 'b', label="Porosity") 

ax1.legend()
ax1.set_ylim([-50, 0])
ax1.set_xlim([0, 1.01])
#ax1.xaxis.set_minor_locator(MultipleLocator(100))

#ax1.yaxis.set_minor_locator(MultipleLocator(10)) 
ax1.set_xlabel("(-)",fontsize=16)
ax1.set_ylabel("Depth (km)",fontsize=16)


	
# axis 2
ax2.plot(cl,depth, 'r', label="melt") 
ax2.plot(cs,depth, 'k', label="solid")  
ax2.plot(cb,depth, 'b', label="bulk") 

ax2.legend()
ax2.set_ylim([-50, 0])
ax2.set_xlim([40,80])
#ax2.xaxis.set_minor_locator(MultipleLocator(100))
#ax2.set_ylim([bot_plot, top_plot])
#ax2.yaxis.set_minor_locator(MultipleLocator(10)) 
ax2.set_xlabel("$SiO_2$ (%)",fontsize=16)
ax2.set_ylabel("Depth (km)",fontsize=16)


# axis 3
ax3.set_title(f'Time = {time} ka' )
ax3.plot(vl,depth, 'r', label="melt") 
ax3.plot(vs,depth, 'k', label="solid")  
ax3.plot(vg,depth, 'g', label="volatile") 
ax3.plot(vb,depth, 'b', label="bulk")

#ax3.legend(loc="lower left")
ax3.set_ylim([-50, 0])
#ax3.set_ylim([-30, 0])
#ax3.xaxis.set_minor_locator(MultipleLocator(100))
#ax3.set_ylim([bot_plot, top_plot])
#ax3.yaxis.set_minor_locator(MultipleLocator(10)) 
ax3.set_xlabel("$H_2O$ concentration",fontsize=16)
ax3.set_ylabel("Depth (km)",fontsize=16)




# axis 4
ax4.plot(Tl,depth, 'r', label="Liquidus") 
ax4.plot(Ts,depth, 'k', label="Solidus")  
ax4.plot(T,depth, 'b', label="T") 

ax4.legend()
ax4.set_ylim([-50, 0])
ax4.set_xlim([0, 1300])
#ax4.xaxis.set_minor_locator(MultipleLocator(100))
#ax4.set_ylim([bot_plot, top_plot])
#ax4.yaxis.set_minor_locator(MultipleLocator(10)) 
ax4.set_xlabel("Temperature (oC)",fontsize=16)
ax4.set_ylabel("Depth (km)",fontsize=16)	    

fig.tight_layout()

fig.savefig("A_figure_output_3Phase_%s.svg" %numb)


