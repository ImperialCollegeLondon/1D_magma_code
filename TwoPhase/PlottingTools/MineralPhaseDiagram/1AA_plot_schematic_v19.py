### import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from celluloid import Camera
import imageio
import ffmpeg
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.patches import RegularPolygon






prompt="What output no. would you like to plot? "
k = int(input(prompt))



### import input file
dtype2 = np.dtype ([('name','str'),('value','f8')])
input = np.loadtxt('1AA_Plot_schematic.txt', dtype=dtype2, delimiter=',')#, skiprows=1)
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



time=[None]*(kEND-kSTART)
### plot initialising
#fig, (ax1, ax2, ax3) = plt.subplots(1,3) #camera method
# camera=Camera(fig) #camera method

#frames1=[]

		


# define plot

	






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

phi = data['phi']
depth = data['z']
ol= data['ol']
opx = data['opx']
cpx = data['cpx']
feld = data['feld']



for i in range(len(depth)):
	

	if depth[i-1]<top_plot and depth[i]>=top_plot:
		top_plot_i=i
		
		
	if depth[i-1]<bot_plot and depth[i]>=bot_plot:
		bot_plot_i=i




ol_sill = ol[bot_plot_i:top_plot_i]
opx_sill = opx[bot_plot_i:top_plot_i]
cpx_sill = cpx[bot_plot_i:top_plot_i]
feld_sill = feld[bot_plot_i:top_plot_i]
phi_sill = phi[bot_plot_i:top_plot_i]
depth_sill = depth[bot_plot_i:top_plot_i]



n_yx = int(len(depth_sill))

n_x=100


data_min = np.zeros((len(depth_sill),n_x), dtype=int)


for i in range(len(depth_sill)):

	ol_sill_int = int(round(ol_sill[i]*100))
	opx_sill_int = int(round(opx_sill[i]*100))
	cpx_sill_int = int(round(cpx_sill[i]*100))
	feld_sill_int = int(round(feld_sill[i]*100))
	phi_sill_int = int(round(phi_sill[i]*100))
	
	if (ol_sill_int+opx_sill_int+cpx_sill_int+feld_sill_int+phi_sill_int)>100:
		tot = ol_sill_int+opx_sill_int+cpx_sill_int+feld_sill_int+phi_sill_int
		A=max(ol_sill_int,opx_sill_int,cpx_sill_int,feld_sill_int, phi_sill_int)
		
		if A==ol_sill_int:
			ol_sill_int=ol_sill_int-(tot-100)
			
		elif A==opx_sill_int:
			opx_sill_int=opx_sill_int-(tot-100)
		
		elif A==cpx_sill_int:
			cpx_sill_int=cpx_sill_int-(tot-100)
			
		elif A==feld_sill_int:
			feld_sill_int=feld_sill_int-(tot-100)
			
		else:
			phi_sill_int=phi_sill_int-(tot-100)
			
			
	if (ol_sill_int+opx_sill_int+cpx_sill_int+feld_sill_int+phi_sill_int)<100:
		tot = ol_sill_int+opx_sill_int+cpx_sill_int+feld_sill_int+phi_sill_int
		A=max(ol_sill_int,opx_sill_int,cpx_sill_int,feld_sill_int, phi_sill_int)
		
		if A==ol_sill_int:
			ol_sill_int=ol_sill_int+(100-tot)
			
		elif A==opx_sill_int:
			opx_sill_int=opx_sill_int+(100-tot)
		
		elif A==cpx_sill_int:
			cpx_sill_int=cpx_sill_int+(100-tot)
			
		elif A==feld_sill_int:
			feld_sill_int=feld_sill_int+(100-tot)
			
		else:
			phi_sill_int=phi_sill_int+(100-tot)
			
					
	
	#print(ol_sill_int+opx_sill_int+cpx_sill_int+feld_sill_int+phi_sill_int, tot, A)
	
	
	counts = {1: ol_sill_int, 2: opx_sill_int, 3: cpx_sill_int, 4:feld_sill_int, 5: phi_sill_int}
	
	arrays = [np.full(count, value) for value, count in counts.items()]
	
	n_z = np.concatenate(arrays)
	
	if phi_sill_int>0:
		np.random.shuffle(n_z)
	else:
		np.random.seed(i)
		np.random.shuffle(n_z)
	
	data_min[i] = n_z
	
	#print(data_min[i])
	
	



#print(data_min)


# calculating time
time[k] = t_int*k/1000


### FIGURE ###

# initialising figure, fig

print('Loading...')

fig, ax = plt.subplots(figsize=(8, 8))

def hex_heatmap(n_y=50, hex_size=1):
    n_x = 100  # fixed as per your requirement
    

    # Spacing between hex centers
    dx = 3/2 * hex_size
    dy = np.sqrt(3) * hex_size / 2

    # Create data (for example, z = some function of x and y)
    z_values = data_min#np.random.rand(n_y, n_x)  # you can replace this with real data
    
    color_map = {1:'#337538' , 2:'#2e2585', 3:'#c26a77', 4:'#dcd57d', 5:'#DDDDDD'}

    # Draw hexagons
    for row in range(n_y):
        for col in range(n_x):
            # Offset every other row
            x = col * dx
            y = row * dy * 2
            if col % 2 == 1:
                y += dy

            color = color_map[z_values[row, col]]  # map z  color
            hex = RegularPolygon((x, y), numVertices=6, radius=hex_size,
                                 orientation=np.radians(30),
                                 facecolor=color, edgecolor=color, lw=0.5)
            ax.add_patch(hex)

    ax.set_aspect('equal')
    ax.autoscale_view()
    ax.axis('off')
    plt.title(f'Hex Heatmap ({n_x} x {n_y})', fontsize=14)
    plt.show()
   


    




# title of figure
print('Once pop up appears, close it to save!')
hex_heatmap(n_y=n_yx)
print('Saving figure...')
fig.savefig("Figure_of_output_Schematic_%s_paper.svg" %numb)
print('Figure saved')



