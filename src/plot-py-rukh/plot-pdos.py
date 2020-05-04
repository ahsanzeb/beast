
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.axes as axx
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import MaxNLocator,FixedLocator,IndexLocator
#matplotlib.use('MacOSX')


#------------------------------------
def readdata(filename):
	a=[]; 
	fin = open(filename,'r')
	lines = fin.readlines();
	for line in lines[1:]: # skip the first line that has comments
		a.append( [ float (x) for x in line.split()] )
	return np.array(a)
#------------------------------------
def getxy(iu,isp,spin):
	ialist = [1,2];
	# select range of columns
	if (spin==1):
		rang = range(2,18);
	elif(spin==-1):
		rang = range(18,33);
	# loop over atoms and sum over columns in rang
	for ia in ialist:
		# filename?
		filename = str(iu)+"/PDOS_S0" + str(isp) + "_A000"+str(ia) + ".OUT.columns"
		# read data 
		data = readdata(filename);
		if (ia == ialist[0]): # initialise output array
			y = np.zeros(data.shape[0])
		y += np.sum(data[:,rang], axis=1) # isp=1: Ni d orbitals
	x0 = data[:,0];
	return x0,y
#------------------------------------


#------------ make figure and a grid of axes -----------------
fig = plt.figure(1,figsize=(10,12), dpi=80)

plt.axes([0.0, 0.0, 0.85, 0.85])
plt.title('Impulse response')
newax=plt.axes()
newax.set_frame_on(False)
newax.patch.set_visible(False)

#newax.patch.set_alpha(0.0001)
newax.set_xticks([])
newax.set_yticks([])
newax.set_xlabel(r'$E (eV)$', fontsize=20, fontweight='bold')
newax.xaxis.set_label_coords(0.6, 0)

newax.set_ylabel(r'$DOS (states/eV/unit cell)$', fontsize=20, fontweight='bold')



grid = AxesGrid(fig, (0.2,0.1,0.75,0.75),  # similar to subplot(111)
                    nrows_ncols=(3, 2),
                    axes_pad=0.0,
                    share_all=True,
                    label_mode="L"
)







# axes ranges, and ticks etc
xi,xf = -4,4
yi,yf = -10,10

xtcs=range(xi+1,xf,1);
ytcs=range(yi+2,yf-1,4);

#--------------------loop over u values and ploting ----------------
for iu in range(1,6,1):

	# get axes
	ax = grid[iu-1]; #.get_axes( ) # select axes from the grid of axes to plot on
	for spin in [1,-1]:
		# obtain data 
		x,y1 = getxy(iu,1,spin) # Ni; 
		x,y2 = getxy(iu,2,spin) # O
		x = x*27.2122 # Hartree to eV
		y1 = y1/27.2122;# per Hartre to per eV
		y2 = y2/27.2122;
		ax.set_xlim(xi,xf)
		ax.set_ylim(yi,yf)
		ax.set_aspect(0.7*(xf-xi)/(yf-yi))
		if iu==1 and spin==1:
			ax.plot(x,y1+y2,lw=3,color="green",label="TDOS")
			ax.plot(x,y1,lw=3,color="red",label="Ni")
			ax.plot(x,y2,lw=3,color="blue",label="O")
		else:
			ax.plot(x,y1+y2,lw=3,color="green")
			ax.plot(x,y1,lw=3,color="red")
			ax.plot(x,y2,lw=3,color="blue")

		# labels , U values
		ax.text(0,8,"U = ${"+str(iu)+"}$ eV",color='k',fontsize=16)


	if ( iu%2 == 1 ):
		ax.yaxis.set_major_locator(FixedLocator(ytcs))
		ax.set_yticklabels(ytcs,fontsize=18)
	else:
		ax.yaxis.set_major_locator(FixedLocator(ytcs))
		ax.set_yticklabels([],fontsize=18)
	# x labels and ticks
	if ( iu >3 ):
		ax.xaxis.set_major_locator(FixedLocator(xtcs))
		ax.set_xticklabels(xtcs,fontsize=18)
	else:
		ax.xaxis.set_major_locator(FixedLocator(xtcs))
		ax.set_xticklabels([],fontsize=18)


ax = grid[5]; #.get_axes( ) # select axes from the grid of axes to plot on
ax.xaxis.set_major_locator(FixedLocator(xtcs))
ax.set_xticklabels(xtcs,fontsize=18)


handles, labels = grid[0].get_legend_handles_labels()
newax.legend(handles, labels, frameon=False, ncol=3,fontsize=18,loc=(0.3,.9))


plt.savefig('pdos-combined.pdf', format='pdf', dpi=200)

plt.show()
