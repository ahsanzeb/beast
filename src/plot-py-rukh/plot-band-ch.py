
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.axes as axx
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import MaxNLocator,FixedLocator,IndexLocator
#matplotlib.use('MacOSX')
from matplotlib.collections import LineCollection

#------------------------------------
def readdata(filename):
	a=[]; b=[]; mkblock=True
	fin = open(filename,'r')
	lines = fin.readlines();
	for line in lines:
		if (len(line.split()) == 0):
			if(mkblock):
				b.append(a);
				a=[];
				mkblock=False
		else:
			mkblock=True
			a.append( [ float (x) for x in line.split()] )
			
	return np.array(b) #! iband, irow, icol
#------------------------------------
def getxy(iu,isp,ia):
	# filename?
	filename = str(iu)+"/BAND_S0" + str(isp) + "_A000"+str(ia) + ".OUT"
	# read data 
	data = readdata(filename);
	print(data.shape)
	n,m,p = data.shape
	return int(n/2), data[0:int(n/2),:,0:3] # return only first half number of bands and first 3 columns
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
yi,yf = -3,3
ytcs=np.arange(yi+2.5,yf-0.1,0.5);

#--------------------loop over u values and ploting ----------------
x = 2; # scale all weights to see a bigger effect

col = ['red','blue'];
lbl = ["Ni","O"];
for iu in range(1,6,1):
	# get axes
	ax = grid[iu-1]; #.get_axes( ) # select axes from the grid of axes to plot on
	for isp in [1,2]:
		for ia in [1,2]:
			# obtain data 
			n,dat = getxy(iu,isp,ia) # Ni; 
			dat[:,:,1] = dat[:,:,1]*27.2122 # Hartree to eV
			xi = dat[0,0,0]; # second index for rows: 0th row, first val of k, min
			xf = dat[0,-1,0]; # second index for rows: -1th row, last val of k, max
			ax.set_xlim(xi,xf)
			ax.set_ylim(yi,yf)
			ax.set_aspect(0.7*(xf-xi)/(yf-yi))

			xtcs=np.arange(xi+1,xf,1);

			for ib in range(n):
				lwidths=dat[ib,:-1, 2]*x # scale weight
				points = dat[ib,:,0:2].reshape(-1, 1, 2)
				segments = np.concatenate([points[:-1], points[1:]], axis=1)
				lc = LineCollection(segments, linewidths=lwidths,color=col[isp-1])
				ax.add_collection(lc)

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
