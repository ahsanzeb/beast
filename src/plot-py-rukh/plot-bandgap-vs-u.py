
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
	for line in lines:
		if (len(line.split())>0):
			a.append( [ float (x) for x in line.split()] )
		else:
			return np.array(a) # only a single spin data required
	return 
#------------------------------------
def getbandgap(striu):
	# filename?
	filename = striu+"/TDOS.OUT";
	# read data for up spin
	data = readdata(filename);
	# some index close to E=0
	i0 = 0
	for [e,dos] in data:
		if (abs(e) <= 0.005):
			break
		else:
			i0 = i0 +1	
	# find top of VB
	for i in range(1000):
		if(data[i0-i,1] >= tol):
			evb = data[i0-i,0];
			break
	
	# find bottom of CB
	for i in range(1000):
		if(data[i0+i,1] >= tol):
			ecb = data[i0+i,0];
			break
	# find bandgap in eV
	gap = (ecb-evb)*27.22
	#print("ecb, evb = ",ecb, evb)

	return gap
#------------------------------------

unames = ["00","01","02","03","04","1",'2','3','4','5'];
uval = [0,0.2,0.4,0.6,0.8,1,2,3,4,5];
tol = 0.01;

out= [];i=0;
for x in unames:
	gap = getbandgap(x)
	out.append([uval[i],gap])
	#print("U, gap = ",uval[i],gap)
	print(uval[i],gap)
	i = i +1;


#exit()

out = np.array(out)


plt.plot(out[:,0],out[:,1],'-X',ms=7,lw=2)


plt.show()


plt.savefig('bandgaps-vs-u.pdf', format='pdf', dpi=200)

