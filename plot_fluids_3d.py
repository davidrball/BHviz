import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#first check for existence of saved numpy arrays
file_path = "np_arr_data/test_"

try:
    ne_array=np.load(file_path+"ne.npy")
    te_array=np.load(file_path+"te.npy")
    b_array=np.load(file_path+"b.npy")
    params=np.load(file_path+"params.npy")
    phys_size = params[0]
    mylen = params[1]
    print('loading from numpy files')
except:
    print('loading from CSV')
    fp = open("CSV_out/3D_test_out.txt" , "r")
    phys_size = float(fp.readline().split(": ")[1][:-3]) #get physical info on first line
    print(phys_size)
    mylen = int(fp.readline().split(": ")[1])
    print(mylen)
    b_array = np.zeros((mylen,mylen,mylen))
    ne_array = np.zeros((mylen,mylen,mylen))
    te_array = np.zeros((mylen,mylen,mylen))
    fp.readline()
    for line in fp.readlines():
        strline = line.split(" , ")
        i = int(strline[0])
        j = int(strline[1])
        k = int(strline[2])
        b = float(strline[3])
        ne = float(strline[4])
        te = float(strline[5])
        b_array[j][i] += b 
    #print(b)
        if ne > 0:
            ne_array[k][j][i] += ne
        te_array[k][j][i] += te 
    fp.close()




extent_list = [-phys_size/2., phys_size/2., -phys_size/2., phys_size/2.]
xaxlabel="$x \; (GM/c^2)$"
zaxlabel = "$z \; (GM/c^2)$"
yaxlabel="$y \; (GM/c^2)$"
half = int(mylen/2.)
constx_slice = ne_array[:,:,half]
consty_slice = ne_array[:,half,:]
constz_slice = ne_array[half,:,:]
mymax=max(np.max(constx_slice), np.max(consty_slice), np.max(constz_slice))

mymax = np.log10(mymax)
fig, (ax0, ax1, ax2) = plt.subplots(1,3,sharey=True)
ax0.imshow(np.log10(constx_slice), origin='lower',vmin=-2,vmax=mymax,extent=extent_list)
ax1.imshow(np.log10(consty_slice),origin='lower',vmin=-2,vmax=mymax,extent=extent_list)
ax2.imshow(np.log10(constz_slice),origin='lower',vmin=-2, vmax=mymax,extent=extent_list)
#ax0.set_xlabel(xaxlabel)
plt.savefig("fluid_plots/3dtest.png",dpi=300,bbox_inches='tight')

