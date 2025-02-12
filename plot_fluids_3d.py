import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#first check for existence of saved numpy arrays
pathext = "arr50_M200"

file_path = "np_arr_dataasdf/"+pathext+"_"

try:
    ne_array=np.load(file_path+"ne.npy")
    te_array=np.load(file_path+"te.npy")
    b_array=np.load(file_path+"b.npy")
    params=np.load(file_path+"params.npy")
    phys_size = params[0]
    mylen = params[1]
    print('loading from numpy files')
    tup=100
    tlow=0
    te_array=np.clip(te_array,tlow,tup)
except:
    print('loading from CSV')
    fp = open("CSV_out/3D_arr50_M50_nohead_nonan.csv" , "r")
    #phys_size = float(fp.readline().split(": ")[1][:-3]) #get physical info on first line
    #print(phys_size)
    #mylen = int(fp.readline().split(": ")[1])
    #print(mylen)
    mylen=50
    b_array = np.zeros((mylen,mylen,mylen))
    ne_array = np.zeros((mylen,mylen,mylen))
    te_array = np.zeros((mylen,mylen,mylen))
    params_array = np.zeros(2)
    params_array[0]=phys_size
    params_array[1]=mylen
    fp.readline()
    for line in fp.readlines():
        strline = line.split(",")
        i = int(strline[0])
        j = int(strline[1])
        k = int(strline[2])
        b = float(strline[3])
        ne = float(strline[4])
        te = float(strline[5])
        b_array[j][i] += b 
    #print(b)
        if True:
            ne_array[k][j][i] += ne
        te_array[k][j][i] += te 
    
    
    #clean up the data a bit
    #ne_array = np.nan_to_num(ne_array)
    #te_array = np.nan_to_num(te_array)
    #b_array = np.nan_to_num(b_array)
    
    tup=100
    tlow=0
    #te_array = np.clip(te_array,tlow,tup)
    
    
    
    
    #np.save(file_path+"ne.npy",ne_array)
    #np.save(file_path+"te.npy",te_array)
    #np.save(file_path+"b.npy",b_array)
    #np.save(file_path+"params.npy",params_array)
    fp.close()



print(np.min(te_array),np.max(te_array))

extent_list = [-phys_size/2., phys_size/2., -phys_size/2., phys_size/2.]
xaxlabel="$x \; (GM/c^2)$"
zaxlabel = "$z \; (GM/c^2)$"
yaxlabel="$y \; (GM/c^2)$"
half = int(mylen/2.)
constx_slice = te_array[:,:,half]
consty_slice = te_array[:,half,:]
constz_slice = te_array[half,:,:]
mymax=max(np.max(constx_slice), np.max(consty_slice), np.max(constz_slice))

mymax = np.log10(mymax)

fig, (ax0, ax1, ax2) = plt.subplots(1,3,sharey=True)
im0=ax0.imshow(np.log10(constx_slice), origin='lower',extent=extent_list)
im1=ax1.imshow(np.log10(consty_slice),origin='lower',extent=extent_list)
im2=ax2.imshow(np.log10(constz_slice),origin='lower',extent=extent_list)
#ax0.set_xlabel(xaxlabel)
fig.colorbar(im0)
#plt.savefig("fluid_plots/"+pathext+".png",dpi=300,bbox_inches='tight')
#fig.close()
#plt.close()

'''
te_array = np.nan_to_num(te_array)
flatarr = te_array.flatten()
plt.hist(flatarr,bins=100)
plt.yscale('log')
plt.savefig('fluid_plots/tmp_hist.png')
'''
