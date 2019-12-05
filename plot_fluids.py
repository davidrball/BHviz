import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

fp = open("CSV_out/test_out.txt" , "r")
fp.readline()



mysize = 0 # get dimensions
for line in fp.readlines():
    mysize += 1

mylen = np.sqrt(mysize)

if mylen / int(mylen) != 1:
    print("grid dimension not square, terminating plotting")
    sys.exit(0)

fp.seek(0)
fp.readline()


mylen = int(mylen)
b_array = np.zeros((mylen,mylen))
ne_array = np.zeros((mylen,mylen))
te_array = np.zeros((mylen,mylen))

xmax = 0
mycount = 0
for line in fp.readlines():
    strline = line.split(" , ")
    i = int(strline[0])
    j = int(strline[1])
    x = float(strline[2])
    z = float(strline[3])
    b = float(strline[4])
    ne = float(strline[5])
    te = float(strline[6])
    ur = float(strline[7])
    uphi = float(strline[8])
    b_array[j][i] += b 
    #print(b)
    if ne > 0:
        ne_array[j][i] += ne
    te_array[j][i] += te 
    mycount += 1

    if mycount == mylen:

        #once we're at one edge of the domain (which happens to be along z in our particular order)
        xmax=round(z)

extent_list = [0,xmax, 0, xmax]
fp.close()

xaxlabel="$x \; (GM/c^2)$"
yaxlabel = "$y \; (GM/c^2)$"

fig, (ax0, ax1, ax2) = plt.subplots(1,3,sharey=True)
ax0.imshow(np.log10(ne_array), origin='lower',extent=extent_list)
ax1.imshow(np.log10(te_array),origin='lower', extent= extent_list,cmap="plasma")
ax2.imshow(np.log10(b_array),origin='lower',extent= extent_list,cmap ="Blues")
ax0.set_title('Density')
ax1.set_title('Temperature')
ax2.set_title('$B^2$')

ax0.set_xlabel(xaxlabel)
ax1.set_xlabel(xaxlabel)
ax2.set_xlabel(xaxlabel)
ax0.set_ylabel(yaxlabel)


plt.savefig("fluid_plots/test.png",dpi=300,bbox_inches='tight')
