

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#how many random seeds to initialized
seednum=5
fileout = "blines_csv/test_blines.csv"


fp = open("blines_csv/test_bout.csv" , "r")
phys_size = float(fp.readline().split(": ")[1][:-3]) #get physical info on first line
print(phys_size)
mylen = int(fp.readline().split(": ")[1])
print(mylen)
mylen=50
bx_array = np.zeros((mylen,mylen,mylen))
by_array = np.zeros((mylen,mylen,mylen))
bz_array = np.zeros((mylen,mylen,mylen))
params_array = np.zeros(2)
params_array[0]=phys_size
params_array[1]=mylen
fp.readline()
for line in fp.readlines():
    strline = line.split(",")
    i = int(strline[0])
    j = int(strline[1])
    k = int(strline[2])
    bx = float(strline[3])
    by = float(strline[4])
    bz = float(strline[5])
    bx_array[k][j][i] += bx
    by_array[k][j][i] += by
    bz_array[k][j][i] += bz
    
bx_array = np.nan_to_num(bx_array)
by_array = np.nan_to_num(by_array)
bz_array = np.nan_to_num(bz_array)    

#testing to see if it actually works
'''
for i in range(mylen):
    for j in range(mylen):
        for k in range(mylen):
            
            myx = i-25
            myy=j-25
            if myx==0:
                myx=1e-16
            phi = np.arctan(myy/myx)
            
            #bx_array[k][j][i] = np.sin(phi)
            #by_array[k][j][i] = -np.cos(phi)
            #bz_array[k][j][i] = 0
            bx_array[k][j][i] = -1
            by_array[k][j][i] = 1
            bz_array[j][k][i] = .45
'''


def randseed(seednum):
    seedlist = []
    for i in range(seednum):
        tmpseed_arr = np.round(50*np.random.rand(3))
        seedlist.append(tmpseed_arr)
    return seedlist
#seedlist = randseed(seednum)


def orderedseed(seednum): #not currently working, may want to implement this though
    seedlist = []
    myrange = np.arange(0,50,seednum)
    for i in range(seednum):
        seedlist.append([myrange[i],myrange[i],myrange[i]])
    return seedlist        


def trilin_interpolate(myarray, myx, myy, myz):
    
    x0 = int(np.floor(myx))
    x1 = x0+1
    y0 = int(np.floor(myy))
    y1 = y0+1
    z0 = int(np.floor(myz))
    z1 = z0+1
            
    #now define fractional coordinate
    
    #print(myx,x1, x0)
    xd = (myx-x0)/(x1-x0)
    yd = (myy-y0)/(y1-y0)
    zd = (myz-z0)/(z1-z0)
    
    #interpolating along x
            
    c000 = myarray[z0][y0][x0]
    c001 = myarray[z1][y0][x0]
    c010 = myarray[z0][y1][x0]
    c011 = myarray[z1][y1][x0]
    c100 = myarray[z0][y0][x1]
    c101 = myarray[z1][y0][x1]
    c110 = myarray[z0][y1][x1]
    c111 = myarray[z1][y1][x1]
            
            
            
    c00 = c000*(1-xd)+c100*xd
    c01 = c001*(1-xd)+c101*xd
    c10 = c010*(1-xd)+c110*xd
    c11 = c011*(1-xd)+c111*xd
            
    #interpolate long y
    c0 = c00*(1-yd) + c10*yd
    c1 = c01*(1-yd) + c11*yd
    
    #interpolate along z
    c=c0*(1-zd) + c1*zd #this is our value of b at the interpolated point!
    return c

def return_fieldlines(seedlist,filepath):
    tot_xlist = []
    tot_ylist = []
    tot_zlist = []
    for seed in seedlist:
        xlist = []
        ylist = []
        zlist = []
        x=int(seed[0])
        y=int(seed[1])
        z=int(seed[2])
     
        '''think more about stopping conditions, number of steps, too close to BH, wrap back around on itself, etc.'''
        
        tcount =0
        
        boundary_stop = x<mylen-1 and y<mylen-1 and z<mylen-1 and x>=0 and y>=0 and z>=0
        while (boundary_stop and tcount <=100000): 
            xlist.append(x)
            ylist.append(y)
            zlist.append(z)
            
            #let's instead try to interpolate the xyz values from the grid
            #first find the upper and lower values for each coordinate
            #make this a func so we can just pass it to our individual arrays
            bx = trilin_interpolate(bx_array, x,y,z)
            by = trilin_interpolate(by_array,x,y,z)
            bz = trilin_interpolate(bz_array,x,y,z)
            
            
            
            #bx = bx_array[z][y][x]
            #by = by_array[z][y][x]
            #bz = bz_array[z][y][x]
            
            
            
            #normalize b to a meaning value in terms of steps
            
            
            btot = np.sqrt(bx**2 + by**2 + bz**2)
            #print(btot)
            if np.abs(btot)>0:
                bx/=btot
                by/=btot
                bz/=btot #picking out unit vecs
                #print(np.round(bx), np.round(by), np.round(bz))
            
            else: #break if we hit this weird spot
                break
            
            
            #update xyz positions
            
            #x+=int(np.round(bx))
            #y+=int(np.round(by))
            #z+=int(np.round(bz))
            
            x+=bx
            y+=by
            z+=bz
            
            tcount +=1
            
        #update stopping criteria
            boundary_stop = x<mylen-1 and y<mylen-1 and z<mylen-1 and x>=0 and y>=0 and z>=0
        tot_xlist.append(xlist)
        tot_ylist.append(ylist)
        tot_zlist.append(zlist)
        
        
        
        
    return tot_xlist, tot_ylist, tot_zlist


def write_fieldlines(out_directory,xlist,ylist,zlist):
    #xlist,ylist,zlist are a list of lists, the sublists are individual coordinates of a specific fieldline
    numflds = len(xlist)
    #for each fieldline, write individual file
    for i in range(numflds):
        #open file associated with this fieldline
        filename = out_directory+"fieldline"+"%03d"%i+".csv"
        fp = open(filename,'w')
        myx = xlist[i]
        myy = xlist[i]
        myz = zlist[i]
        for j in range(len(myx)):
            fp.write('{},{},{}\n'.format(myx[j],myy[j],myz[j]))
            
        fp.close()
    
    
seedlist = randseed(5)

xlist, ylist, zlist = return_fieldlines(seedlist,'tmp')  

write_fieldlines("blines_csv/fieldlines/test/",xlist,ylist,zlist)




'''
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.view_init(0,0)
c_list = ["Red","Blue","Orange","Green","Magenta","Yellow","Cyan"]
for i in range(len(xlist)):
    ax.plot(xlist[i],ylist[i],zlist[i])
    #plt.savefig('prtl_plots/prtl_{}_test.png'.format(sparse))
#plt.xlim(0,50)
#plt.ylim(0,50)
plt.show()
'''


