#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:00:08 2019

@author: dball
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#first check for existence of saved numpy arrays
#pathext = "arr50_M200"

print('loading from CSV')
fp = open("CSV_out/3D_arr50_M200_nospace.csv" , "r")
phys_size = float(fp.readline().split(": ")[1][:-3]) #get physical info on first line
print(phys_size)
mylen = int(fp.readline().split(": ")[1])
print(mylen)
b_array = np.zeros((mylen,mylen,mylen))
ne_array = np.zeros((mylen,mylen,mylen))
te_array = np.zeros((mylen,mylen,mylen))
params_array = np.zeros(2)
params_array[0]=phys_size
params_array[1]=mylen
print(fp.readline())
for line in fp.readlines():
    strline = line.split(",")
    i = int(strline[0])
    j = int(strline[1])
    k = int(strline[2])
    b = float(strline[3])
    ne = float(strline[4])
    te = float(strline[5])
    b_array[k][j][i] += b 
    #print(b)
   
    ne_array[k][j][i] += ne
    te_array[k][j][i] += te 
    
fp.close
    #clean up the data a bit
ne_array = np.nan_to_num(ne_array)
te_array = np.nan_to_num(te_array)
b_array = np.nan_to_num(b_array)
    
fp = open('cleaned_CSV_flds/3D_arr50_M200.csv','w')
mylen = np.shape(ne_array)[0]
fp.write('i,j,k,b,ne,te\n')
for i in range(mylen):
    for j in range(mylen):
        for k in range(mylen):
            ne=ne_array[k][j][i]
            te=te_array[k][j][i]
            b=b_array[k][j][i]
            fp.write("{},{},{},{},{},{}\n".format(i,j,k,b,ne,te))
fp.close()



'''
tup=100
tlow=.1
#te_array = np.clip(te_array,tlow,tup)
    
denslow=0.0001
densup=5

#ne_array = np.clip(ne_array,denslow,densup)
 
ne_slice = ne_array[25,:,:]
te_slice = te_array[:,25,:]
b_slice = b_array[:,25,:]

print('dens max min')
print(np.max(ne_array),np.min(ne_array))

print('temp max min')
print(np.max(te_slice),np.min(te_slice))

print('b max min')
print(np.max(b_slice), np.min(b_slice))

fig, (ax0, ax1, ax2) = plt.subplots(1,3)
ax0.imshow(np.log10(ne_slice), origin='lower')
ax1.imshow(np.log10(te_slice),origin='lower')
ax2.imshow(np.log10(b_slice),origin='lower')
plt.show()
'''

    

