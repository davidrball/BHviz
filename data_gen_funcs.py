import numpy as np
import time
def get_vals(snapshot,write_bool=0): #generate arrays of values
    fp = open(snapshot , "r")
    phys_size = float(fp.readline().split(": ")[1][:-3]) #get physical info on first line
    print(phys_size)
    mylen = int(fp.readline().split(": ")[1])
    print(mylen)
    b_array = np.zeros((mylen,mylen,mylen))
    ne_array = np.zeros((mylen,mylen,mylen))
    te_array = np.zeros((mylen,mylen,mylen))

    param_array = np.array([phys_size,mylen])

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
    if write_bool==1: #should always save, way faster to read arrays
        np.save("np_arr_data/test_params.npy", param_array)
        np.save("np_arr_data/test_ne.npy",ne_array)
        np.save("np_arr_data/test_te.npy",te_array)
        np.save("np_arr_data/test_b.npy",b_array)
    return ne_array, te_array, b_array, param_array
#timing get_vals func to see if it's faster than loading numpy arrays

#start_time = time.time()
get_vals("CSV_out/3D_test_out.txt",1)
#print("{} seconds to create arrays".format(time.time()-start_time))

#testing loading capabilities
#start_time = time.time()
ne=np.load("np_arr_data/test_ne.npy")
te=np.load("np_arr_data/test_te.npy")
b=np.load("np_arr_data/test_b.npy")
params=np.load("np_arr_data/test_params.npy")
print(params)
#print("{} seconds to read arrays".format(time.time() - start_time))

def write_vals(valdict):

    return 0
    #write csv file in format required by unity