import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt


#need functions to read in csv, clean it,


#testing tofile
#Nbins=100
#testcube = np.zeros((Nbins,Nbins,Nbins))
#testcube+=0.5

testfilename = "3D_arr50_M50_nohead_nonan.csv"

#output 0 corresponds to dens, 1 to te, and 2 to b
def CSV_to_datacube(filename,densclip_low,densclip_high,output):
    print('loading from CSV '+filename)
    fp = open("CSV_out/"+filename , "r")
    #phys_size = float(fp.readline().split(": ")[1][:-3]) #get physical info on first line
    #print(phys_size)
    #mylen = int(fp.readline().split(": ")[1])
    #print(mylen)
    mylen=50
    phys_size=50
    b_array = np.zeros((mylen,mylen,mylen))
    ne_array = np.zeros((mylen,mylen,mylen))
    te_array = np.zeros((mylen,mylen,mylen))
    params_array = np.zeros(2)
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
    #data should be cleaned already but just in case
    ne_array = np.nan_to_num(ne_array)
    te_array = np.nan_to_num(te_array)
    b_array = np.nan_to_num(b_array)

    #now take log
    
    ne_array = np.nan_to_num(np.log10(ne_array))
    te_array = np.nan_to_num(np.log10(te_array))
    b_array = np.nan_to_num(np.log10(b_array))

    #artifically clip array to reasonable values
    
    if output==0:
        input_arr = ne_array
    elif output==1:
        input_arr = te_array
    elif output==2:
        input_arr = b_array


    #if both are 0, don't clip
    if densclip_low==0 and densclip_high==0:
        out_array = input_arr
    else:
        out_array = np.clip(input_arr,densclip_low,densclip_high)
    
    print(np.min(out_array))
    print(np.max(out_array))

    mymin = np.min(out_array)
    mymax = np.max(out_array)    
    
    
    scaled_array = (ne_array-mymin)/(mymax-mymin)

    return scaled_array


denslow,denshigh = -1.5,1
datacube_dens = CSV_to_datacube(testfilename,denslow,denshigh,0)

templow,temphigh=0,2
datacube_temp = CSV_to_datacube(testfilename,templow,temphigh,1)

blow=-3
bhigh=0
datacube_b=CSV_to_datacube(testfilename,blow,bhigh,2)




#give it a data cube normalized to 0-1 (cube_norm), and it puts out binary file
def tofile(cube_norm,name='binary_out_test',format=8,T=(3,2,1)): #working
    if format==8:
        cube_fmt = (cube_norm*255).astype('uint8')
    elif format==32:
        cube_fmt = cube_norm.astype('float32')
    else:
        print("unsupported format")
        return
    filename="binary_outputs/{}".format(name)
    print("writing " + filename)
    cube_fmt.transpose(T[0]-1,T[1]-1,T[2]-1).tofile(filename)


tofile(datacube_dens,name="dens_binary")
tofile(datacube_temp,name="temp_binary")
tofile(datacube_temp,name="b_binary")
