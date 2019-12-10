import numpy as np
import math

def split_dens(nparr_path,sparse_fac=1.):
    #take density, tmp as an input, returns particles randomly spaced in cell (w/ size deltax, deltay, deltaz)
    ne_path=nparr_path + "ne.npy"
    te_path=nparr_path + "te.npy"
    params_path=nparr_path+"params.npy"
    #load relevant files
    ne = np.load(ne_path)
    te = np.load(te_path)
    te=np.nan_to_num(te)
    print('max, min, avg of temp are : {}, {}, {}'.format(np.max(te), np.min(te), np.mean(te)))
    #there are some pretty ridiculous numerical artifacst yielding very large values of te
    #let's just apply some upper and lower bounds

    
    params = np.load(params_path)
    totdens=np.sum(ne)
    print('sum of density is {}'.format(totdens))
    phys_size, arr_len = params[0], params[1]
    arr_len=int(arr_len)
    xmin = -phys_size/2.
    xmax = phys_size/2. 
    deltax = phys_size / arr_len #width in which to distribute the particles
    #construct coordinate spacing grid
    #coord_arr = np.linspace(xmin,xmax,arr_len)
    #print(coord_arr)

    filepath = "prtl_csv/limlogtest_prtls_sparse{}.csv".format(sparse_fac)
    fp = open(filepath,'w')

    ne =ne*sparse_fac
    print(np.max(ne))
    fp.write("x,y,z,temperature\n")
    totprts = 0
    avgtmp = np.sum(te)
    med_tmp = np.median(te)
    for i in range(arr_len):
        x=i*deltax - phys_size/2.
        for j in range(arr_len):
            y=j*deltax - phys_size/2.
            for k in range(arr_len):
                z=k*deltax - phys_size/2. 
                dens=ne[k][j][i]
                temp = te[k][j][i] #maybe turn this into maxwellian to draw particle energies
                if temp<3:
                    temp=np.log10(10)
                elif temp > 30:
                    temp=np.log10(20)
                else:
                    temp=np.log10(temp)
                numprts=dens
                numprts=math.floor(numprts)
                leftover=dens-numprts
                #now add particle in MC-like way
                myrand = np.random.rand(1)[0]
              
                if myrand < leftover: #add particle if our random number is less than fractional particle density
                    numprts +=1
                #print("numprts : {}".format(numprts))
                for prtnum in range(int(numprts)):
                    myx = x+((np.random.rand(1)[0]-0.5)*deltax) #give random spread to particles w/in cell of width deltax
                    myy = y+((np.random.rand(1)[0]-0.5)*deltax)
                    myz = z+((np.random.rand(1)[0]-0.5)*deltax)
                    totprts+=1
                    fp.write("{},{},{},{}\n".format(myx,myy,myz,temp))
    fp.close()
    print('wrote {} particles'.format(totprts))
    return 0

split_dens("np_arr_data/test_",.4)
                



    



   

tmppath = "np_arr_data/test_"

