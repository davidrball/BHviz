import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

sparse=0.3
myf = open("prtl_csv/test_prtls_sparse{}.csv".format(sparse),'r')

x_list = []
y_list = []
z_list = []
t_list = []
c_list = []
myf.readline()

for line in myf.readlines():
    tmplist = line.split(',')
    #print(tmplist)
    x,y,z,tmp = tmplist[0],tmplist[1],tmplist[2],tmplist[3]
    x_list.append(float(x))
    y_list.append(float(y))
    z_list.append(float(z))
    t_list.append(float(tmp))

t_list = np.array(t_list)
tmin = np.min(t_list)
tmax = np.max(t_list)
t_list = (t_list-tmin)/(tmax-tmin)
print(max(t_list),min(t_list))

for t in t_list:
    c_list.append((t,0,1-t))


fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.view_init(0,0)
ax.scatter(x_list,y_list,z_list,c=c_list,marker='o',s=.05)
plt.savefig('prtl_plots/prtl_{}_test.png'.format(sparse))