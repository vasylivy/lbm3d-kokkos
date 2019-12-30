import numpy as np
import matplotlib.pyplot as plt
import sys
import glob

dir = sys.argv[1];
NX = int(sys.argv[2]);
NY = int(sys.argv[3]);
NZ = int(sys.argv[4]);
umax = float(sys.argv[5]);
cuda = int(sys.argv[6]);
re = float(sys.argv[7]);
fn_data = sys.argv[8];

files = glob.glob1(dir,"*.bin");

frame = 0;
for file in files:
    frame = max(frame,int(file.split(".")[0].split("_")[1]));

data = np.loadtxt(fn_data, delimiter=',');
    
tag = '_' + '%03d' % frame + '.bin';

fn_u = dir + '/' + 'u' + tag;
   
u = np.fromfile(fn_u, dtype="double");
    
if (cuda) : 
    u = np.reshape(u, (NY,NX,NZ), order='F');
else:
    u = np.reshape(u, (NY,NX,NZ));
    
u = u[1:-1,1:-1,1:-1] / umax;
               
ucenter = np.zeros(NY);
ucenter[NY-1] = 1.0; 

if ((NX % 2 == 0) and (NZ % 2 == 0)):
    j = (NX - 2)//2 - 1;
    k = (NZ - 2)//2 - 1;
    ucenter[1:-1] = 0.25*(u[:,j,k] + u[:,j + 1,k] + u[:,j + 1,k + 1] + u[:,j,k + 1]);
else: 
    j = (NX - 1)//2 - 1
    k = (NZ - 1)//2 - 1;    
    ucenter[1:-1] = u[:,j,k];

y = np.zeros(NY);
y[NY-1] = 1.0;
y[1:-1] = np.linspace(0.5, (NY-2) - 0.5, (NY-2)) / (NY - 2);

plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')

fig1 = plt.figure(figsize=(8,8), dpi=100);
plt.plot(ucenter,y,'b-',linewidth=2.5);

if (re == 100) : 
    plt.plot(data[:,1],data[:,0],'ko--')
elif (re == 400) : 
    plt.plot(data[:,2],data[:,0],'ko--')
elif (re == 1000) :
    plt.plot(data[:,3], data[:,0],'ko--')

plt.grid(axis='both',which='both');
plt.ylabel('y-position', fontsize=18);
plt.xlabel('u-component',fontsize=18);
plt.ylim([0, 1]);
plt.legend(['lbm', 'Jiang et al. (1994)'],fontsize=18)
fig1.savefig("re" + str(int(re)) + "_" + str(NY) + "x" + str(NX) + "x" + str(NZ) + ".png");

