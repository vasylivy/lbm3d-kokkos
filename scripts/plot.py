import numpy as np
import matplotlib.pyplot as plt
import sys

dir = sys.argv[1];
frame = int(sys.argv[2]);
NX = int(sys.argv[3]);
NY = int(sys.argv[4]);
NZ = int(sys.argv[5]);
umax = float(sys.argv[6]);
cuda = int(sys.argv[7]);
re = float(sys.argv[8]);
fn_data = sys.argv[9];
fn_out = sys.argv[10];

data = np.loadtxt(fn_data, delimiter=',');
    
tag = '_' + '%03d' % frame + '.bin';

fn_u = dir + '/' + 'u' + tag;
   
u = np.fromfile(fn_u, dtype="double");
    
if (cuda) : 
    u = np.reshape(u, (NY,NX,NZ), order='F');
else:
    u = np.reshape(u, (NY,NX,NZ));
    
u = u[1:-1,1:-1,1:-1];
               
ucenter = np.zeros(NY);
ucenter[NY-1] = 1.0; 
ucenter[1:-1] = 0.25*(u[:,NX//2,NZ//2] + u[:,NX//2 + 1,NZ//2] + u[:,NX//2,NZ//2 + 1] + u[:,NX//2 + 1,NZ//2 + 1]) / umax;

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
fig1.savefig(fn_out + ".png");

