import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
        "-reynolds",
        nargs="*",
        type=int,
    )

parser.add_argument(
        "-nx",
        type=int,
    )
 
parser.add_argument(
        "-ny",
        type=int,
    )

parser.add_argument(
        "-nz",
        type=int,
    )
 
parser.add_argument(
        "-umax",
        type=float,
    )
 
parser.add_argument(
        "-cuda",
        type=float,
        default=0,
    )
 
parser.add_argument(
        "-data",
        type=str,
        default = "jiang_data.csv",
    )

args = parser.parse_args();
 
reynolds = args.reynolds;
NX = args.nx;
NY = args.ny;
NZ = args.nz;
umax = args.umax;
cuda = args.cuda;
fn_data = args.data;

plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
fig1 = plt.figure(1,figsize=(8,8), dpi=100);
plt.grid(axis='both',which='both');
plt.ylabel('y-position', fontsize=18);
plt.xlabel('u-component',fontsize=18);
plt.ylim([0, 1]);

for re in reynolds:
    
    dir = "re" + str(re) + "_" + str(NY) + "x" + str(NY) + "x" + str(NZ) + "/output"

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
    
    plt.figure(1)
    plt.plot(ucenter,y,'b-',linewidth=2.5);
    
    if (re == 100) : 
        plt.plot(data[:,1],data[:,0],'ko--')
    elif (re == 400) : 
        plt.plot(data[:,2],data[:,0],'ko--')
    elif (re == 1000) :
        plt.plot(data[:,3], data[:,0],'ko--')

plt.figure(1)
plt.legend(['lbm', 'Jiang et al. (1994)'],fontsize=18)
fig1.savefig(str(NY) + "x" + str(NX) + "x" + str(NZ) + ".png");

