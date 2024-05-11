import os
import h5py
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser()
parser.add_argument('--x_min', required=True, type=int, help='minimum crop x')
parser.add_argument('--x_max', required=True, type=int, help='maximum crop x')
parser.add_argument('--y_min', required=True, type=int, help='minimum crop y')
parser.add_argument('--y_max', required=True, type=int, help='maximum crop y')
parser.add_argument('--z_min', required=True, type=int, help='minimum crop z')
parser.add_argument('--z_max', required=True, type=int, help='maximum crop z')
parser.add_argument('--res', required=True, type=float, help='maximum resolution')
parser.add_argument('--n_x', required=True, type=int, help='number of cells in x direction')
parser.add_argument('--n_y', required=True, type=int, help='number of cells in y direction')
parser.add_argument('--n_z', required=True, type=int, help='number of cells in z direciton')
parser.add_argument('--n_level', required=True, type=int, help='level of refinement')



opt = parser.parse_args()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

x_min = opt.x_min
x_max = opt.x_max
y_min = opt.y_min
y_max = opt.y_max
z_min = opt.z_min
z_max = opt.z_max

dimX = x_max - x_min  
dimY = y_max - y_min
dimZ = z_max - z_min

res=opt.res

n_x = opt.n_x
n_y = opt.n_y
n_z = opt.n_z

n_level = opt.n_level

ndir=0
for dir in os.listdir(os.getcwd()):
  if (is_number(dir)):
    ndir +=1

time    = np.zeros(ndir);
Species = np.zeros((n_x, n_y, n_z, ndir))


x = np.zeros(dimX*dimY*dimZ)
y = np.zeros(dimX*dimY*dimZ)
z = np.zeros(dimX*dimY*dimZ)

file = open("constant/polyMesh/cellCenters","r")
Lines = file.readlines()
count =0
wbool=0


for line in Lines:
    ls = line.strip()
    if (ls==")"):
      break
    if (wbool==1):
      x[count]=float(ls.split("(")[1].split(")")[0].split()[0])
      y[count]=float(ls.split("(")[1].split(")")[0].split()[1])
      z[count]=float(ls.split("(")[1].split(")")[0].split()[2])
      count +=1
    if (ls=="("):
      wbool=1


resX = res*dimX/n_x
resY = res*dimY/n_y
resZ = res*dimZ/n_z


tcount =0
for dir in os.listdir(os.getcwd()): 
  if (is_number(dir)):
    print(dir)
    time[tcount]=float(dir)

    file = open(dir+"/Species","r")
    Lines = file.readlines()
    count =0
    wbool=0
    for line in Lines:
      ls = line.strip()
      if (ls==")"):
        break
      if (wbool==1):
        a = np.floor((x[count])/resX)
        b = np.floor((y[count])/resY)
        c = np.floor((z[count])/resZ)
        Species[a.astype(int), b.astype(int),c.astype(int),tcount]=float(ls) 
        count +=1
      if (ls=="("):
        wbool=1
    tcount +=1

p=1
for i in range(0,n_level):
    resX=resX/2
    resY=resY/2
    resZ=resZ/2
    p0=p
    p=2*p0
    Species0=Species
    Species=np.zeros((n_x*p, n_y*p, n_z*p,ndir));
    for t in range(0,ndir):
        for i in range(0,n_x*p0):
            for ii in range(0,2):
                for j in range(0,n_y*p0):
                    for jj in range(0,2):
                        for k in range(0,n_z*p0):
                            for kk in range(0,2):
                                Species[2*i+ii,2*j+jj,2*k+kk,t]=Species0[i,j,k,t]

tcount =0
for dir in os.listdir(os.getcwd()):
  if (is_number(dir)):
    print(dir)


    file = open(dir+"/Species","r")
    Lines = file.readlines()
    count =0
    wbool=0
    for line in Lines:
      ls = line.strip()
      if (ls==")"):
        break
      if (wbool==1):
        a = np.floor((x[count])/resX)
        b = np.floor((y[count])/resY)
        c = np.floor((z[count])/resZ)
        Species[a.astype(int), b.astype(int),c.astype(int),tcount] = float(ls) 
        count +=1
      if (ls=="("):
        wbool=1
    tcount +=1


f = h5py.File("voxelisedResults.hdf5","r+")

for name in list(f.keys()):
    if (name=='Species'):
        del f['Species']

Species=np.swapaxes(Species,0,2)

v = np.argsort(time)
time = time[v]
Species=Species[:,:,:,v]

f.create_dataset('Species', data=Species, dtype="float", compression="gzip")




f.close()

