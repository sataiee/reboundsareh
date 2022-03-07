# This file read a fargo output and convert it to a two columns style which can be read by rebound
import numpy as np
import os
import matplotlib.pyplot as plt

def GiveBasics(dirc):
  __, __, __, __, __, __, nr, ns = np.loadtxt(dirc+'/dims.dat')
  Nr = int(nr); Ns = int(ns)
  rad = np.loadtxt(dirc+'/used_rad.dat')
  if os.path.isfile(dirc+'/used_azi.dat'):
    thet = np.loadtxt(dirc+'/used_azi.dat')
  else:
    thet = np.arange(ns) * 2*pi/ns
  return rad, thet, Nr, Ns

def Give2DArray(dirc, field, nout, Nr, Ns):
  arr = np.fromfile('{}/gas{}{}.dat'.format(dirc,field, nout))
  arr.resize(Nr,Ns)
  return arr
  
def MakeSigmaFile(dir_source, dir_aim, nout):  
    rad, _, Nr, Ns = GiveBasics(dir_source)
    arr = Give2DArray(dir_source, 'dens', nout, Nr, Ns)
    rmed = 0.5*(rad[:-1]+rad[1:])
    arr_av = np.average(arr, axis=1)
    with open(dir_aim+'/sigma.dat', 'w') as f:
        for a, b in zip(rmed, arr_av):
            f.write('%e\t%e\n'%(a, b))

dir_source = str(input("Please give the source directory (fargo):  "))
dir_aim    = str(input("Please give the n-body directory (rebound):  "))
nout       = int(input("Which fargo output you want to be used?:  "))
MakeSigmaFile(dir_source, dir_aim, nout)
print ("Your sigma.dat file is made in %s"%dir_aim)
