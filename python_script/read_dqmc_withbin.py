import sys
import numpy as np
#from pathlib import Path
from numpy import *
from scipy import integrate
import math
from matplotlib.pyplot import *


filenames = sys.argv[1:]
if len(filenames)==0:
   print("No input file.")
   sys.exit(0)

    
def readintpara(file):
    temp = ''
    for s in range (0,len(file)):
        if file[s].isdigit():
            temp += file[s]
    return int(temp)

def readfloatpara(file):
    temp = ''
    for s in range (0,len(file)):
        if (file[s].isdigit() or file[s]=='.'):
            temp += file[s]
    return float(temp)
        
    #[float(s) for s in file.find()]
    #L = int(file[52:56])
    #for s in file.split():
    #    if s.isdigit():
    #        read = int(s)
    #        print('read=',read)
    #return read
    
def geterr(dat, sgn, N, mean, std):
    dat_blocked = dat[:].reshape(10,int(N/10))
    sgn_blocked = sgn[:].reshape(10,int(N/10))
    dat_avgs = np.sum(dat_blocked,axis=0)
    sgn_avgs = np.sum(sgn_blocked,axis=0)
    sum_dat = np.sum(dat_avgs)
    sum_sgn = np.sum(sgn_avgs)
    x_j = (sum_dat - dat_avgs)/(sum_sgn - sgn_avgs)
    mean = sum(x_j)/float(N/10)
    std = np.sqrt(np.dot((x_j-mean),(x_j-mean)) * float((N/10)-1)/float(N/10))

def printdata(filename,ii,G4mean):
    mylines = []                                    # Declare an empty list.    
    with open(filename, 'rt') as myfile:            # Open the file for reading text ('rt')
        for eachline in myfile:                     # For each line in the file,
            mylines.append(eachline.rstrip('\n'))   # strip newline and add to list.
    
    for i in range(4,len(mylines),5):
        if((i+1)%5 == 0):
            ti = (i+1)/5 - 1
            k1 = int(ti/Nc/band/band/band/band/nbin)
            i1 = ti - k1*Nc*band*band*band*band*nbin
            k2 = int(i1/band/band/band/band/nbin)
            i2 = i1 - k2*band*band*band*band*nbin
            a1 = int(i2/band/band/band/nbin)
            i3 = i2 - a1*band*band*band*nbin
            a2 = int(i3/band/band/nbin)
            i4 = i3 - a2*band*band*nbin
            a3 = int(i4/band/nbin)
            i5 = i4 - a3*band*nbin
            a4 = int(i5/nbin)
            nb = int(i5 - a4*nbin)
            #print("k1=",k1,"k2=",k2,"a1=",a1,"a2=",a2,"a3=",a3,"a4=",a4,"nb=",nb)
            #print(mylines[i][3:12], mylines[i][19:28])
            G4mean[k1,a1,a3,k2,a2,a4,nb] = float(mylines[i][3:12])
            sgnt[nb] = float(mylines[i][19:28])
    #for i in range(0,Nx*Ny):
    #    geterr(chippt0dat[i,:], chippt0sgn[i,:], 1080, chippt0mean[i], chippt0std[i])
    
    print("Done reading data")
    
    
    # old integration that includes endpoint
    #Xsc_old = mylines[len(mylines)-1]

    #print('\n')
    #print(Nx)
    #print(mu)
    #print(beta_str)
    #print(navg)
    #print(Xavg)
    #print('\n')
    #print(Xsc_integrand.shape)
    #print(dXsc_integrand.shape)
    #print('Copy line for Xsc and dXsc:')
    #print('\n')
    #print("%.6f" % Xsc, "%.6f" % dXsc)
    #print('\n')

    #print('Old values of Xsc +- dXsc:', Xsc_old[41:68])
    #print('New values of Xsc +- dXsc:   ', "%.6f" % Xsc,'+-    ',"%.6f" % dXsc)
mu = -2.18
Nx = 4
Ny = 4
Nc = Nx*Ny
t = 1
tp = -0.25
U = 6
beta = 8
band = 3
len1 = len(filenames)
nbin = 1080
print('Nx = ',Nx)
print('Ny = ',Ny)
print('t = ',t)
print('tp = ',tp)
print('U = ',U)
print('beta = ',beta)
print('mu = ',mu)
    
G4rmean  = np.zeros((Nx*Ny,band,band,Nx*Ny,band,band,nbin),dtype='double')
G4kmean  = np.zeros((Nx*Ny,band,band,Nx*Ny,band,band,nbin),dtype='complex')
sgnt  = np.zeros(nbin,dtype='double')


print('|------------------------------------------------------------------------------|')
print('|------------------------- Read DQMC data -------------------------------------|')
print('|------------------------------------------------------------------------------|')
print('\n')
print('Filenames: ', filenames)
print('\n')

ii = 0
for fnum in filenames:
    printdata(fnum,ii,G4rmean)
    ii = ii + 1



Rvec  = np.zeros((Nc,2),dtype='int')
for i in range(-1,3):
    for j in range(-1,3):
        Rvec[(i+2-1)*4+j+2-1,0] = i
        Rvec[(i+2-1)*4+j+2-1,1] = j

Kvec = Rvec.copy()
Kvec = Kvec * np.pi/2

bandphase = np.zeros((band,band,2),dtype='double')
bandphase[0,1,0]=1/2
bandphase[1,0,0]=-1/2
bandphase[0,2,1]=1/2
bandphase[2,0,1]=-1/2
bandphase[1,2,0]=-1/2
bandphase[1,2,1]=1/2
bandphase[2,1,0]=1/2
bandphase[2,1,1]=-1/2

for i1 in range (Nc):
    for a1 in range (band):
        for a2 in range (band):
            for i2 in range (Nc):
                for a3 in range (band):
                    for a4 in range (band):
                        for k1 in range (Nc):
                            for k2 in range (Nc):
                                kx1 = Kvec[k1,0]; ky1 = Kvec[k1,1];
                                kx2 = Kvec[k2,0]; ky2 = Kvec[k2,1];
                                rx1 = Rvec[i1,0]; ry1 = Rvec[i1,1];
                                rx2 = Rvec[i2,0]; ry2 = Rvec[i2,1];
                                G4kmean[k1,a1,a2,k2,a3,a4,:] += G4rmean[i1,a1,a2,i2,a3,a4,:] * exp(-1.0 * 1j * (kx1 * rx1 + ky1 * ry1 + kx1 * bandphase[a1,a2,0] + ky1 * bandphase[a1,a2,1])) * exp(1.0 * 1j * (kx2 * rx2 + ky2 * ry2 + kx2 * bandphase[a3,a4,0] + ky2 * bandphase[a3,a4,1]))

print("Done Fourier Transform")
#PS=0.0; PScc=0.0; PSoxox=0.0; PSoyoy=0.0; PScox=0.0; PScoy=0.0; PSoxoy=0.0;
#for ik1 in range(Nc):
#    for ik2 in range(Nc):
#        PScc += G4kmean[ik1,0,0,ik2,0,0]
#        PSoxox += G4kmean[ik1,1,1,ik2,1,1]
#        PSoyoy += G4kmean[ik1,2,2,ik2,2,2]
#        PScox += G4kmean[ik1,0,0,ik2,1,1] + G4kmean[ik1,1,1,ik2,0,0]
#        PScoy += G4kmean[ik1,0,0,ik2,2,2] + G4kmean[ik1,2,2,ik2,0,0]
#        PSoxoy += G4kmean[ik1,1,1,ik2,2,2] + G4kmean[ik1,2,2,ik2,1,1]
        
#PScc = PScc/(float(Nc))
#print("G4 s-wave Pairfield susceptibility for Cu-Cu is ",PScc)
#PSoxox = PSoxox/float(Nc)
#print("G4 s-wave Pairfield susceptibility for Ox-Ox is ",PSoxox)
#PSoyoy = PSoyoy/float(Nc)
#print("G4 s-wave Pairfield susceptibility for Oy-Oy is ",PSoyoy)
#PScox = PScox/float(Nc)
#print("G4 s-wave Pairfield susceptibility for Cu-Ox is ",PScox)
#PScoy = PScoy/float(Nc)
#print("G4 s-wave Pairfield susceptibility for Cu-Oy is ",PScoy)
#PSoxoy = PSoxoy/float(Nc)
#print("G4 s-wave Pairfield susceptibility for Ox-Oy is ",PSoxoy)
#PS = PScc+PSoxox+PSoyoy+PScox+PScoy+PSoxoy
#print("G4 s-wave Pairfield susceptibility is ",PS)
            

rightmatrix = zeros((Nc,band,band),dtype='complex')
leftmatrix = zeros((Nc,band,band),dtype='complex')

leftmatrix[:,0,0]=rightmatrix[:,0,0]=1.0

for k in range(Nc):
    gammak = sqrt(sin(Kvec[k,0]/2)**2+sin(Kvec[k,1]/2)**2)
    if (abs(gammak) < 0.0000000001):
        rightmatrix[k,1,1]= 1j/sqrt(2)
        rightmatrix[k,1,2]= -1j/sqrt(2)
        rightmatrix[k,2,1]= -1j/sqrt(2)
        rightmatrix[k,2,2]= -1j/sqrt(2)
        leftmatrix[k,1,1]= -1j/sqrt(2)
        leftmatrix[k,1,2]= 1j/sqrt(2)
        leftmatrix[k,2,1]= 1j/sqrt(2)
        leftmatrix[k,2,2]= 1j/sqrt(2)
    else:
        kx = Kvec[k,0]; ky = Kvec[k,1];
        rightmatrix[k,1,1]= 1j * sin(kx/2)/gammak
        rightmatrix[k,1,2]= -1j * sin(ky/2)/gammak
        rightmatrix[k,2,1]= -1j * sin(ky/2)/gammak
        rightmatrix[k,2,2]= -1j * sin(kx/2)/gammak
        leftmatrix[k,1,1]= -1j * sin(kx/2)/gammak
        leftmatrix[k,1,2]= 1j * sin(ky/2)/gammak
        leftmatrix[k,2,1]= 1j * sin(ky/2)/gammak
        leftmatrix[k,2,2]= 1j * sin(kx/2)/gammak
G4kmeanLLbar = zeros((Nc,band,band,Nc,band,band,nbin),dtype='complex')
def changeG4LLbar(Tchi0LLbar,Tchi0):
    Tchi0LLbartemp1 = zeros((Nc,band,band,Nc,band,band),dtype='complex')
    for l2 in range(band):
        for l1 in range(band):
            for k1 in range(Nc):
                for k2 in range(Nc):
                    Tchi0LLbartemp1[k1,:,:,k2,l1,l2]=np.dot(np.dot(leftmatrix[k1,:,:],Tchi0[k1,:,:,k2,l1,l2]),rightmatrix[k1,:,:])

    for l4 in range(band):
        for l3 in range(band):
            for k1 in range(Nc):
                for k2 in range(Nc):
                    Tchi0LLbar[k1,l3,l4,k2,:,:]=np.dot(np.dot(rightmatrix[k2,:,:],Tchi0LLbartemp1[k1,l3,l4,k2,:,:]),leftmatrix[k2,:,:])
        
for nb in range(nbin):
    changeG4LLbar(G4kmeanLLbar[:,:,:,:,:,:,nb],G4kmean[:,:,:,:,:,:,nb])


gkd = cos(Kvec[:,0]) - cos(Kvec[:,1])
Pd = zeros((band,band,band,band,nbin),dtype='complex')
for iK1 in range(Nc):
    for iK2 in range(Nc):
        Pd[:,:,:,:,:] += gkd[iK1]*G4kmeanLLbar[iK1,:,:,iK2,:,:,:]*gkd[iK2]
        
Pd /= Nc*beta
Pdall = sum(sum(sum(sum(Pd,axis=0),axis=0),axis=0),axis=0)
Pddddd = Pd[0,0,0,0,:]

print("Pdall = ",mean," +- ",std)
geterr(Pdall, sgnt, nbin, mean, std)
print("Pdall = ",mean," +- ",std)
geterr(Pddddd, sgnt, nbin, mean, std)
print("Pddddd = ",mean," +- ",std)
