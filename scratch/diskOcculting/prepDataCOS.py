# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 10:22:58 2025

@author: ecruzaguirre
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob #leave this and astropy for if I need it later
from astropy.io import fits

def openCOSRecovered(csvfile):
    data=pd.read_csv(csvfile,encoding='latin-1')
    keys=data.keys()
    wavelength=data[keys[0]].to_numpy()
    recovflux=data[keys[1]].to_numpy()
    recoverr=data[keys[2]].to_numpy()
    return wavelength,recovflux,recoverr 

def openULLYSES(uname):
    #open ULLYSES preview-spec files
    hdu=fits.open(uname)
    dat=hdu[1].data
    wave=dat['WAVELENGTH'][0]
    flux=dat['FLUX'][0]
    errr=dat['ERROR'][0]
    return wave,flux,errr

def plotPower(array):
    pltpow=np.floor(np.log10(np.nanmax(array)))
    pltdiv=10**(pltpow)
    return int(pltpow),pltdiv

def waveVrad(wav,cen):
    c=2.998e5 #speed of light in km/s
    v=c*((wav-cen)/cen)
    return v

def vradWave(rad,cen):
    c=2.998e5 #speed of light in km/s
    w=cen*(1+rad/c)
    return w



#files=glob('*_x1d.fits')
files=['C:/Users/ecruzaguirre/Documents/youngStars/diskOcculting/EP_Cha/hlsp_ullyses_hst_cos-stis_recx-11_uv-opt_dr7_preview-spec.fits']
auxes=['C:/Users/ecruzaguirre/Documents/youngStars/Modeling/RECX 11/RECX11_Recovered_1PA.csv']

namS,wavS,flxS,errS=[],[],[],[]
for file,aux in zip(files,auxes):
    #ww,ff,ee=openCOSRecovered(file) #open COS data that has had airglow subtracted from LyA (OI subtraction not performed)
    ww,ff,ee=openULLYSES(file)
    namS.append(aux.split('/')[-2])
    wavS.append(ww)
    flxS.append(ff)
    errS.append(ee)
    
#plot the data
ppow,pscl=plotPower(flxS[0][8400:])
fig,ax=plt.subplots()
for i in range(len(wavS)):
    plt.plot(wavS[i],flxS[i]/pscl,linewidth=2,label=namS[i])
    plt.fill_between(wavS[i],(flxS[i]+errS[i])/pscl,(flxS[i]-errS[i])/pscl,alpha=0.35)
plt.xlabel('Wavelength ($\\AA$)',fontsize=18)
plt.ylabel('Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',fontsize=18)
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)
plt.show()

#prepare the data for use with Nicole's model
waveRange=[1150,1650] #from PDS 70 proposal
mask=(waveRange[0]<=wavS[0])&(wavS[0]<=waveRange[1])
wavC=wavS[0][mask]
flxC=flxS[0][mask]
errC=errS[0][mask]
namC='EP Cha'

'''
notes from France et al. 2012 (https://ui.adsabs.harvard.edu/abs/2012ApJ...756..171F/abstract)
-----------------------------
Lines in several progressions were identified for EP Cha (RECX 11)
Av is 0.0 (Table 1)
Inclination is 70 deg (Table 1)
Mass is 0.8 Msun (Table 1)
'''

#data from Table 3 of F+12
progres=['[1,4]',
         '[1,7]',
         '[3,13]',
         '[4,13]',
         '[3,16]',
         '[4,4]',
         '[3,0]',
         '[0,1]',
         '[0,2]',
         '[0,3]'] #progression
lyaWave=[1216.07,
         1215.73,
         1213.36,
         1213.68,
         1214.47,
         1214.78,
         1217.04,
         1217.21,
         1217.64,
         1219.09] #LyA pumping wavelength
labWave=[[1431.01,1489.57,1446.12,1504.76],
         [1467.08,1500.45,1524.65,1556.87],
         [1608.33,1615.43],
         [1415.33,1509.45,1613.99],
         [1418.23,1513.99,1593.26,1621.12],
         [1477.05,1526.55,1613.72],
         [1435.05,1591.32,1636.34],
         [1398.95,1460.17,1521.59],
         [1402.65,1463.83,1525.15],
         [1395.20,1407.29,1468.39]] #Ang
vExcite=[[1,      1,      1,      1],      
         [1,      1,      1,      1],
         [3,      3],
         [4,      4,      4],
         [3,      3,      3,      3],
         [4,      4,      4],
         [3,      3,      3],
         [0,      0,      0],
         [0,      0,      0],
         [0,      0,      0]]
vGround=[[6,      7,      6,      7],      
         [6,      7,      7,      8],
         [9,      10],
         [6,      8,      11],
         [5,      7,      9,      10],
         [8,      9,      11],
         [7,      10,     11],
         [5,      6,      2],
         [5,      6,      7],
         [5,      5,      6]]
jExcite=[[4,      4,      4,      4],      
         [7,      7,      7,      7],
         [13,     13],
         [13,     13,     13],
         [16,     16,     16,     16],
         [4,      4,      4],
         [0,      0,      0],
         [1,      1,      1],
         [2,      2,      2],
         [3,      3,      3]]
jGround=[[3,      3,      5,      5],      
         [8,      6,      8,      6],
         [14,     12],
         [12,     12,     12],
         [15,     15,     15,     15],
         [5,      5,      5],
         [1,      1,      1],
         [2,      2,      2],
         [3,      3,      3],
         [2,      4,      4]]

#other relevant data
radVel=18.0 #km/s, from SIMBAD
ismVel=-1.35 #km/s, from LISM Kinematic Calculator
distPc=98.816 #pc, from SIMBAD
extinc=0.0 #Av, extinction coefficient, from K+12
stMass=0.8 #M_sun, from K+12
dskInc=70 #degrees, from K+12

#plot COS data and emission lines
ppow,pscl=plotPower(flxC[7500:])
fig,ax=plt.subplots(figsize=(16,8))
plt.plot(wavC,flxC/pscl,linewidth=2,label=namC,color='m')
plt.fill_between(wavC,(flxC+errC)/pscl,(flxC-errC)/pscl,alpha=0.35,color='m')
for p in range(len(progres)):
    plt.plot([],[],color=f'C{p}',linewidth=1,label=progres[p])
    for w in labWave[p]: 
        plt.axvline(w,color=f'C{p}',linewidth=1)
    plt.axvline(lyaWave[p],color=f'C{p}',linewidth=1,linestyle='--')
plt.plot([],[],color='k',linewidth=1,linestyle='--',label='Ly$\\alpha$ Pump')
plt.xlabel('Wavelength ($\\AA$)',fontsize=18)
plt.ylabel('Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',fontsize=18)
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)
plt.xlim(waveRange)
plt.ylim(-0.03,1.25)
plt.show()

goodProg=[True,True,False,False,False,
          False,False,True,False,False] #if there is at least one good line,mask with true
goodLine=[[True,True,True,True],[True,True,True,True],[],[],[],
          [],[],[True,True,False],[],[]] #for good progressions, explicitly call out which lines are to be extracted

#isolate individual emission lines that are usable
#wavelength arrays are typically +/-150 km/s, or if off-center, cover a ~300 km/s range (for COS)
#the G140L data only covers ~3 points in 300 km/s, need to extend to capture the full line
width=300 #km/s width
labS,prgS,lyaS=[],[],[]
velI,flxI,errI=[],[],[]
for p in range(len(progres)):
    if goodProg[p]:
        for w in range(len(labWave[p])):
            if goodLine[p][w]:
                velC=waveVrad(wavC,labWave[p][w]) #convert to km/s centered on lab wavelength
                velMask=((-width/2)<=velC)&(velC<=(width/2))
                veli=velC[velMask]
                flxi=flxC[velMask]
                erri=errC[velMask]
                labS.append(labWave[p][w])
                prgS.append(progres[p])
                lyaS.append(lyaWave[p])
                velI.append(veli)
                flxI.append(flxi)
                errI.append(erri)
    
#inspect the extracted emission lines, make sure width is good
fig,ax=plt.subplots()
for i in range(len(velI)):
    plt.plot(velI[i],flxI[i],label=f'{prgS[i]}: {labS[i]}')
    plt.fill_between(velI[i],flxI[i]+errI[i],flxI[i]-errI[i],alpha=0.15)
plt.legend()
plt.show()

#store everything into a dictionary
saveData={}
for i in range(len(labS)):
    key=str(labS[i])
    saveData[key]=[labS[i],radVel,prgS[i],lyaS[i],velI[i],flxI[i],errI[i]]
    
#save the data to be used for modeling
#np.savez('EP_Cha_COS_ULLYSES_H2_Data.npz',data=saveData)
#save laboratory centroid (A), velocity shift (km/s), progression, velocity array (km/s), flux array (erg/cm^2/s/A), error array (erg/cm^2/s/A)
#need to save as dictionary to deal with ragged arrays