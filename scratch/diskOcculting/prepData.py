# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 10:32:10 2025

@author: ecruzaguirre
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from astropy.io import fits

def openSTIS(fitsfile,order='all',scale=False):
    hdu=fits.open(fitsfile)
    data=hdu[1].data
    
    if order=='all':
        wave=data['wavelength'].flatten()
        flux=data['flux'].flatten()
        errr=data['error'].flatten()  
    
        sort_key=np.argsort(wave)
        wave=wave[sort_key]
        flux=flux[sort_key]
        errr=errr[sort_key] 
        
    else: #assume a valid order is picked
        wave=data['wavelength'][order]
        flux=data['flux'][order]
        errr=data['error'][order]
        
        if not isinstance(scale,bool):
            flux*=scale
            errr*=scale

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



files=glob('*_x1d.fits')
#oeby01010 is G230L
#oeby01020 is G140L

grats=['G230L','G140L']
namS,wavS,flxS,errS=[],[],[],[]
for file,grat in zip(files,grats):
    ww,ff,ee=openSTIS(file)
    namS.append(file.split('_')[0]+' ('+grat+')')
    wavS.append(ww)
    flxS.append(ff)
    errS.append(ee)
    
#plot both raw files
ppow,pscl=plotPower(flxS[1])
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
mask=(waveRange[0]<=wavS[1])&(wavS[1]<=waveRange[1])
wavG=wavS[1][mask]
flxG=flxS[1][mask]
errG=errS[1][mask]
namG='PDS 70 (G140L)'

'''
notes from Skinner and Audard 2022 (https://ui.adsabs.harvard.edu/abs/2022ApJ...938..134S/abstract)
----------------------------------
PDS 70 is a WTTS
H2 lines in WTTSs are rare
5 weak H2 lines were identified, 4 from the [1,4] progression and one from the [1,7] progression
these are pumped by 1216.07 and 1215.73, respectively
there are other weaker H2 lines, but they are not well resolved from stronger atomic emission lines
Extinction value is Av=0.05, the extinction law from Whittet+2004 is used to deredden
'''

#data from Table 4 of SA22
labWave=[1431.01,1489.57,1446.12,1504.76,1467.08] #Ang
obsWave=[1431.1, 1489.6, 1445.8, 1504.6, 1467.3] #Ang
progres=['[1,4]','[1,4]','[1,4]','[1,4]','[1,7]'] #progression
lyaWave=[1216.07,1216.07,1216.07,1216.07,1215.73] #LyA pumping wavelength
vExcite=[1,      1,      1,      1,      1]
vGround=[6,      7,      6,      7,      6]
jExcite=[4,      4,      4,      4,      7]
jGround=[3,      3,      5,      5,      8]

#other relevant data
radVel=3.13 #km/s, from SIMBAD
ismVel=-22.69 #km/s, from LISM Kinematic Calculator
ismErr=1.20 #km/s, LoS is only through the NGP Cloud 
distPc=112 #pc, from proposal
extinc=0.05 #Av, extinction coefficient, from SA22/proposal
stMass=0.76 #M_sun, from proposal
dskInc=52 #degrees, from proposal

#plot only G140L
ppow,pscl=plotPower(flxG)
fig,ax=plt.subplots()
plt.plot(wavG,flxG/pscl,linewidth=2,label=namG)
plt.fill_between(wavG,(flxG+errG)/pscl,(flxG-errG)/pscl,alpha=0.35)
for w in obsWave:
    plt.axvline(w,color='k',linewidth=2)
plt.axvline(w,color='k',linewidth=2,label='H$_2$ Line')
for a in lyaWave:
    plt.axvline(a,color='r',linewidth=2)
plt.axvline(a,color='r',linewidth=2,label='Ly$\\alpha$ Pumping')
plt.xlabel('Wavelength ($\\AA$)',fontsize=18)
plt.ylabel('Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',fontsize=18)
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)
plt.show()

#isolate individual emission lines
#wavelength arrays are typically +/-150 km/s, or if off-center, cover a ~300 km/s range (for COS)
#the G140L data only covers ~3 points in 300 km/s, need to extend to capture the full line
width=2000 #km/s width
velI,flxI,errI,offV=[],[],[],[]
for w in range(len(labWave)):
    velG=waveVrad(wavG,labWave[w]) #convert to km/s centered on lab wavelength
    offsetV=np.sum(waveVrad(np.array([labWave[w],obsWave[w]]),labWave[w]))
    velMask=((offsetV-width/2)<=velG)&(velG<=(offsetV+width/2))
    veli=velG[velMask]
    flxi=flxG[velMask]
    erri=errG[velMask]
    velI.append(veli)
    flxI.append(flxi)
    errI.append(erri)
    offV.append(offsetV)
    
fig,ax=plt.subplots()
for i in range(len(velI)):
    plt.plot(velI[i],flxI[i],label=labWave[w])
    plt.fill_between(velI[i],flxI[i]+errI[i],flxI[i]-errI[i],alpha=0.15)
plt.legend()
plt.show()

#store everything into a dictionary
saveData={}
for i in range(len(labWave)):
    key=str(labWave[i])
    saveData[key]=[labWave[i],obsWave[i],offV[i],progres[i],lyaWave[i],velI[i],flxI[i],errI[i]]
    
#save the data to be used for modeling
#np.savez('PDS70_STIS_G140L_MAST_Data.npz',data=saveData)
#save laboratory centroid (A), obsreved centroid (A), velocity shift (km/s), velocity array (km/s), flux array (erg/cm^2/s/A), error array (erg/cm^2/s/A)
#need to save as dictionary to deal with ragged arrays