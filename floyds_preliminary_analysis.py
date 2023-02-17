# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:36:51 2023

@author: prera
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
#%%
red_lamp_files = glob('ttflat*red*.fits', recursive=True)
blue_lamp_files = glob('ttflat*blue*.fits', recursive=True)
red_spectrum_files = glob('ttGD71*red*.fits', recursive=True)
blue_spectrum_files = glob('ttGD71*red*.fits', recursive=True)

files = blue_lamp_files
images = [fits.open(f)[0].data for f in files]
data_type = 'blue lampflats'
#%%
for i, im in enumerate(images):
    plt.imshow(im)
    plt.colorbar()
    plt.title(f'{data_type}, {i}')
    plt.show()
#%% Mean linecut
mean_linecut = np.mean(images[0],axis=0)
width = np.shape(images[0])[1]
plt.plot(np.linspace(0,width,width), mean_linecut)
plt.title(f'Mean along the length of the {data_type} spectrum')
plt.xlabel('x pixel')
plt.ylabel('Mean pixel value')
#%% Non-mean linecut
linecut = images[0][45]
plt.plot(np.linspace(0,width,width), linecut)
plt.title(f'linecut along the length of the {data_type} spectrum')
plt.xlabel('x pixel')
plt.ylabel('Mean pixel value')
#%% Fourier Transform of the linecut 
from scipy.fft import fft, fftfreq
import numpy.ma as ma
#N = SAMPLE_RATE * DURATION
N = 1*width
yf = fft(mean_linecut)
xf = fftfreq(N)

#Mask out the region around zero
masked_xf = ma.masked_inside(xf, -0.015, 0.015)
plt.title(f'Power Spectrum of the {data_type}')
plt.plot(masked_xf, np.abs(yf))
plt.xlim(-0.1,0.1)
plt.xlabel('Frequency')
plt.ylabel('Power')
plt.show()