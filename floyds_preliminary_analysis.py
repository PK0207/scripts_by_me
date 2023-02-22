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
red_lamp_files = sorted(glob('ttflat*red*.fits', recursive=True))
blue_lamp_files = sorted(glob('ttflat*blue*.fits', recursive=True))
red_spectrum_files = sorted(glob('nttGD71*red*.fits', recursive=True))
blue_spectrum_files = sorted(glob('nttGD71*red*.fits', recursive=True))

files = red_spectrum_files
images = [fits.open(f)[0].data for f in files]
data_type = 'red spectrum'
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
plt.plot(np.abs(masked_xf), yf)
#plt.xlim(-0.1,0.1)
plt.xlabel('Frequency')
plt.ylabel('Power')
plt.show()
#%% For every image create the power spectrum
width = np.shape(images[0])[1]
fringe_freq = []
for image in images:
    mean_linecut = np.mean(image,axis=0)
    N = 1*width
    yf = fft(mean_linecut)
    xf = fftfreq(N)
    
    #Mask out the region around zero
    masked_freqs = ma.masked_inside(xf,-0.015, 0.015)
    frequency = np.abs(xf[~masked_freqs.mask])
    power = np.abs(yf[~masked_freqs.mask])
    plt.title(f'Power Spectrum of the {data_type}')
    plt.plot(frequency, power)
    plt.xlim(0,0.1)
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.show()
    
    fringe_freq.append(frequency[np.argmax(power)])
    print(f'best frequency {frequency[np.argmax(power)]}')

#%% Telescope parameters
header_values = ['UTSTART', 'EXPTIME', 'AZIMUTH', 'ALTITUDE', 'ROTANGLE', 'APERWID', 'CCDATEMP']
headers = [fits.open(f)[0].header for f in files]
sorting_idx = np.argsort([header['UTSTART'] for header in headers])

env_parameters = []
for header in headers:
    env_parameters.append([header[key] for key in header_values])

env_parameters = np.array(env_parameters)

#plt.style.use('seaborn')
fig, axes = plt.subplots(2, 3, sharex=True, figsize=(30,20))
fig.tight_layout(pad = 10)
x = np.arange(len(env_parameters))
#x = np.sort([header['UTSTART'] for header in headers])
for ax, i in zip(axes.flatten(), range(len(env_parameters[0]))):
    ax.scatter(x, env_parameters[:,i][sorting_idx], s=15)
    ax.set_ylabel(header_values[i], fontsize=20)
    ax.set_xlabel('Image Number', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=14, size=10)

plt.show()
#%% Fringing residuals
def normalize_2d(matrix):
    norm = np.linalg.norm(matrix)
    matrix = matrix/norm  # normalized matrix
    return matrix

template = normalize_2d(images[16])
plt.style.use('default')
for i, im in enumerate(images):
    plt.imshow(np.subtract(template,normalize_2d(im)))
    plt.colorbar()
    plt.title(f'{data_type} {i}, residual')
    plt.show()
#%%
fig = plt.figure(figsize=(10,10), dpi=100)
ax = fig.add_subplot(111, projection='3d')
azimuth = np.array(env_parameters[:,2], dtype='float')
altitude = np.array(env_parameters[:,3], dtype='float')
rotangle = np.array(env_parameters[:,4], dtype='float')
img = ax.scatter(azimuth, altitude, rotangle, c=fringe_freq, vmin=0.0354, vmax=0.0355)
ax.set_xlabel('Azimuth')
ax.set_ylabel('Altitude')
ax.set_zlabel('Rotangle')
cb = fig.colorbar(img)
cb.ax.set_ylabel('fringe frequency', rotation=270)
plt.show()
