# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 10:42:59 2023

@author: prera
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from scipy.fft import fft, fftfreq
import numpy.ma as ma
#%%
files = sorted(glob('lco_data-20230216-33/lco_data-20230216-33/flat*red_*.fits', recursive=True))
#files = sorted(glob('New_AltAz_data/ogg2m001-*w00.fits.fz', recursive=True))
images = [fits.open(f)[0].data for f in files]
data_type = 'red flat'
#%% Select a subregion of this non-rectified flat
for i, im in enumerate(images):
    plt.imshow(im[175:225, 1200:1500], origin='lower')
    plt.colorbar()
    plt.title(f'{data_type}, {i} fringe region')
    plt.show()
#Now look at the linecut
fringe_rms = []
for i, im in enumerate(images):
    plt.plot(np.arange(300), im[200, 1200:1500])
    plt.title(f'{data_type}, {i} fringe region linecut')
    plt.show()
    fringe_rms.append(np.std(im[200, 1200:1500]))
#%% Plot std with image number
plt.figure(dpi=200)
plt.scatter(np.arange(len(images)), fringe_rms, s=5)
plt.title('Fringe STD for each image')
plt.xlabel('Image Number')
plt.ylabel('Standard Deviation of the sub-region')
plt.show()
#%% Plot power spectrum
width = np.shape(images[0])[1]
fringe_freq = []
fringe_power = []
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
    fringe_power.append(np.max(power))
    #print(f'best frequency {frequency[np.argmax(power)]}')
    #print(f'Power at that frequency {np.max(power)}')
#%% Plot power with image number
plt.figure(dpi=200)
plt.scatter(np.arange(len(images)), fringe_power, s=5, label=f'fringe frequency = {fringe_freq[0]}')
plt.title('Power of fringing from fft for each image')
plt.xlabel('Image Number')
plt.ylabel('Power of highest peak in fft')
plt.legend()
plt.show()
#%% Now shift the spectrum a little and remeasure
def xy_shift(data, x_shift=0, y_shift=0, x_up=True, y_up=True):
    if x_up==True and y_up==True:
        data = np.roll(data, x_shift, axis=1)
        data = np.roll(data, y_shift, axis=0)
    if x_up==True and y_up==False:
        data = np.roll(data, x_shift, axis=1)
        data = np.roll(data, y_shift*(-1), axis=0)
    if x_up==False and y_up==True:
        data = np.roll(data, x_shift*(-1), axis=1)
        data = np.roll(data, y_shift, axis=0)
    if x_up==False and y_up==False:
        data = np.roll(data, x_shift*(-1), axis=1)
        data = np.roll(data, y_shift*(-1), axis=0)
    return data

from PIL import Image
def rotate(data, angle=0, clockwise=False):
    pil_image = Image.fromarray(data)
    if clockwise==False:
        new_im = np.array(pil_image.rotate(360-angle))
    else:
        new_im = np.array(pil_image.rotate(angle))
    return new_im

im1_num = 0
im2_num = 1
div_im = np.divide(images[im1_num],images[im2_num])

div_im_std = np.std(div_im[200, 1200:1500])
def std_fringe(im):
    return np.std(im[200, 1200:1500])

def power_fringe(im):
    mean_linecut = im[200, 1200:1500]
    N = 1*300
    yf = fft(mean_linecut)
    xf = fftfreq(N)

    #Mask out the region around zero
    masked_freqs = ma.masked_inside(xf,-0.015, 0.015)
    frequency = np.abs(xf[~masked_freqs.mask])
    power = np.abs(yf[~masked_freqs.mask])
    #Quantify
    return frequency[np.argmax(power)], np.max(power)


divim_fringe_power, divim_fringe_freq = power_fringe(div_im)
#%% Shift the images in xy
im_new = xy_shift(images[im1_num], 1, 0, x_up=False, y_up=True)
div_im = np.divide(im_new,images[im2_num])

shifted_fringing_std = std_fringe(div_im)
shifted_fringing_power, shifted_fringing_freq = power_fringe(div_im)
#%% Rotate the images
im_new = rotate(images[0], 1, True)
div_im = np.divide(im_new,images[im2_num])
rotated_fringing_std = std_fringe(div_im)
rotated_fringing_power, rotated_fringing_freq = power_fringe(div_im)
#%% Shift AND Rotate
im_new = xy_shift(images[im1_num], 1, 0, x_up=True, y_up=True)
im_new = rotate(im_new, 0.01, True)
div_im = np.divide(im_new,images[im2_num])

new_fringing_std = std_fringe(div_im)
new_fringing_power, rotated_fringing_freq = power_fringe(div_im)
#%% Make a map of fringing based on altitude-rotation coordinates
header_values = ['AZIMUTH', 'ALTITUDE', 'ROTANGLE']
headers = [fits.open(f)[0].header for f in files]
sorting_idx = np.argsort([header['UTSTART'] for header in headers])

altitude = [header['ALTITUDE'] for header in headers]
azimuth = [header['AZIMUTH'] for header in headers]
rotangle = [header['ROTANGLE'] for header in headers]
#%% altitude-rotangle plot
fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(10,10))
plot = ax.scatter(np.array(rotangle)*np.pi/180, altitude, c = fringe_rms, alpha=0.8)
ax.set_title(r'Fringe $\sigma$ with altitude-rotangle', fontsize=20)
ax.set_xlabel(r'Altitude ($\degree$)', labelpad=5, fontsize=20)
#---raxis formatting---
ax.set_rgrids((20,40,60,80))
#---theta axis formatting---
ax.set_theta_zero_location("N", offset=-36)
ax.set_thetamax(40)
ax.set_thetamin(30)
ax.set_thetagrids((30, 35, 40))
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'Fringing $\sigma$', rotation=270, labelpad=25, fontsize=16)
plt.show()
#%% alt-az plot
fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(10,10))
plot = ax.scatter(np.array(azimuth)*np.pi/180, altitude, c = fringe_rms, s = rotangle, alpha=0.8)
ax.set_title(r'Fringe $\sigma$ with altitude-azimuth', fontsize=20)
ax.set_xlabel(r'Altitude ($\degree$)', labelpad=5, fontsize=20)
ax.set_rgrids((20,40,60,80))
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'Fringing $\sigma$', rotation=270, labelpad=25, fontsize=16)
plt.show()
#%% rotangle-azimuth plot
fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(10,10))
plot = ax.scatter(np.array(rotangle)*np.pi/180, azimuth, c = fringe_rms, alpha=0.8)
ax.set_title(r'Fringe $\sigma$ with rotangle-azimuth', fontsize=20)
ax.set_xlabel(r'Azimuth ($\degree$)', labelpad=5, fontsize=20)
#ax.set_rgrids((20,40,60,80))
#---theta axis formatting---
ax.set_theta_zero_location("N", offset=-36)
ax.set_thetamax(40)
ax.set_thetamin(30)
ax.set_thetagrids((30, 35, 40))
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'Fringing $\sigma$', rotation=270, labelpad=25, fontsize=16)
plt.show()