#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 10:29:22 2023

@author: pkottapalli
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from datetime import datetime
#%%
files = sorted(glob('New_AltAz_data/ogg2m001-*w00.fits.fz', recursive=True))
images = [fits.open(f)['SCI'].data for f in files]
headers = [fits.open(f)['SCI'].header for f in files]
data_type = 'red flat'
#%%
def normalize(data):
    min_data = min(data)
    max_data = max(data)
    new_data = [(i-min_data)/(max_data-min_data) for i in data]
    return new_data

for image in images[:10]:
    red_fringe_region = image[140:190, 1100:1300]
    red_fringe_linecut = np.mean(red_fringe_region, axis=0)
    plt.figure(dpi=200)
    plt.plot(np.linspace(1100,1300, 200), red_fringe_linecut)
    plt.title('linecut showing fringing in the red region')
    plt.xlabel('x pixel')
    plt.ylabel('pixel value')
    plt.show()
    blue_region = image[50:125, 250:500]
    blue_linecut = np.mean(blue_region, axis=0)
    plt.figure(dpi=200)
    plt.plot(np.linspace(250, 500, 250), blue_linecut)
    plt.title('linecut in the blue region')
    plt.xlabel('x pixel')
    plt.ylabel('pixel value')
    plt.show()
#%% Sort data into aperture widths and divide by exposure time
aperwid = np.array([header['APERWID'] for header in headers])
exptime = np.array([header['EXPTIME'] for header in headers])
times = np.array([datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f') for header in headers])
aperwid_2 = np.where(aperwid == 2)
aperwid_16 = np.where(aperwid == 1.6)
aperwid_12 = np.where(aperwid == 1.2)
aperwid_6 = np.where(aperwid == 6)
aperwid_2_images = np.array(images)[aperwid_2]
aperwid_2_images = [image/exp for image, exp in zip(aperwid_2_images, exptime)]
#%%
aperwid_2_times = times[aperwid_2]

red_mins = []
red_maxs = []
red_var = []

blue_mins = []
blue_maxs = []
blue_var = []

for image in aperwid_2_images:
    red_fringe_region = image[140:190, 1100:1300]
    red_fringe_linecut = np.mean(red_fringe_region, axis=0)
    red_var.append(np.std(red_fringe_linecut))
    red_maxs.append(max(red_fringe_linecut))
    red_mins.append(min(red_fringe_linecut))
    blue_region = image[50:125, 250:500]
    blue_linecut = np.mean(blue_region, axis=0)
    blue_var.append(np.std(blue_linecut))
    blue_maxs.append(max(blue_linecut))
    blue_mins.append(min(blue_linecut))

plt.figure(dpi=200)
plt.title('Maximum and minimum value of red fringing area in every image')
plt.scatter(aperwid_2_times, red_maxs, label='Maxima', cmap = 'tab20', s=5)
plt.scatter(aperwid_2_times, red_mins, label='Minima', cmap = 'tab20', s=5)
plt.legend()
# cax = plt.colorbar()
# cax.set_ticks([1.2, 1.6, 2.0, 6.0])
plt.xlabel('Date Observation')
plt.ylabel('Pixel value')
plt.xticks(rotation=45)
plt.show()

plt.figure(dpi=200)
plt.title('Maximum and minimum value of blue part of spectrum in every image')
plt.scatter(aperwid_2_times, blue_maxs, label='Maxima', cmap = 'tab20', s=5)
plt.scatter(aperwid_2_times, blue_mins, label='Minima', cmap = 'tab20', s=5)
plt.legend()
# cax = plt.colorbar()
# cax.set_ticks([1.2, 1.6, 2.0, 6.0])
plt.xlabel('Date Observation')
plt.ylabel('Pixel value')
plt.xticks(rotation=45)
plt.show()

plt.figure(dpi=200)
plt.title('Variation of blue part of spectrum in every image')
plt.scatter(aperwid_2_times, blue_var, s=5)
plt.xlabel('Date Observation')
plt.ylabel('Pixel value')
plt.xticks(rotation=45)
plt.show()
#%%
red_lamp_files = sorted(glob('lco_data-20230216-33/ttflat*red*.fits', recursive=True))
rectified_flats = [fits.open(f)[0].data for f in red_lamp_files]
input_fringe = rectified_flats[0]
correction_fringe_frame = rectified_flats[-1]

import sep
input_fringe_background = sep.Background(input_fringe.byteswap().newbyteorder()).back()

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, dpi=200)
inp = ax1.imshow(input_fringe)
ax1.set_title('Input Frame')
ax1.axis('off')
plt.colorbar(inp, location='bottom', ax=ax1)
corr = ax2.imshow(correction_fringe_frame)
plt.colorbar(corr, location='bottom', ax=ax2)
ax2.set_title('Fringe Frame to correct input')
ax2.axis('off')
back = ax3.imshow(input_fringe_background)
ax3.set_title('Background of Input')
plt.colorbar(back, location='bottom', ax=ax3)
ax3.axis('off')
fig.show()
#%%
from scipy.optimize import minimize

def iraf_defringe(s):
    result = np.mean((input_fringe-input_fringe_background - s*(correction_fringe_frame-np.mean(correction_fringe_frame)))*(correction_fringe_frame-np.mean(correction_fringe_frame)))
    return result

res = minimize(iraf_defringe, 0.2816)