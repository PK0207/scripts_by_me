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
rotangle = np.array([header['ROTANGLE'] for header in headers])
altitude = np.array([header['ALTITUDE'] for header in headers])
azimuth = np.array([header['AZIMUTH'] for header in headers])

times = np.array([datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f') for header in headers])
exp_80 = np.where([np.isclose(time, 80, atol = 1) for time in exptime])
aperwid_2 = np.where(aperwid == 2)
aperwid_16 = np.where(aperwid == 1.6)
aperwid_12 = np.where(aperwid == 1.2)
aperwid_6 = np.where(aperwid == 6)
exp80_and_aperwid2 = np.intersect1d(aperwid_2, exp_80)
aperwid_2_images = np.array(images)[exp80_and_aperwid2]
#%%
aperwid_2_times = times[exp80_and_aperwid2]

red_mins = []
red_maxs = []
red_var = []

blue_means = []
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
    blue_means.append(np.mean(blue_linecut))

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
plt.title('Mean value of blue part of spectrum in every image')
plt.scatter(aperwid_2_times, blue_means, label='Maxima', cmap = 'tab20', s=5)
#plt.scatter(aperwid_2_times, blue_mins, label='Minima', cmap = 'tab20', s=5)
#plt.legend()
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
fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(10,10))
plot = ax.scatter(np.array(rotangle[exp80_and_aperwid2])*np.pi/180, altitude[exp80_and_aperwid2], c = red_var, s=7, alpha=0.8)
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
plot = ax.scatter(np.array(azimuth[exp80_and_aperwid2])*np.pi/180, altitude[exp80_and_aperwid2], c = red_var, s = rotangle[exp80_and_aperwid2], alpha=0.8)
ax.set_title(r'Fringe $\sigma$ with altitude-azimuth', fontsize=20)
ax.set_xlabel(r'Altitude ($\degree$)', labelpad=5, fontsize=20)
ax.set_rgrids((20,40,60,80))
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'Fringing $\sigma$', rotation=270, labelpad=25, fontsize=16)
plt.show()
#%% rotangle-azimuth plot
fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(10,10))
plot = ax.scatter(np.array(rotangle[exp80_and_aperwid2])*np.pi/180, azimuth[exp80_and_aperwid2], c = red_var, s=7, alpha=0.8)
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
#%% IRAF defringing
red_lamp_files = sorted(glob('lco_data-20230216-33/lco_data-20230216-33/ttflat*red*.fits', recursive=True))
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
from scipy.optimize import minimize_scalar

def iraf_defringe(s):
    result = np.mean((input_fringe-input_fringe_background - s*(correction_fringe_frame-np.mean(correction_fringe_frame)))*(correction_fringe_frame-np.mean(correction_fringe_frame)))
    return result

res = minimize_scalar(iraf_defringe, bounds = (0.28, 0.3), method='bounded')