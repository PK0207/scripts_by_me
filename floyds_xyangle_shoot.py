#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:07:44 2023

@author: pkottapalli
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
#%%
#Non-rectified spectrum
#files = sorted(glob('lco_data-20230216-33/lco_data-20230216-33/flat*red_*.fits', recursive=True))
files = sorted(glob('New_AltAz_data/*.fits.fz', recursive=True))
images = []
for f in files:
    try:
        images.append(fits.open(f))
    except OSError:
        print(f)
data_type = 'red flat'
#%%
for i, im in enumerate(images):
    plt.imshow(im, origin='lower')
    plt.colorbar()
    plt.title(f'{data_type}, {i}')
    plt.show()
#%%
#alt az
import matplotlib.dates as mdates
from datetime import datetime

header_values = ['AZIMUTH', 'ALTITUDE', 'ROTANGLE']
headers = [fits.open(f)['SCI'].header for f in files]
sorting_idx = np.argsort([header['DATE-OBS'] for header in headers])

env_parameters = []
for header in headers:
    env_parameters.append([header[key] for key in header_values])

env_parameters = np.array(env_parameters)

#plt.style.use('seaborn')
fig, axes = plt.subplots(1, 3, sharex=True, figsize=(35,10))
fig.tight_layout(pad = 12)
times = np.sort([datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f') for header in headers])
t0 = times[0]
t1 = times[-1]
for ax, i in zip(axes.flatten(), range(len(env_parameters[0]))):
    ax.scatter(times, env_parameters[:,i][sorting_idx], s=15)
    ax.set_ylabel(header_values[i], fontsize=20)
    ax.set_xlabel('Obs time', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=14, size=10)
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    #ax.xaxis.set_major_formatter(formatter)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
fig.show()
#%%
im1_num = 0
im2_num = 1
alt1, az1, rot1 = env_parameters[im1_num,1], env_parameters[im1_num,0], env_parameters[im1_num,2]
alt2, az2, rot2 = env_parameters[im2_num,1], env_parameters[im2_num,0], env_parameters[im2_num,2]
#divide one flat by another without shifting
div_im = np.divide(images[im1_num],images[im2_num])
median = np.median(div_im)
std = np.std(div_im)
fig = plt.figure(figsize=(7.5,5), dpi=200)
plt.imshow(div_im, origin='lower', vmin=median-2*std, vmax = median+2*std)
plt.title('Two flats divided')
plt.xlabel(f'ALTAZROT1 : ({alt1:.2f}, {az1:.2f}, {rot1:0.2f}), ALTAZROT2 : ({alt2:.2f}, {az2:.2f}, {rot2:0.2f})')
plt.colorbar()
fig.show()
#%%
#And then shift it a little
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

im1_num = 0
im2_num = 1
alt1, az1, rot1 = env_parameters[im1_num,1], env_parameters[im1_num,0], env_parameters[im1_num,2]
alt2, az2, rot2 = env_parameters[im2_num,1], env_parameters[im2_num,0], env_parameters[im2_num,2]

im_new = xy_shift(images[im1_num], 1, 0, x_up=False, y_up=True)
div_im = np.divide(im_new,images[im2_num])

median = np.median(div_im)
std = np.std(div_im)
fig = plt.figure(figsize=(7.5,5), dpi=200)
plt.imshow(div_im, origin='lower', vmin=median-2*std, vmax = median+2*std)
plt.title('Two flats divided, shifted left 1 in x')
plt.xlabel(f'ALTAZROT1 : ({alt1:.2f}, {az1:.2f}, {rot1:0.2f}), ALTAZROT2 : ({alt2:.2f}, {az2:.2f}, {rot2:0.2f})')
plt.colorbar()
fig.show()
#%%
#Rotate
# from scipy.spatial.transform import Rotation as R
# r = R.from_euler('z', 1, degrees=True).as_matrix()
# #r.as_euler()
# im_new = np.matmul([images[0]], r)
# plt.imshow(im_new)
# plt.show()
from PIL import Image
def rotate(data, angle=0, clockwise=False):
    pil_image = Image.fromarray(data)
    if clockwise==False:
        new_im = np.array(pil_image.rotate(360-angle))
    else:
        new_im = np.array(pil_image.rotate(angle))
    return new_im

im_new = rotate(images[0], 1, True)
div_im = np.divide(im_new,images[im2_num])

median = np.median(div_im)
std = np.std(div_im)
fig = plt.figure(figsize=(7.5,5), dpi=200)
plt.imshow(div_im, origin='lower', vmin=median-2*std, vmax = median+2*std)
plt.title('Two flats divided, one rotated by 1 degree')
plt.xlabel(f'ALTAZROT1 : ({alt1:.2f}, {az1:.2f}, {rot1:0.2f}), ALTAZROT2 : ({alt2:.2f}, {az2:.2f}, {rot2:0.2f})')
plt.colorbar()
fig.show()
#%% Shift AND Rotate
im_new = xy_shift(images[im1_num], 1, 0, x_up=True, y_up=True)
im_new = rotate(im_new, 0.01, True)
div_im = np.divide(im_new,images[im2_num])

median = np.median(div_im)
std = np.std(div_im)
fig = plt.figure(figsize=(7.5,5), dpi=200)
plt.imshow(div_im, origin='lower', vmin=median-1*std, vmax = median+1*std)
plt.title('Flatted flat, one rotated by 0.01, shifted right 1')
plt.xlabel(f'ALTAZROT1 : ({alt1:.2f}, {az1:.2f}, {rot1:0.2f}), ALTAZROT2 : ({alt2:.2f}, {az2:.2f}, {rot2:0.2f})')
plt.colorbar()
fig.show()
#%% Shift and rotate with very different altazrot
im_new = xy_shift(images[0], 1, 1, x_up=False, y_up=True)
im_new = rotate(im_new, 0.1, True)
div_im = np.divide(im_new,images[-1])

median = np.median(div_im)
std = np.std(div_im)
fig = plt.figure(figsize=(7.5,5), dpi=200)
plt.imshow(div_im, origin='lower', vmin=median-1*std, vmax = median+1*std)
plt.title('Flatted flat, one rotated by 0.1, shifted up and left 1')
plt.xlabel(f'ALTAZROT1 : ({alt1:.2f}, {az1:.2f}, {rot1:0.2f}), ALTAZROT2 : ({alt2:.2f}, {az2:.2f}, {rot2:0.2f})')
plt.colorbar()
fig.show()
#%% Distort image for template matching
