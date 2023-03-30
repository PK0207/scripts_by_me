#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 13:46:51 2023

@author: pkottapalli
"""
from glob import glob
from astropy.io import fits
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy.stats import sigma_clipped_stats
#%% Grab data, separate by aperture width and exposure time
files = glob('xyrot_tweaked_data/ogg/en06/2022*/processed/*.fz', recursive=True)
files = sorted(files)

headers = [fits.open(f)['SCI'].header for f in files]

aperwid = np.array([header['APERWID'] for header in headers])
exptime = np.array([header['EXPTIME'] for header in headers])
aperwid_2 = np.where(aperwid == 2)
exp_80 = np.where([np.isclose(t, 80, atol = 1) for t in exptime])
exp_40 = np.where([np.isclose(t, 40, atol = 1) for t in exptime])
exp80_and_aperwid2 = np.intersect1d(aperwid_2, exp_80)
exp40_and_aperwid2 = np.intersect1d(aperwid_2, exp_40)
use_headers = [headers[i] for i in exp40_and_aperwid2]
use_files = [files[i] for i in exp40_and_aperwid2]
times = np.array([datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f') for header in use_headers])
rotangle = np.array([header['ROTANGLE'] for header in use_headers])
altitude = np.array([header['ALTITUDE'] for header in use_headers])
azimuth = np.array([header['AZIMUTH'] for header in use_headers])
xshift = np.array([header['ORDXSHFT'] for header in use_headers])
yshift = np.array([header['ORDYSHFT'] for header in use_headers])
rotation = np.array([header['ORDROT'] for header in use_headers])
ccdtemp = np.array([header['CCDATEMP'] for header in use_headers])
wmstemp = np.array([header['WMSTEMP'] for header in use_headers])
#%% Altitude Azimuth plot with yshift
fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(10,10))
plot = ax.scatter(np.array(azimuth)*np.pi/180, altitude, c = yshift, s = rotangle, alpha=0.8)
ax.set_title(r'Best shift with altitude-azimuth', fontsize=20)
ax.set_xlabel(r'Altitude ($\degree$)', labelpad=5, fontsize=20)
ax.set_rgrids((20,40,60,80))
ax.set_rmin(90)
ax.set_rmax(0)
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'Y shift', rotation=270, labelpad=25, fontsize=16)
plt.show()
#%% Three plots showing variation of parameters with time
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, dpi=200, sharex=True, figsize=(21,5))
ax1.scatter(times, rotation, s=9)
locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)
ax1.xaxis.set_major_locator(locator)
ax1.xaxis.set_major_formatter(formatter)
ax1.set_xlabel('Time')
ax1.set_ylabel('ORDER ROTATION')

ax2.scatter(times, yshift, s=9)
ax2.xaxis.set_major_locator(locator)
ax2.xaxis.set_major_formatter(formatter)
ax2.set_xlabel('Time')
ax2.set_ylabel('YSHIFT')

rot_plot = ax3.scatter(times, xshift, s=9)
ax3.xaxis.set_major_locator(locator)
ax3.xaxis.set_major_formatter(formatter)
ax3.set_xlabel('Time')
ax3.set_ylabel('XSHIFT')

fig.suptitle('Rotation, Y, and X shift')
plt.figtext(0.7, 0.95, f'N images = {len(times)}, Aperwid = 2, Exptime = 40')
fig.show()
#%% Brightness of the image vs the parameters
brightness = [sigma_clipped_stats(fits.open(f)['SCI'].data, sigma = 5)[1] for f in use_files]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, dpi=200, sharex=True, figsize=(21,5))
ax1.scatter(brightness, rotation, s=9)
ax1.set_xlabel('Median of Image')
ax1.set_ylabel('ORDER ROTATION')

ax2.scatter(brightness, yshift, s=9)
ax2.set_xlabel('Median of Image')
ax2.set_ylabel('YSHIFT')

ax3.scatter(brightness, xshift, s=9)
ax3.set_xlabel('Median of Image')
ax3.set_ylabel('XSHIFT')

fig.suptitle('Rotation, Y, and X shift')
plt.figtext(0.7, 0.95, f'N images = {len(times)}, Aperwid = 2, Exptime = 40')
fig.show()
#%% X and Y shift plotted on the same plot
fig, ax1 = plt.subplots(dpi=200)
ax1.scatter(times, xshift, s=5, label='xshift')
ax1.scatter(times, yshift, s=5, label='yshift')
ax1.hlines(np.mean(yshift), times[0], times[-1], 'r', '--', label=f'Mean y-shift = {np.mean(yshift): 0.2f}')
ax1.hlines(np.mean(xshift), times[0], times[-1], 'r', '--', label=f'Mean x-shift = {np.mean(xshift): 0.2f}')
locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)
ax1.xaxis.set_major_locator(locator)
ax1.xaxis.set_major_formatter(formatter)
ax1.set_xlabel('Time')
ax1.set_ylabel('Shift relative to skyflat (pixels)')
ax1.legend()
fig.show()
#%%
fig, ax1 = plt.subplots(dpi=200)
ax1.scatter(xshift, yshift, s=5)
locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)
ax1.xaxis.set_major_locator(locator)
ax1.set_xlabel('XSHIFT (pixel)')
ax1.set_ylabel('YSHIFT (pixel)')
fig.show()
#%% GIF of shifted flats
import matplotlib.animation as animation
from scipy.ndimage import shift
from PIL import Image
def rotate(data, angle=0, clockwise=False):
    pil_image = Image.fromarray(data)
    if clockwise==False:
        new_im = np.array(pil_image.rotate(360-angle))
    else:
        new_im = np.array(pil_image.rotate(angle))
    return new_im

def std_fringe(im):
    return np.std(im[140:190, 1100:1300])

def std_align(im):
    return np.std(im[10:150, 250:500])

fig, ax = plt.subplots(dpi = 200)
template = fits.open(use_files[0])['SCI'].data
template = shift(template, (yshift[0], 0), mode='wrap')
first_im = fits.open(use_files[1])['SCI'].data
image = shift(first_im, (yshift[1], 0), mode='wrap')
#image = rotate(image, rotation[1]-rotation[0], clockwise=False)
image /= template
im1 = ax.imshow(image, vmin = np.median(first_im)-3*np.std(first_im), vmax = np.median(first_im)+3*np.std(first_im), origin='lower')
plt.colorbar(im1)
ax.set_title(f'Image No. 1 flatfielded by a template after shifting by {yshift[1]:0.2f}', fontsize=8)
def init():
    return im1,
    
def update(i):
    i+=1
    image = fits.open(use_files[i])['SCI'].data
    image = shift(image, (yshift[i], 0), mode='wrap')
    #image = rotate(image, rotation[i]-rotation[0], clockwise=False)
    image /= template
    #median = sigma_clipped_stats(image, sigma = 5)[1]
    #std = np.std(image)
    im1 = ax.imshow(image[140:190, 1100:1300], vmin = np.median(first_im)-3*np.std(first_im), vmax = np.median(first_im)+3*np.std(first_im), origin='lower')
    #ax.set_title(f'Image No. {i} flatfielded by a template shifted: (rot: {rotation[i]:0.2f}, y: {yshift[i]:0.2f}, x: {xshift[i]:0.2f})', fontsize=8)
    ax.set_title(f'Image No. {i} flatfielded by a template shifted: (y: {yshift[i]:0.2f})', fontsize=8)
    ax.set_xlabel(f'Quality of reduction: {std_fringe(image):0.2f}')
    return im1, 

anim = animation.FuncAnimation(fig, update, frames=np.arange(50), init_func=init, interval = 50, blit = True)
anim.save('shift_fringes_animation.gif', fps=2)
#%% GIF of non-shifted flats
fig, ax = plt.subplots(dpi = 200)
template = fits.open(use_files[0])['SCI'].data
first_im = fits.open(use_files[1])['SCI'].data
first_im /= template
im1 = ax.imshow(first_im, vmin = np.median(first_im)-3*np.std(first_im), vmax = np.median(first_im)+3*np.std(first_im), origin='lower')
plt.colorbar(im1)
ax.set_title('Image No. 1 flatfielded by a template without shifting', fontsize=8)
def init():
    return im1,
    
def update(i):
    i+=1
    image = fits.open(use_files[i])['SCI'].data
    image /= template
    #median = sigma_clipped_stats(image, sigma = 5)[1]
    #std = np.std(image)
    im1 = ax.imshow(image[0:150, 250:500], vmin = np.median(first_im)-3*np.std(first_im), vmax = np.median(first_im)+3*np.std(first_im), origin='lower')
    ax.set_title(f'Image No. {i} flatfielded by a template without shifting', fontsize=8)
    ax.set_xlabel(f'Quality of reduction: {std_align(image):0.2f}')
    return im1, 

anim = animation.FuncAnimation(fig, update, frames=np.arange(50), init_func=init, interval = 50, blit = True)
anim.save('non_shift_edges_animation.gif', fps=2)
#%% Shift vs. quality
template = fits.open(use_files[0])['SCI'].data
template = shift(template, (yshift[0], 0), mode='wrap')
no_shift_quality = []
shift_quality = []
for i, f in enumerate(use_files):
    image = fits.open(f)['SCI'].data
    #image = shift(image, (yshift[i]-yshift[0], xshift[i]-xshift[0]), mode='wrap')
    #image = rotate(image, rotation[i]-rotation[0], clockwise=False)
    image /= template
    #no_shift_quality.append(sigma_clipped_stats(image, sigma = 5)[0])
    no_shift_quality.append(std_align(image))
    

for i, f in enumerate(use_files):
    image = fits.open(f)['SCI'].data
    image = shift(image, (yshift[i], 0), mode='wrap')
    #image = rotate(image, rotation[i]-rotation[0], clockwise=False)
    image /= template
    #shift_quality.append(sigma_clipped_stats(image, sigma = 5)[0])
    shift_quality.append(std_align(image))
#shift_dist = np.sqrt((yshift-yshift[0])**2 + (rotation-rotation[0])**2 + (xshift-xshift[0])**2)
shift_dist = np.sqrt((yshift-yshift[0])**2)
#%%
times = np.array([header['MJD-OBS'] for header in use_headers])
#Plot
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, dpi=200, figsize=(10,5))
ax1.scatter(yshift, shift_quality, s=7, c=times, alpha=0.8)
ax1.set_title(r'$\sigma$ of alignment edges vs. shift distance')
#ax1.set_title(r'Mean of image vs. shift distance')
#ax1.set_xlabel(r'Shift Distance $\sqrt{(y-y_0)^2+(x-x_0)^2+(rot-rot_0)^2}$')
ax1.set_xlabel('Y shift')
ax1.set_ylabel('Quality')

plot = ax2.scatter(yshift, no_shift_quality, s=7, c=times, alpha=0.8)
ax2.set_title(r'$\sigma$ of alignment edges without applying shift')
#ax2.set_title(r'Mean of image without applying shift')
ax2.set_xlabel('Y shift')
#ax1.set_ylabel('Quality')
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'Image MJD', rotation=270, labelpad=20)
fig.show()