#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:37:30 2023

@author: pkottapalli
"""
from glob import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.animation as animation
from scipy.ndimage import shift
#%% Show Xshift yshift and rotation with alt az
files = glob('xy_tweaked_data/ogg/en06/2022*/processed/*.fz', recursive=True)
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
#rotation = np.array([header['ORDROT'] for header in use_headers])
ccdtemp = np.array([header['CCDATEMP'] for header in use_headers])
wmstemp = np.array([header['WMSTEMP'] for header in use_headers])
#%%
template = fits.open(use_files[0])['SCI'].data
fig, (ax1, ax2) = plt.subplots(1, 2, dpi=200, figsize=(12, 5))
im1 = ax1.imshow(template[140:190, 1100:1300], origin='lower', vmin=np.median(template)-5*np.std(template), vmax=np.median(template)+5*np.std(template))
ax1.axis('off')
ax1.set_title('Fringe region of template lampflat')
cax = plt.colorbar(im1)
im2 = ax2.imshow(template[10:150, 250:500], origin='lower', vmin=np.median(template)-5*np.std(template), vmax=np.median(template)+5*np.std(template))
ax2.set_title('Edge region of template lampflat')
ax2.axis('off')
cax = plt.colorbar(im2)
fig.show()
#%% Manually shift image in y animation
def std_fringe(im):
    return np.std(im[140:190, 1100:1300])

def std_align(im):
    return np.std(im[10:150, 250:500])

im_num = 550
shift_arr = np.arange(-1, 1, 0.1)

fig, ax = plt.subplots(1, 3, dpi = 200, figsize=(15,5))
template = fits.open(use_files[0])['SCI'].data
template = shift(template, (yshift[0], xshift[0]))
im1 = ax[0].imshow(template, vmin = 0.8, vmax = 1.2, origin='lower')

ax[0].set_title('full template', fontsize=8)
im2 = ax[1].imshow(template[140:190, 1100:1300], vmin = 0.8, vmax = 1.2, origin='lower')
ax[1].set_title('fringe region', fontsize=8)
im3 = ax[2].imshow(template[10:150, 250:500], vmin = 0.8, vmax = 1.2, origin='lower')
ax[2].set_title('edge region', fontsize=8)
plt.colorbar(im3)

def init():
    return im1, im2, im3

def y_update(i):
        first_im = fits.open(use_files[im_num])['SCI'].data
        shifted_image = shift(first_im, (shift_arr[i], 0))
        shifted_image /= template
        im1 = ax[0].imshow(shifted_image, vmin = 0.8, vmax = 1.2, origin='lower')
        ax[0].set_title(f'yshift={shift_arr[i]:0.2f}')
        im2 = ax[1].imshow(shifted_image[140:190, 1100:1300], vmin = 0.8, vmax = 1.2, origin='lower')
        ax[1].set_title(f'fringe region: yshift={shift_arr[i]:0.2f}', fontsize=8)
        im3 = ax[2].imshow(shifted_image[10:150, 250:500], vmin = 0.8, vmax = 1.2, origin='lower')
        ax[2].set_title(f'edge region: yshift={shift_arr[i]:0.2f}', fontsize=8)
        return im1, im2, im3

anim = animation.FuncAnimation(fig, y_update, frames=len(shift_arr), init_func=init, interval = 50, blit = True)
anim.save('manual_yshift.gif', fps=2)
#%% Manually shift in y and plot edge and fringe parameter minimization
fringe_amt = []
edge_amt = []
#edge_region = [10:150, 250:500]
#fringe_region = [140:190, 1100:1300]
shift_arr = np.arange(-2, 2, 0.1)
for i, y_shift_amt in enumerate(shift_arr):
    first_im = fits.open(use_files[im_num])['SCI'].data
    shifted_image = shift(first_im, (y_shift_amt, 0))
    shifted_image /= template
    fringe_amt.append(std_fringe(shifted_image))
    edge_amt.append(std_align(shifted_image))

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=200)
ax1.plot(shift_arr, fringe_amt)
ax1.set_ylabel(r'Fringe region $\sigma$')
ax2.plot(shift_arr, edge_amt)
ax2.set_ylabel(r'Edge region $\sigma$')
ax2.set_xlabel('Shift in y')
fig.show()

print(f'Best Edge removal:{shift_arr[np.argmin(np.abs(edge_amt))]}')
print(f'Best Fringe removal:{shift_arr[np.argmin(np.abs(fringe_amt))]}')
#%% Manually shift image in x
im_num = 550
shift_arr = np.arange(-3, 3, 0.1)

def x_update(i):
        first_im = fits.open(use_files[im_num])['SCI'].data
        shifted_image = shift(first_im, (0, shift_arr[i]))
        shifted_image /= template
        im1 = ax[0].imshow(shifted_image, vmin = 0.8, vmax = 1.2, origin='lower')
        ax[0].set_title(f'xshift={shift_arr[i]:0.2f}')
        im2 = ax[1].imshow(shifted_image[140:190, 1100:1300], vmin = 0.8, vmax = 1.2, origin='lower')
        ax[1].set_title(f'fringe region: xshift={shift_arr[i]:0.2f}', fontsize=8)
        im3 = ax[2].imshow(shifted_image[10:150, 250:500], vmin = 0.8, vmax = 1.2, origin='lower')
        ax[2].set_title(f'edge region: xshift={shift_arr[i]:0.2f}', fontsize=8)
        return im1, im2, im3

anim = animation.FuncAnimation(fig, x_update, frames=len(shift_arr), init_func=init, interval = 50, blit = True)
anim.save('manual_xshift.gif', fps=2)
#%% Manually shift in x and plot edge and fringe parameter minimization
fringe_amt = []
edge_amt = []
#edge_region = [10:150, 250:500]
#fringe_region = [140:190, 1100:1300]
shift_arr = np.arange(-2, 2, 0.1)
for i, x_shift_amt in enumerate(shift_arr):
    first_im = fits.open(use_files[im_num])['SCI'].data
    shifted_image = shift(first_im, (0, x_shift_amt))
    shifted_image /= template
    fringe_amt.append(std_fringe(shifted_image))
    edge_amt.append(std_align(shifted_image))

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=200)
ax1.plot(shift_arr, fringe_amt)
ax1.set_ylabel(r'Fringe region $\sigma$')
ax2.plot(shift_arr, edge_amt)
ax2.set_ylabel(r'Edge region $\sigma$')
ax2.set_xlabel('Shift in x')
fig.show()

print(f'Best Edge removal:{shift_arr[np.argmin(np.abs(edge_amt))]}')
print(f'Best Fringe removal:{shift_arr[np.argmin(np.abs(fringe_amt))]}')
#%% Manually shift image in x and y
im_num = 14
template = fits.open(use_files[0])['SCI'].data
template = shift(template, (yshift[0], xshift[0]))
first_im = fits.open(use_files[im_num])['SCI'].data
fringe_amt = []
edge_amt = []
y_shift_arr = np.arange(-2, 2, 0.1)
x_shift_arr = np.arange(-2, 2, 0.1)
fringe_grid = np.ones((len(y_shift_arr), len(x_shift_arr)))
align_grid = np.ones((len(y_shift_arr), len(x_shift_arr)))
for i, y_shift_amt in enumerate(y_shift_arr):
    for j, x_shift_amt in enumerate(x_shift_arr):
        shifted_image = shift(first_im, (y_shift_amt, x_shift_amt))
        shifted_image /= template
        fringe_grid[i, j] = std_fringe(shifted_image)
        align_grid[i, j] = std_align(shifted_image)
#%% Gridsearch heatmap of minimization values
import matplotlib.colors as colors
fig, (ax1, ax2) = plt.subplots(1, 2, dpi=200, figsize=(12, 5))
im1 = ax1.imshow(fringe_grid, origin='lower', norm=colors.LogNorm(vmin=fringe_grid.min()+0.001, vmax=fringe_grid.max()))
ax1.set_xlabel('X Shift')
ax1.set_ylabel('Y Shift')
ax1.set_xticks(np.arange(0, len(x_shift_arr), 1), labels = [str(round(shift, 2)) for shift in x_shift_arr], rotation=90, fontsize=6)
ax1.set_yticks(np.arange(0, len(y_shift_arr), 1), labels = [str(round(shift, 2)) for shift in y_shift_arr], fontsize=6)
ax1.set_title('Fringe region grid search')
cax = plt.colorbar(im1)
cax.ax.set_ylabel('Fringe region $\sigma$', labelpad=15, rotation=270)
im2 = ax2.imshow(align_grid, origin='lower', norm=colors.LogNorm(vmin=align_grid.min(), vmax=align_grid.max()))
ax2.set_xlabel('X Shift')
ax2.set_ylabel('Y Shift')
ax2.set_xticks(np.arange(0, len(x_shift_arr), 1), labels = [str(round(shift, 2)) for shift in x_shift_arr], rotation=90, fontsize=6)
ax2.set_yticks(np.arange(0, len(y_shift_arr), 1), labels = [str(round(shift, 2)) for shift in y_shift_arr], fontsize=6)
ax2.set_title('Edge region grid search')
cax = plt.colorbar(im2)
cax.ax.set_ylabel('Edge region $\sigma$', labelpad=15, rotation=270)
fig.suptitle(f'Image no.: {im_num}, Altitude dist.:{altitude[im_num]- altitude[0]: 0.2f}, Azimuth dist.:{azimuth[im_num]- azimuth[0]: 0.2f}', fontsize=15)
fig.show()
#%% Show what best shift looks like
fringe_best_y, fringe_best_x = y_shift_arr[np.where(fringe_grid == np.min(fringe_grid))[0]], x_shift_arr[np.where(fringe_grid == np.min(fringe_grid))[1]]
edge_best_y, edge_best_x = y_shift_arr[np.where(align_grid == np.min(align_grid))[0]], x_shift_arr[np.where(align_grid == np.min(align_grid))[1]]
best_y_idx = np.where(align_grid == np.min(align_grid))[0]
best_y, best_x = y_shift_arr[best_y_idx], x_shift_arr[np.where(fringe_grid[best_y_idx, :] == np.min(fringe_grid[best_y_idx, :]))[1]]

first_im = fits.open(use_files[im_num])['SCI'].data
fringe_shifted_image = shift(first_im, (fringe_best_y, fringe_best_x))
fringe_shifted_image /= template

edge_shifted_image = shift(first_im, (edge_best_y, edge_best_x))
edge_shifted_image /= template

best_of_both_worlds = shift(first_im, (best_y, best_x))
best_of_both_worlds /= template

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, dpi=200, figsize=(20, 10))
ax1.imshow(fringe_shifted_image, origin='lower', vmin=0.8, vmax=1.2)
ax1.set_title('Best fringe removal')
ax1.set_ylabel('Full Image')
ax2.imshow(edge_shifted_image, vmin=0.8, vmax=1.2, origin='lower')
ax2.set_title('Best Edge removal')
im3 = ax3.imshow(best_of_both_worlds, vmin=0.8, vmax=1.2, origin='lower')
ax3.set_title('best x for fringe along best y for edge')
plt.colorbar(im3)

ax4.imshow(fringe_shifted_image[140:190, 1100:1300], origin='lower', vmin=0.8, vmax=1.2)
ax4.set_ylabel('Fringe Region')
ax5.imshow(edge_shifted_image[140:190, 1100:1300], vmin=0.8, vmax=1.2, origin='lower')
im6 = ax6.imshow(best_of_both_worlds[140:190, 1100:1300], vmin=0.8, vmax=1.2, origin='lower')
plt.colorbar(im6)

ax7.imshow(fringe_shifted_image[10:150, 250:500], origin='lower', vmin=0.8, vmax=1.2)
ax7.set_ylabel('Edge Region')
ax7.set_xlabel(f'xshift: {fringe_best_x[0]:0.2f}, yshift: {fringe_best_y[0]:0.2f}', fontsize=20)
ax8.imshow(edge_shifted_image[10:150, 250:500], vmin=0.8, vmax=1.2, origin='lower')
ax8.set_xlabel(f'xshift: {edge_best_x[0]:0.2f}, yshift: {edge_best_y[0]:0.2f}', fontsize=20)
im9 = ax9.imshow(best_of_both_worlds[10:150, 250:500], vmin=0.8, vmax=1.2, origin='lower')
ax9.set_xlabel(f'xshift: {best_x[0]:0.2f}, yshift: {best_y[0]:0.2f}', fontsize=20)
plt.colorbar(im9)

fig.suptitle(f'Image no.: {im_num}, Altitude dist.:{altitude[im_num]- altitude[0]: 0.2f}, Azimuth dist.:{azimuth[im_num]- azimuth[0]: 0.2f}', fontsize=20)
fig.show()
#%%
plt.figure(dpi=200)
shifted_image = shift(first_im, (yshift[im_num], xshift[im_num]))
shifted_image /= template
plt.imshow(shifted_image, origin='lower', vmin = 0.8, vmax = 1.2)
plt.show()