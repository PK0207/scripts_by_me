# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 10:10:29 2023

@author: prera
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from scipy.fft import fft, fftfreq
import numpy.ma as ma
from skimage.feature import match_template
#%%
files = sorted(glob('New_AltAz_data/ogg2m001-*w00.fits.fz', recursive=True))
images = [fits.open(f)['SCI'].data for f in files]
data_type = 'red flat'
#%% Get headers fringing based on altitude-rotation coordinates
header_values = ['AZIMUTH', 'ALTITUDE', 'ROTANGLE']
headers = [fits.open(f)['SCI'].header for f in files]
sorting_idx = np.argsort([header['UTSTART'] for header in headers])

altitude = [header['ALTITUDE'] for header in headers]
azimuth = [header['AZIMUTH'] for header in headers]
rotangle = [header['ROTANGLE'] for header in headers]
#%% shift in x and y until fringing is minimized, then map out shift on sky
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
    return np.std(im[165, 1100:1300])

def power_fringe(im):
    mean_linecut = im[165, 1100:1300]
    N = 1*300
    yf = fft(mean_linecut)
    xf = fftfreq(N)

    #Mask out the region around zero
    masked_freqs = ma.masked_inside(xf,-0.015, 0.015)
    frequency = np.abs(xf[~masked_freqs.mask])
    power = np.abs(yf[~masked_freqs.mask])
    #Quantify
    return frequency[np.argmax(power)], np.max(power)
#Shift in every direction by 1 pixel
shift_num = 4

shift_dict = {'north': xy_shift(images[im1_num], 0, shift_num, y_up=True), 
              'northeast': xy_shift(images[im1_num],shift_num,shift_num, x_up=True, y_up=True), 
              'east': xy_shift(images[im1_num],shift_num, 0, x_up=True), 
              'southeast': xy_shift(images[im1_num],shift_num,shift_num, x_up=True, y_up=False),
              'south': xy_shift(images[im1_num], 0,shift_num, y_up=False),
              'southwest': xy_shift(images[im1_num],shift_num,shift_num, x_up=False, y_up=False),
              'west': xy_shift(images[im1_num],shift_num, 0, x_up=False),
              'northwest': xy_shift(images[im1_num],shift_num,shift_num, x_up=False, y_up=True), 
              'none': images[im1_num]}

alt1 = altitude[0]
az1 = azimuth[0]
rot1 = rotangle[0]
#Shift every image relative to the first image, in every direction
best_shift = []
for i in range(len(images[:1])):
    im1_num = i
    alt2 = altitude[i]
    az2 = azimuth[i]
    rot2 = rotangle[i]
    fringe_amount = []
    for shift in shift_dict:
        shifted_image = shift_dict[shift]
        div_im = np.divide(images[0], shifted_image)
        fringe_amount.append(std_fringe(div_im))
        #Plot fringe section
        plt.imshow(div_im[140:190, 1200:1500], origin='lower')
        plt.colorbar()
        plt.title(f'{data_type}, {i} fringe region')
        plt.show()
        #Plot linecut
        plt.plot(np.arange(200), div_im[165, 1100:1300])
        plt.title(f'{data_type}, {i} fringe region linecut')
        plt.xlabel('x pixel')
        plt.ylabel('pixel value')
        plt.show()
        #Plot divided images
        median = np.median(div_im)
        std = np.std(div_im)
        fig = plt.figure(figsize=(7.5,5), dpi=200)
        plt.imshow(div_im, origin='lower', vmin=median-2*std, vmax = median+2*std)
        plt.title(f'Two flats divided, shifted {shift}')
        plt.xlabel(f'ALTAZROT1 : ({alt1:.2f}, {az1:.2f}, {rot1:0.2f}), ALTAZROT2 : ({alt2:.2f}, {az2:.2f}, {rot2:0.2f})')
        plt.colorbar()
        plt.show()
    best_shift.append(list(shift_dict)[np.argmin(fringe_amount)])
#%% sklearn template matching
fringe_region1 = images[0][140:190, 1200:1500]
similarity = []
for image in images:
    fringe_region2 = image[140:190, 1200:1500]
    result = match_template(fringe_region1, fringe_region2)
    similarity.append(result)
    print(f'Similarity: {result.flatten()}')
#Plot similarity on alt az
fig, ax = plt.subplots(dpi=200, subplot_kw={'projection':'polar'}, figsize=(10,10))
plot = ax.scatter(np.array(azimuth)*np.pi/180, altitude, c = similarity, s = rotangle, alpha=0.8)
ax.set_title(r'Fringe $\sigma$ with altitude-azimuth', fontsize=20)
ax.set_xlabel(r'Altitude ($\degree$)', labelpad=5, fontsize=20)
ax.set_rgrids((20,40,60,80))
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'Fringing $\sigma$', rotation=270, labelpad=25, fontsize=16)
plt.show()
# ij = np.unravel_index(np.argmax(result), result.shape)
# x, y = ij[::-1]

# fig = plt.figure(figsize=(8, 3))
# ax1 = plt.subplot(1, 3, 1)
# ax2 = plt.subplot(1, 3, 2)
# ax3 = plt.subplot(1, 3, 3, sharex=ax2, sharey=ax2)

# ax1.imshow(fringe_region2, cmap=plt.cm.gray)
# ax1.set_axis_off()
# ax1.set_title('template')

# ax2.imshow(fringe_region1, cmap=plt.cm.gray)
# ax2.set_axis_off()
# ax2.set_title('image')
# # highlight matched region
# hcoin, wcoin = fringe_region2.shape
# rect = plt.Rectangle((x, y), wcoin, hcoin, edgecolor='r', facecolor='none')
# ax2.add_patch(rect)

# ax3.imshow(result)
# ax3.set_axis_off()
# ax3.set_title('`match_template`\nresult')
# # highlight matched region
# ax3.autoscale(False)
# ax3.plot(x, y, 'o', markeredgecolor='r', markerfacecolor='none', markersize=10)

# plt.show()
#%% Image registration to determine sub-pixel offset between images
from skimage.registration import phase_cross_correlation
from skimage.registration._phase_cross_correlation import _upsampled_dft
from scipy.ndimage import fourier_shift

image=fringe_region1
offset_image=fringe_region2

# pixel precision first
shift, error, diffphase = phase_cross_correlation(image, offset_image)

fig = plt.figure(figsize=(8, 3))
ax1 = plt.subplot(1, 3, 1)
ax2 = plt.subplot(1, 3, 2, sharex=ax1, sharey=ax1)
ax3 = plt.subplot(1, 3, 3)

ax1.imshow(image, cmap='gray')
ax1.set_axis_off()
ax1.set_title('Reference image')

ax2.imshow(offset_image.real, cmap='gray')
ax2.set_axis_off()
ax2.set_title('Offset image')

# Show the output of a cross-correlation to show what the algorithm is
# doing behind the scenes
image_product = np.fft.fft2(image) * np.fft.fft2(offset_image).conj()
cc_image = np.fft.fftshift(np.fft.ifft2(image_product))
ax3.imshow(cc_image.real)
ax3.set_axis_off()
ax3.set_title("Cross-correlation")

plt.show()

print(f'Detected pixel offset (y, x): {shift}')

# subpixel precision
shift, error, diffphase = phase_cross_correlation(image, offset_image,
                                                  upsample_factor=100)

fig = plt.figure(figsize=(8, 3))
ax1 = plt.subplot(1, 3, 1)
ax2 = plt.subplot(1, 3, 2, sharex=ax1, sharey=ax1)
ax3 = plt.subplot(1, 3, 3)

ax1.imshow(image, cmap='gray')
ax1.set_axis_off()
ax1.set_title('Reference image')

ax2.imshow(offset_image.real, cmap='gray')
ax2.set_axis_off()
ax2.set_title('Offset image')

# Calculate the upsampled DFT, again to show what the algorithm is doing
# behind the scenes.  Constants correspond to calculated values in routine.
# See source code for details.
cc_image = _upsampled_dft(image_product, 150, 100, (shift*100)+75).conj()
ax3.imshow(cc_image.real)
ax3.set_axis_off()
ax3.set_title("Supersampled XC sub-area")


plt.show()

print(f'Detected subpixel offset (y, x): {shift}')