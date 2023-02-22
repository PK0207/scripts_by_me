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
files = sorted(glob('flat*red_2.0*.fits', recursive=True))
images = [fits.open(f)[0].data for f in files]
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

header_values = ['DATE-OBS', 'AZIMUTH', 'ALTITUDE']
headers = [fits.open(f)[0].header for f in files]
sorting_idx = np.argsort([header['UTSTART'] for header in headers])

env_parameters = []
for header in headers:
    env_parameters.append([header[key] for key in header_values])

env_parameters = np.array(env_parameters)

#plt.style.use('seaborn')
fig, axes = plt.subplots(1, 3, sharex=True, figsize=(35,10))
fig.tight_layout(pad = 12)
times = np.sort([header['DATE-OBS'] for header in headers])
t0 = datetime.strptime(times[0], '%Y-%m-%dT%z')
t1 = datetime.strptime(times[-1], '%Y%m%d')
for ax, i in zip(axes.flatten(), range(len(env_parameters[0]))):
    ax.scatter(times, env_parameters[:,i][sorting_idx], s=15)
    ax.set_ylabel(header_values[i], fontsize=20)
    ax.set_xlabel('Obs time', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=14, size=10)
    locator = mdates.DateLocator()
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_minor_locator(locator)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

plt.show()
#%%
