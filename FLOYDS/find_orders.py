#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 14:53:54 2023

@author: pkottapalli
"""

from glob import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime
import sys

from banzai_floyds.orders import OrderLoader, OrderTweaker, OrderSolver
from banzai_floyds import settings
from banzai import dbs
from banzai.logs import set_log_level
from banzai.context import Context
import logging
from banzai_floyds.frames import FLOYDSFrameFactory
from banzai.bias import OverscanSubtractor                                                                
from banzai.trim import Trimmer
from banzai.gain import GainNormalizer
from banzai.uncertainty import PoissonInitializer
from banzai.data import DataProduct
#%%
class FLOYDSPipeline():
    def __init__(self):
        return

    def setup_pipeline(self, processed_path, db_path=None):
        print('Please do not run this more than one time. Database is set up in your current directory if not provided')     
        set_log_level('DEBUG')
        logger = logging.getLogger('banzai')
        
        self.db_path = db_path
        if self.db_path:
            if self.db_path[-1] == '/':
                self.db_path = self.db_path[:-1]
                os.environ['DB_ADDRESS'] = f'sqlite:///{db_path}/test.db'
        else:
            os.environ['DB_ADDRESS'] = 'sqlite:///test.db'
        os.environ['CONFIGDB_URL'] = 'http://configdb.lco.gtn/sites'
        os.environ['OPENTSDB_PYTHON_METRICS_TEST_MODE'] = 'True'
        
        settings.processed_path = os.path.join(os.getcwd(), 'test_data')
        settings.fpack=True
        settings.db_address = os.environ['DB_ADDRESS']
        settings.reduction_level = 92 #This doesn't matter for sandboxing
        
        # set up the context object.
        import banzai.main
        context = banzai.main.parse_args(settings, parse_system_args=False)
        context = vars(context)
        context['no_bpm'] = True
        context['processed_path'] = processed_path
        context['post_to_archive'] = False
        context['no_file_cache'] = False
        self.context = Context(context)
        
        # initialize the DB with some instruments from ConfigDB
        
        os.system(f'banzai_create_db --db-address={os.environ["DB_ADDRESS"]}')
        os.system(f'banzai_update_db --db-address={os.environ["DB_ADDRESS"]} --configdb-address={os.environ["CONFIGDB_URL"]}')
        
        # wow, after all that you can actually open an image!
        
        self.frame_factory = FLOYDSFrameFactory()
    
    def run_pipeline(self, lampflats_path=None, skyflats_path=None):
        #Load flats in
        try:
            self.lampflats_path
        except AttributeError:
            if not lampflats_path:
                self.lampflats_path = input('Please provide a path to your lampflats: ')
            else:
                self.lampflats_path = lampflats_path
        
        try:
            self.skyflats_path
        except AttributeError:
            if not lampflats_path:
                self.skyflats_path = input('Please provide a path to your lampflats: ')
            else:
                self.skyflats_path = skyflats_path
        
        try:
            self.context
        except AttributeError:
            print('Please run the .setup_pipeline method before running the pipeline')
            sys.exit(-1)
            
        skyflat = glob(os.path.join(self.skyflats_path, '*.fz'), recursive=True)
        
        overscan_subtract = OverscanSubtractor(self.context)
        trim = Trimmer(self.context)
        gain_norm = GainNormalizer(self.context)
        uncertainty = PoissonInitializer(self.context)
        load_orders = OrderLoader(self.context)
        solve_orders = OrderSolver(self.context)
        tweak_orders = OrderTweaker(self.context)
        
        # for image_path in skyflat:
        #     print(image_path)
        #     cal_image = self.frame_factory.open({'path': image_path}, self.context)
        #     cal_image.is_master = True
        #     cal_image = overscan_subtract.do_stage(cal_image)
        #     cal_image = trim.do_stage(cal_image)
        #     cal_image = gain_norm.do_stage(cal_image)
        #     cal_image = uncertainty.do_stage(cal_image)
        #     cal_image = solve_orders.do_stage(cal_image)
        #     cal_image.write(self.context)
        skyflat_path = os.path.join('ogg2m001-en06-20190329-0018-x92.fits.fz')
        skyflat_image = self.frame_factory.open({'path': skyflat_path}, self.context)
        skyflat_image.is_master = True
        dbs.save_calibration_info(skyflat_image.to_db_record(DataProduct(None, filename=os.path.basename(skyflat_path),
                                                                            filepath=os.path.dirname(skyflat_path))),
                                                                            os.environ['DB_ADDRESS'])
        files = glob(os.path.join(self.lampflats_path, '*.fz'), recursive=True)
        files = sorted(files)
        for path in files[721:]:
            image = self.frame_factory.open({'path': path}, self.context)
            print(image.instrument)
            image = overscan_subtract.do_stage(image)
            image = trim.do_stage(image)
            image = gain_norm.do_stage(image)
            print('uncertainty')
            image = uncertainty.do_stage(image)
            print('load orders')
            image = load_orders.do_stage(image) #Load order estimates from calibration image
            #print('solve orders')
            #image = solve_orders.do_stage(image)
            print('tweak orders')
            image = tweak_orders.do_stage(image)
            image.write(self.context)
            print('DONE')
            
#%%
pipeline = FLOYDSPipeline()
pipeline.setup_pipeline(processed_path = 'all_tweaked_data')
pipeline.run_pipeline(lampflats_path = 'All_AltAz_data', skyflats_path='skyflats')
#%% Show Xshift yshift and rotation with alt az
files = glob('all_tweaked_data/ogg/en06/2022*/processed/*.fz', recursive=True)

headers = [fits.open(f)['SCI'].header for f in files]

aperwid = np.array([header['APERWID'] for header in headers])
exptime = np.array([header['EXPTIME'] for header in headers])
aperwid_2 = np.where(aperwid == 2)
exp_80 = np.where([np.isclose(t, 80, atol = 1) for t in exptime])
exp_40 = np.where([np.isclose(t, 40, atol = 1) for t in exptime])
exp80_and_aperwid2 = np.intersect1d(aperwid_2, exp_80)
exp40_and_aperwid2 = np.intersect1d(aperwid_2, exp_40)
use_headers = [headers[i] for i in exp40_and_aperwid2]
times = np.array([datetime.strptime(header['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f') for header in use_headers])
rotangle = np.array([header['ROTANGLE'] for header in use_headers])
altitude = np.array([header['ALTITUDE'] for header in use_headers])
azimuth = np.array([header['AZIMUTH'] for header in use_headers])
xshift = np.array([header['ORDXSHFT'] for header in use_headers])
yshift = np.array([header['ORDYSHFT'] for header in use_headers])
rotation = np.array([header['ORDROT'] for header in use_headers])
ccdtemp = np.array([header['CCDATEMP'] for header in use_headers])
#%%
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
#%%
import matplotlib.dates as mdates
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, dpi=200, sharex=True, figsize=(21,5))
ax1.scatter(times, xshift, s=9)
#ax1.set_ylim(ymin=200, ymax=750)
locator = mdates.AutoDateLocator()
formatter = mdates.ConciseDateFormatter(locator)
ax1.xaxis.set_major_locator(locator)
ax1.set_xlabel('Time')
ax1.set_ylabel('XSHIFT')

ax2.scatter(times, yshift, s=9)
#ax2.set_ylim(ymin=-200)
ax2.xaxis.set_major_locator(locator)
ax2.set_xlabel('Time')
ax2.set_ylabel('YSHIFT')

rot_plot = ax3.scatter(times, rotation, s=9)
#ax3.set_ylim(ymin=-10)
ax3.xaxis.set_major_locator(locator)
ax3.set_xlabel('Time')
ax3.set_ylabel('ROTATION')
#cbar = plt.colorbar(rot_plot, ax=ax3)
#cbar.ax.set_ylabel(r'Exposure time', rotation=270, labelpad=25)
fig.show()
#%%
fig, ax1 = plt.subplots(dpi=200)
ax1.scatter(times, xshift, s=5, label='xshift')
ax1.scatter(times, yshift, s=5, label='yshift')
ax1.hlines(np.mean(yshift), times[0], times[-1], 'r', '--', label=f'Mean y-shift = {np.mean(yshift)}')
ax1.hlines(np.mean(xshift), times[0], times[-1], 'r', '--', label=f'Mean x-shift = {np.mean(xshift)}')
#ax1.set_ylim(ymin=200, ymax=750)
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