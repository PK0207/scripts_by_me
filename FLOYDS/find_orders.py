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
            
        #skyflat = glob(os.path.join(self.skyflats_path, '*.fz'), recursive=True)
        files = glob(os.path.join(self.lampflats_path, '*.fz'), recursive=True)
        files = sorted(files)
        
        overscan_subtract = OverscanSubtractor(self.context)
        trim = Trimmer(self.context)
        gain_norm = GainNormalizer(self.context)
        uncertainty = PoissonInitializer(self.context)
        load_orders = OrderLoader(self.context)
        solve_orders = OrderSolver(self.context)
        tweak_orders = OrderTweaker(self.context)
        
        #for image_path in skyflat:
            # print(image_path)
            # cal_image = self.frame_factory.open({'path': image_path}, self.context)
            # cal_image = overscan_subtract.do_stage(cal_image)
            # cal_image = trim.do_stage(cal_image)
            # cal_image = gain_norm.do_stage(cal_image)
            # cal_image = uncertainty.do_stage(cal_image)
            # cal_image = solve_orders.do_stage(cal_image)
            # cal_image.write(self.context)
            # dbs.save_calibration_info(cal_image.to_db_record(DataProduct(None, filename=os.path.basename(image_path),
            #                                                                     filepath=os.path.dirname(image_path))),
            #                                                                     os.environ['DB_ADDRESS'])
        
        skyflat_path = os.path.join('lampflat_tweaked_data','ogg','en06','20190329','processed', 'ogg2m001-en06-20190329-0018-x92.fits.fz')
        #skyflat_path = os.path.join('all_tweaked_data','ogg','en06','20220105','processed', 'ogg2m001-en06-20220105-0003-w92.fits.fz')
        
        skyflat_image = self.frame_factory.open({'path': skyflat_path}, self.context)
        skyflat_image.is_master = True
        skyflat_image.meta['OBSTYPE']='SKYFLAT'
        dbs.save_calibration_info(skyflat_image.to_db_record(DataProduct(None, filename=os.path.basename(skyflat_path),
                                                                            filepath=os.path.dirname(skyflat_path))),
                                                                            os.environ['DB_ADDRESS'])
        
        # skyflat_image = self.frame_factory.open({'path': 'test_data/ogg/en06/20220105/processed/ogg2m001-en06-20220105-0003-w92.fits.fz'}, self.context)
        # skyflat_image.is_master = True
        # skyflat_image.meta['OBSTYPE']='SKYFLAT'
        # skyflat_image.grouping_criteria = ['filter']
        # dbs.save_calibration_info(skyflat_image.to_db_record(DataProduct(None, filename=os.path.basename('test_data/ogg/en06/20220105/processed/ogg2m001-en06-20220105-0003-w92.fits.fz'),
        #                                                                     filepath=os.path.dirname('test_data/ogg/en06/20220105/processed/ogg2m001-en06-20220105-0003-w92.fits.fz'))),
        #                                                                     os.environ['DB_ADDRESS'])
        files = glob(os.path.join(self.lampflats_path, '*.fz'), recursive=True)
        files = sorted(files)
        
        for path in files:
            print(path)
            image = self.frame_factory.open({'path': path}, self.context)
            print(image.instrument)
            image = overscan_subtract.do_stage(image)
            image = trim.do_stage(image)
            image = gain_norm.do_stage(image)
            print('uncertainty')
            image = uncertainty.do_stage(image)
            print('load orders')
            image = load_orders.do_stage(image) #Load order estimates from calibration image
            # print('solve orders')
            # image = solve_orders.do_stage(image)
            print('tweak orders')
            image = tweak_orders.do_stage(image)
            image.write(self.context)
            print('DONE')
            
#%%
pipeline = FLOYDSPipeline()
pipeline.setup_pipeline(processed_path = 'yrot_tweaked_data')
pipeline.run_pipeline(lampflats_path = 'All_AltAz_data', skyflats_path='skyflats')
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
skyflat = fits.open('/home/pkottapalli/FLOYDS/data/skyflats/ogg2m001-en06-20190329-0018-x00.fits.fz')

sky_alt = skyflat['SCI'].header['ALTITUDE']
sky_az = skyflat['SCI'].header['AZIMUTH']
sky_rot = skyflat['SCI'].header['ROTANGLE']
alt_az_dist = [(altitude[0]-alt)**2 + (azimuth[0]-az)**2 for alt, az in zip(altitude, azimuth)]
plt.figure(dpi=200)
plt.scatter(wmstemp, yshift, s=7)
#plt.xlim(xmin=-50, xmax=2500)
plt.ylabel('Y shift (pixel)')
#plt.xlabel('Distance in alt-rot space')
plt.xlabel('WMSTEMP')
plt.title('Y shift vs. environment temperature')
plt.show()
#%% Solve just skyflat
pipeline = FLOYDSPipeline()
pipeline.setup_pipeline(processed_path = 'lampflat_tweaked_data')

overscan_subtract = OverscanSubtractor(pipeline.context)
trim = Trimmer(pipeline.context)
gain_norm = GainNormalizer(pipeline.context)
uncertainty = PoissonInitializer(pipeline.context)
load_orders = OrderLoader(pipeline.context)
solve_orders = OrderSolver(pipeline.context)
tweak_orders = OrderTweaker(pipeline.context)
skyflat = glob(os.path.join('skyflats', '*.fz'), recursive=True)
for image_path in skyflat:
    print(image_path)
    cal_image = pipeline.frame_factory.open({'path': image_path}, pipeline.context)
    cal_image = overscan_subtract.do_stage(cal_image)
    cal_image = trim.do_stage(cal_image)
    cal_image = gain_norm.do_stage(cal_image)
    cal_image = uncertainty.do_stage(cal_image)
    cal_image = solve_orders.do_stage(cal_image)
    cal_image.write(pipeline.context)
#%%Compare data to region of fit
from numpy.polynomial.legendre import Legendre
from banzai_floyds.orders import order_region

def normalize_2d(matrix):
    norm = np.linalg.norm(matrix)
    matrix = matrix/norm  # normalized matrix
    return matrix

path = 'lampflat_tweaked_data/ogg/en06/20190329/processed/ogg2m001-en06-20190329-0018-x92.fits.fz'
#path = files[0]
image = pipeline.frame_factory.open({'path': path}, pipeline.context)
order_estimates = [(coeff, height, domain)
                   for coeff, height, domain in
                   zip(image.orders.coeffs, image.orders.order_heights, image.orders.domains)]

for i, (coeff, height, domain) in enumerate(order_estimates):

    x2d, y2d = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
    #region = np.logical_and(x2d <= domain[1], x2d >= domain[0])
    center_model = Legendre(coeff, domain=domain)
    region = order_region(height, center_model, image.shape)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, dpi=200)
    #im1 = ax1.imshow(image.data[0:200, 0:1700], origin='lower')
    im1 = ax1.imshow(image.data, origin='lower')
    ax1.set_title('Skyflat Image data', fontsize=7)
    #ax1.axis('off')
    plt.colorbar(im1)
    #im2 = ax2.imshow(region[0:200, 0:1700], origin='lower')
    im2 = ax2.imshow(region, origin='lower')
    ax2.set_title('Order Region', fontsize=7)
    #ax2.axis('off')
    plt.colorbar(im2)
    diffim = normalize_2d(image.data[0:200, 0:1700]) - region[0:200, 0:1700]
    median = np.median(diffim)
    std = np.std(diffim)
    im3 = ax3.imshow(diffim, origin='lower', vmin = median-3*std, vmax = median+3*std)
    ax3.set_title('Image data-order region', fontsize=7)
    #ax3.axis('off')
    plt.colorbar(im3)
    fig.show()