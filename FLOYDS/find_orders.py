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
import scipy.stats as stats
import os
import time
from datetime import datetime
import sys
import matplotlib.pyplot as plt

from banzai_floyds.orders import Orders, OrderLoader, OrderSolver
from banzai.calibrations import make_master_calibrations
from banzai_floyds import settings
from banzai import dbs
from banzai.utils.stage_utils import run_pipeline_stages
from banzai.logs import set_log_level
from banzai.context import Context
import logging
from banzai_floyds.frames import FLOYDSFrameFactory
from banzai.bias import OverscanSubtractor                                                                
from banzai.trim import Trimmer
from banzai.gain import GainNormalizer
from banzai.uncertainty import PoissonInitializer
from banzai_floyds.wavelengths import WavelengthSolutionLoader
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
            
        frame_factory = FLOYDSFrameFactory()
        skyflat = glob(f'{self.skyflats_path}/*.fz', recursive=True)
        skyflat_hdu = fits.open(skyflat[0])
        skyflat_hdu['SCI'].header['OBSTYPE'] = 'SKYFLAT'
        skyflat_hdu.writeto(skyflat[0], overwrite=True)
        skyflat_hdu.close()
        
        for image_path in skyflat:
            print(image_path)
            cal_image = self.frame_factory.open({'path': image_path}, self.context)
            cal_image.is_master = True
            dbs.save_calibration_info(cal_image.to_db_record(DataProduct(None, filename=os.path.basename(image_path),
                                                                                filepath=os.path.dirname(image_path))),
                                                                                os.environ['DB_ADDRESS'])
        
        overscan_subtract = OverscanSubtractor(self.context)
        trim = Trimmer(self.context)
        gain_norm = GainNormalizer(self.context)
        uncertainty = PoissonInitializer(self.context)
        load_orders = OrderLoader(self.context)
        solve_orders = OrderSolver(self.context)
        
        files = glob(f'{self.lampflats_path}/*.fz', recursive=True)
        files = sorted(files)
        for path in files:
            image = frame_factory.open({'path': path}, self.context)
            print(image.instrument)
            image = overscan_subtract.do_stage(image)
            image = trim.do_stage(image)
            image = gain_norm.do_stage(image)
            print('uncertainty')
            image = uncertainty.do_stage(image)
            print('load orders')
            image = load_orders.do_stage(image)
            print('solve orders')
            image = solve_orders.do_stage(image)
            image.write(self.context)
#%%
pipeline = FLOYDSPipeline()
pipeline.setup_pipeline(processed_path = '/home/pkottapalli/FLOYDS/data/test_data/')
pipeline.run_pipeline(lampflats_path = '/home/pkottapalli/FLOYDS/data/New_AltAz_data/', skyflats_path='/home/pkottapalli/FLOYDS/data/skyflats/')