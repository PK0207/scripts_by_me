#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 11:02:53 2023

@author: pkottapalli
"""

"""
Get all unique alt az data. Including red and blue lampflats,
red and blue science frames
all rectified, flat corrected,
and also raw
"""
#%% Imports
import requests
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from urllib.parse import urljoin
#%% Authenticate
requests.post(
    'https://observe.lco.global/api/api-token-auth/',
    data = {
        'username': 'pkottapalli@lco.global',
        'password': '&2m9X<n2'
    }
).json()

#%%
def _get_observations_offset(query_params: dict):
    observations_endpoint = urljoin('https://archive-api.lco.global', 'frames')
    try:
        url = requests.Request('GET', observations_endpoint, params=query_params).prepare().url
        response = requests.get(url, headers={'Authorization': 'Token 9b605ef2229d317ba1b031b54f2a0115aec69b9f'})
        response.raise_for_status()
    except requests.HTTPError:
        logging.warning("observations endpoint {} failed with status code {}: {}".format(response.url,
                                                                                         response.status_code,
                                                                                         response.text))
        return []
    return response.json()
def get_observations(query_params):
    query_params['offset'] = 0
    response = _get_observations_offset(query_params)
    observations = response['results']
    while response['count'] > len(observations):
        query_params['offset'] = len(observations)
        response = _get_observations_offset(query_params)
        observations.extend(response['results'])

    return observations
query_params={'offset':0}#, 'instrument_id':'en06', 'site_id':'ogg', 'start':'2022-01-01', 'end':'2023-01-01', 'configuration_type':'LAMPFLAT', 'public':True}
#Get all headers from a year
altitude = []
azimuth = []
rotangle = []
frameid = []
archive_record = requests.get('https://archive-api.lco.global/frames/?instrument_id=en06&site_id=ogg&start=2022-01-01&end=2023-01-01&configuration_type=LAMPFLAT&public=True&limit=1000', headers={'Authorization': 'Token 9b605ef2229d317ba1b031b54f2a0115aec69b9f'}).json()
#while archive_record['count'] > len(altitude):
for rec in archive_record['results']:    
    query_params['offset'] = len(altitude)
    #response = _get_observations_offset(query_params)
    frame_id = rec["id"]
    frameid.append(frame_id)
    header = requests.get(f'https://archive-api.lco.global/frames/{frame_id}/headers/', headers={'Authorization': 'Token 9b605ef2229d317ba1b031b54f2a0115aec69b9f'}).json()
    altitude.append(header['data']['ALTITUDE'])
    azimuth.append(header['data']['AZIMUTH'])
    rotangle.append(header['data']['ROTANGLE'])
df = pd.DataFrame({'Altitude':altitude, 'Azimuth':azimuth, 'Rotangle': rotangle}, index=frameid)
#%% Plot the locations of each coordinate
plt.figure(dpi=200)
plt.scatter(altitude, azimuth, s=2, alpha = 0.3)
plt.xlabel('Altitude')
plt.ylabel('Azimuth')
plt.show()
#%% Find observations that were close to each other in the alt-rot space
#Take a circular approximation for the distance between points in alt-rot space
unique_coords = []
for i in range(len(altitude)):
    for k in range(len(altitude)):
        dist = (altitude[i]-altitude[k])**2 + (rotangle[i]-rotangle[k])**2
    if dist > 1:
        unique_coords.append(df.index[i])
#%% download the frames that are unique
import multiprocessing

doubles_path = 'New_AltAz_data'
def download_file(records):
    for frame_id in records:
        archive_record = requests.get(f'https://archive-api.lco.global/frames/?instrument_id=en06&site_id=ogg&start=2022-01-01&end=2023-01-01&configuration_type=LAMPFLAT&id{frame_id}&public=True&limit=1000', headers={'Authorization': 'Token 9b605ef2229d317ba1b031b54f2a0115aec69b9f'}).json()['results']
        for rec in archive_record:
            #Give path to write files to
            with open(f'{doubles_path}/{rec["filename"]}', 'wb') as f:
                f.write(requests.get(rec['url']).content)
            
N_PROCESSES = multiprocessing.cpu_count()
with multiprocessing.Pool(N_PROCESSES) as pool:
    pool.map(download_file, unique_coords)
    pool.close()
    pool.join()