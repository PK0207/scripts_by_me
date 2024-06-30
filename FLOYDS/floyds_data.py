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
import os
#%% Authenticate
requests.post(
    'https://observe.lco.global/api/api-token-auth/',
    data = {
        'username': 'eng@lco.global',
        'password': 'sbatoo1'
    }
).json()
#%%
import csv
csv_path = 'On_demand_report_2023-03-17T22_37_35.945Z_5087d390-c514-11ed-a268-bd5dde369fb5.csv'
filenames = []
with open(csv_path, newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        filename = row[0]
        filenames.append(filename)

#%%
# def _get_observations_offset(query_params: dict):
#     observations_endpoint = urljoin('https://archive-api.lco.global', 'frames')
#     try:
#         url = requests.Request('GET', observations_endpoint, params=query_params).prepare().url
#         response = requests.get(url, headers={'Authorization': 'Token 9b605ef2229d317ba1b031b54f2a0115aec69b9f'})
#         response.raise_for_status()
#     except requests.HTTPError:
#         logging.warning("observations endpoint {} failed with status code {}: {}".format(response.url,
#                                                                                          response.status_code,
#                                                                                          response.text))
#         return []
#     return response.json()
# def get_observations(query_params):
#     query_params['offset'] = 0
#     response = _get_observations_offset(query_params)
#     observations = response['results']
#     while response['count'] > len(observations):
#         query_params['offset'] = len(observations)
#         response = _get_observations_offset(query_params)
#         observations.extend(response['results'])

#     return observations

#query_params={'offset':0}#, 'instrument_id':'en06', 'site_id':'ogg', 'start':'2022-01-01', 'end':'2023-01-01', 'configuration_type':'LAMPFLAT', 'public':True}
#Get all headers from a year
altitude = []
azimuth = []
rotangle = []
frameid = []
for f in filenames[1:]:
    filename = os.path.splitext(os.path.splitext(f)[0])[0]
    print(filename)
    archive_record = requests.get(f'https://archive-api.lco.global/frames/?site_id=ogg&basename={filename}&start=2022-01-01 00%3A00&end=2023-01-01 00%3A00&public=true&limit=1000', headers={'Authorization': 'Token efc8c22ed48db4962008085fc4af4bfa5354fd7d'}).json()
    #while archive_record['count'] > len(altitude):
    if len(archive_record['results']) > 0:
        rec = archive_record['results'][0]
        #query_params['offset'] = len(altitude)
        #response = _get_observations_offset(query_params)
        frame_id = rec["id"]
        frameid.append(frame_id)
        header = requests.get(f'https://archive-api.lco.global/frames/{frame_id}/headers/', headers={'Authorization': 'Token efc8c22ed48db4962008085fc4af4bfa5354fd7d'}).json()
        altitude.append(header['data']['ALTITUDE'])
        azimuth.append(header['data']['AZIMUTH'])
        rotangle.append(header['data']['ROTANGLE'])
    else:
        print('could not find ' + filename)
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

doubles_path = 'All_AltAz_data'
def download_file(records):
    for frame_id in records:
        print(frame_id)
        archive_record = requests.get(f'https://archive-api.lco.global/frames/?instrument_id=en06&site_id=ogg&start=2022-01-01&end=2023-01-01&configuration_type=LAMPFLAT&id{frame_id}&public=true&limit=1000', headers={'Authorization': 'Token efc8c22ed48db4962008085fc4af4bfa5354fd7d'}).json()['results']
        for rec in archive_record:
            #Give path to write files to
            with open(f'{doubles_path}/{rec["filename"]}', 'wb') as f:
                f.write(requests.get(rec['url']).content)

N_PROCESSES = multiprocessing.cpu_count()
with multiprocessing.Pool(N_PROCESSES) as pool:
    pool.map(download_file, unique_coords)
    pool.close()
    pool.join()

#%%
download_file(unique_coords)
