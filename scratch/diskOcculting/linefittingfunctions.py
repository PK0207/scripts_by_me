#Nicole Arulanantham
#May 2016
#Functions to fit H2 (or other) emission lines in HST-COS spectra

import numpy as np
import math
from scipy.io import readsav
from scipy.optimize import curve_fit
import scipy.ndimage as scind
#import matplotlib.pyplot as plt
import os, sys
import pickle
import platform
import pandas as pd

# #if platform.python_version()[0] == 2: #Pidly is only available in Python 2's - won't work in Python 3, which is what GUI has to be run in
# import pidly
#
# idl = pidly.IDL('/Applications/exelis/idl85/bin/idl') #Call to IDL functions
# idl('.compile /Users/Nicole/Documents/CU/Research/RY_Lup/Scripts/H2_PPD_modeling_code/ccm_unred.pro') #Keri's modification for the PPDs has been edited out here

#Need to get line spread functions (for COS data only - STIS lines are Gaussian shaped)
file_table = readsav("cos_lsf_ltp1.idl")
lsf_list = file_table["lsf"]
lsfchan = file_table["lsfchan"]
lsfpixscale = file_table["lsfpixscale"]
lsfpix = file_table["lsfpix"]
lsfwave = file_table["lsfwave"]

def gaus(wavelengths, amplitude, x0, sigma):
    return amplitude*np.exp(-np.square(wavelengths-x0) / (2.*sigma**2.))

def cos_lsf_keeney(mean_wave, pixel): #Original IDL code from B. Keeney; Don't need to include xarray input, because it hasn't been initialized but will be in function; Also don't need "new" input (for now anyway - we'll always use the files from the later lifetime position)

    if mean_wave < 1800.0:
        if mean_wave > 1450.0:
            chan="g160m"
        else:
            chan="g130m"
    else:
        chan="g225m"

    chan_loc = np.where(lsfchan == chan) #Which index in dictionary represents the desired channel

    chan = chan_loc[0]

    #Pick nearest wavelength point (LSF varies slowly enough with lambda that this should be good enough)
    wave_to_minimize = np.absolute(lsfwave[:, chan] - mean_wave) #Need to find where the wavelength we want is in the file - this way can enter exact wavelength instead of closest integer
    lamind = np.where(wave_to_minimize == np.amin(wave_to_minimize))[0]

    lsf = lsf_list[:, lamind, chan] #I don't think this needs to be transposed, but it does need to be flattened into a 1-D array

    if pixel == True:
        xarray = lsfpix #Pixel values from -100 to +100
    else:
        xarray = mean_wave + (lsfpix*0.001*lsfpixscale[chan])

    return xarray, lsf

def cos_lsf_arulanantham(mean_wave, ltp, pixel):

    #Converted original IDL code from B. Keeney to Python, modified to include LSFs from all three COS lifetime positions

    #Arguments:
    #mean_wave = central wavelength of desired LSF; dtype = float
    #ltp = desired lifetime position (enter "LTP1", "LTP2", "LTP3"); dtype = string
    #pixel = return x-values as pixels or wavelength in Angstroms (enter "True" or "False"); dtype = string

    #LTP2 and LTP3 dictionaries created in LSF_dictionaries.py

    #Don't need to include xarray input (from Keeney version), because it hasn't been initialized but will be in function

    #Also don't need "new" input (from Keeney version) - "new" forces code to use the empirical LSFs from LTP1, which we'll always use anyway

    #Tested against E. Tilton's Python conversion for LTP1 and LTP2 - the LSFs agree!

    ### DETERMINE WHICH COS CHANNEL TO PULL LSF FROM ######################
    if mean_wave < 1800.:
        if mean_wave > 1450.:
            chan=b'g160m'
        else:
            chan= b'g130m'
    else:
        chan= b'g225m'

    #### READ IN LSF FILE, BASED ON DESIRED LIFETIME POSITION #############
    if ltp == "LTP1":
        file_table = readsav("cos_lsf_ltp1.idl")

    elif ltp == "LTP2":
        with open("cos_lsf_ltp2.pickle", "rb") as handle:
            file_table = pickle.load(handle)

    elif ltp == "LTP3":
        with open("cos_lsf_ltp3.pickle", "rb") as handle:
            file_table = pickle.load(handle)

    ##### SORT DICTIONARY INTO ARRAYS BY KEYWORDS ########################
    lsf_list = file_table["lsf"]
    lsfchan = file_table["lsfchan"]
    lsfpixscale = file_table["lsfpixscale"]
    lsfpix = file_table["lsfpix"]
    lsfwave = file_table["lsfwave"]

    #### DETERMINE WHICH INDEX REPRESENTS DESIRED CHANNEL ################
    chan_loc = np.where(lsfchan == chan)[0]
    chan = chan_loc[0]

    #### LTP1 ONLY HAS ONE FILE WITH LSFs ################################
    if ltp == "LTP1":
        #Pick nearest wavelength point (LSF varies slowly enough with lambda that this should be good enough)
        wave_to_minimize = np.absolute(lsfwave[:, chan] - mean_wave) #Need to find where the wavelength we want is in the file
        lamind = np.where(wave_to_minimize == np.amin(wave_to_minimize))[0]

        lsf = lsf_list[:, lamind, chan] #I don't think this needs to be transposed, but it does need to be flattened into a 1-D array

    #### LTP2 AND LTP3 HAVE MULTIPLE FILES, SO DICTIONARY HAS EXTRA WAVELENGTH DIMENSION ##
    else:
        lsffiletopick = file_table["lsffiletopick"]
        file_to_pick = np.absolute(lsffiletopick[chan][0] - mean_wave)
        laminwave = np.argmin(file_to_pick)

        wave_to_minimize = np.absolute(lsfwave[chan][0] - mean_wave)
        lamind = np.where(wave_to_minimize == np.amin(wave_to_minimize))[0]

        lsf = lsf_list[chan][0][laminwave, lamind, :].flatten()
        lsf = lsf / np.sum(lsf) #Normalize LSF

    #### RETURN EITHER ARRAY OF PIXEL VALUES OR ARRAY OF WAVELENGTH VALUES ################
    if pixel == True:
        xarray = lsfpix #Pixel values from -100 to +100

    else:
        xarray = mean_wave + (lsfpix*0.001*lsfpixscale[chan])

    return xarray, lsf

def cos_lsf_arulanantham_ullyses(mean_wave, ltp, pixel):

    #Converted original IDL code from B. Keeney to Python, modified to include LSFs from all three COS lifetime positions

    #Arguments:
    #mean_wave = central wavelength of desired LSF; dtype = float
    #ltp = desired lifetime position (enter "LTP1", "LTP2", "LTP3"); dtype = string
    #pixel = return x-values as pixels or wavelength in Angstroms (enter "True" or "False"); dtype = string

    #LTP2 and LTP3 dictionaries created in LSF_dictionaries.py

    #Don't need to include xarray input (from Keeney version), because it hasn't been initialized but will be in function

    #Also don't need "new" input (from Keeney version) - "new" forces code to use the empirical LSFs from LTP1, which we'll always use anyway

    #Tested against E. Tilton's Python conversion for LTP1 and LTP2 - the LSFs agree!

    ### DETERMINE WHICH COS CHANNEL TO PULL LSF FROM ######################
    if mean_wave < 1800.:
        if mean_wave > 1470.: #ULLYSES data switches to G160M at 1470
            chan=b'g160m'
        else:
            chan= b'g130m'
    else:
        chan= b'g225m'

    #### READ IN LSF FILE, BASED ON DESIRED LIFETIME POSITION #############
    if ltp == "LTP1":
        file_table = readsav("cos_lsf_ltp1.idl")

    elif ltp == "LTP2":
        with open("cos_lsf_ltp2.pickle", "rb") as handle:
            file_table = pickle.load(handle)

    elif ltp == "LTP3":
        with open("cos_lsf_ltp3.pickle", "rb") as handle:
            file_table = pickle.load(handle)

    ##### SORT DICTIONARY INTO ARRAYS BY KEYWORDS ########################
    lsf_list = file_table["lsf"]
    lsfchan = file_table["lsfchan"]
    lsfpixscale = file_table["lsfpixscale"]
    lsfpix = file_table["lsfpix"]
    lsfwave = file_table["lsfwave"]

    #### DETERMINE WHICH INDEX REPRESENTS DESIRED CHANNEL ################
    chan_loc = np.where(lsfchan == chan)[0]
    chan = chan_loc[0]

    #### LTP1 ONLY HAS ONE FILE WITH LSFs ################################
    if ltp == "LTP1":
        #Pick nearest wavelength point (LSF varies slowly enough with lambda that this should be good enough)
        wave_to_minimize = np.absolute(lsfwave[:, chan] - mean_wave) #Need to find where the wavelength we want is in the file
        lamind = np.where(wave_to_minimize == np.amin(wave_to_minimize))[0]

        lsf = lsf_list[:, lamind, chan] #I don't think this needs to be transposed, but it does need to be flattened into a 1-D array

    #### LTP2 AND LTP3 HAVE MULTIPLE FILES, SO DICTIONARY HAS EXTRA WAVELENGTH DIMENSION ##
    else:
        lsffiletopick = file_table["lsffiletopick"]
        file_to_pick = np.absolute(lsffiletopick[chan][0] - mean_wave)
        laminwave = np.argmin(file_to_pick)

        wave_to_minimize = np.absolute(lsfwave[chan][0] - mean_wave)
        lamind = np.where(wave_to_minimize == np.amin(wave_to_minimize))[0]

        lsf = lsf_list[chan][0][laminwave, lamind, :].flatten()
        lsf = lsf / np.sum(lsf) #Normalize LSF

    #### RETURN EITHER ARRAY OF PIXEL VALUES OR ARRAY OF WAVELENGTH VALUES ################
    if pixel == True:
        xarray = lsfpix #Pixel values from -100 to +100

    else:
        xarray = mean_wave + (lsfpix*0.001*lsfpixscale[chan])

    return xarray, lsf

def fcosx_function(wavelengths, nlines, polyorder, ltp, *a_set):

    #IDL version from Kevin; Don't need underscores in front of variable names
    #INPUT VARIABLES:
    #x - wavelength (independant variable) array

    #a_set - array of parameters to generate the model with:
    #a[0+3*i] - amplitude of gaussian i
    #a[1+3*i] - central wavelength of gaussian i
    #a[2+3*i] - sigma of gaussian i
    #a[last polyorder+1 parameters] = poly params

    a_set = list(a_set[0])

    n = np.size(a_set)

    acont = a_set[n-polyorder-1:n] #a_set[0][n-polyorder-1:n] #Estimates for polynomial fit
    a = a_set[0:nlines*3] #a_set[0][0:nlines*3] #Estimates for Gaussian parameters

    widths = []
    for val in np.arange(2, nlines*3, 3):
        widths.append(a[val])

    #model over dataset
    minwave = min(wavelengths)
    maxwave = max(wavelengths)
    pixelspacing = np.median(wavelengths - np.roll(wavelengths, 1))

    x = wavelengths #np.array([min(wavelengths) + (pixelspacing*idx) for idx in range(0, len(wavelengths))]).astype(np.float64)

    #Generate background
    f = np.zeros(np.size(x))
    for i in range(0, polyorder+1):
        f = f + acont[i]*(x - np.mean(x))**i

    #Loop through the rest of the a array, taking 3 parameters at a time and generate Gaussians
    for i in range(0, nlines):
        z2=(x - a[1+3*i]) / a[2+3*i]
        MASK = z2**2
        for j in range(0, np.size(MASK)):
            if MASK[j] >= 1000.0: #Ask Kevin about this - currently putting a big dip in the middle of the Gaussian
                MASK[j] = 1.0
        gauss = a[0+3*i]*np.exp(-np.square(z2) / (2.0*1.0))*1.0
        f = f+gauss

    if ltp == "LTP0":
        f_nolsf = np.interp(wavelengths, x, f)
        return f_nolsf

    else:
        lsfx, lsfy = cos_lsf_arulanantham(a_set[1], ltp, False)
        lsfx = lsfx - a_set[1] #Center line profile at zero
        lsfy = lsfy.flatten() #Shape of lsfx = (201,); shape of lsfy = (201, 1)
        #USER LSF
        #pixel spacing of data
        pix = np.median(x - np.roll(x, 1))

        if len(lsfx) < len(wavelengths): #Need to increase length of LSF to match data
            #Number of final pixels to use based on original LSF
            nfinelsf = (np.ceil(np.amax(np.absolute(lsfx)) / pix)+1)*2 + 1
            #Interpolate to data scale
            xfinelsf = (np.arange(1, nfinelsf+1, 1) - (nfinelsf - 1.0)/2.0)*pix
            yfinelsf = np.interp(xfinelsf, lsfx, lsfy)

            lsf = yfinelsf / np.sum(yfinelsf) #Normalize LSF

        else: #Need to decrease length of LSF (LTP2 and LTP3 have more values, so length of LSF is longer than length of data)
            reduction_factor = 1. - ((float(np.size(lsfy)) - float(np.size(x))) / float(np.size(lsfy)))
            yfinelsf = scind.interpolation.zoom(lsfy, reduction_factor)

            lsf = yfinelsf / np.sum(yfinelsf)

        fa = np.convolve(f, lsf, mode="same") #Convolution is commutative
        f = np.interp(wavelengths, x, fa) #This is the model!

        return f

def get_range(wavelengths, fluxes, errors, min_wave, max_wave):

    """Get subset of spectrum with min_wave < wavelengths < max_wave
    :param wavelengths: array of wavelengths in spectrum
    :param fluxes: array of flux values corresponding to "wavelengths"
    :param errors: array of flux uncertainties correpsonding to "fluxes"
    :param min_wave: float, minimum wavelength in subset of spectrum
    :param max_wave: float, maximum wavelength in subset of spectrum
    :return: subsets of "wavelengths", "fluxes", and "errors"
    """
    subset_indices = (wavelengths >= min_wave) & (wavelengths <= max_wave) & (errors > 1.e-16)

    return wavelengths[subset_indices], fluxes[subset_indices], errors[subset_indices]

def fit_model(wavelengths, fluxes, flux_errors, a, num_lines, poly, ltp, bounds=(-np.inf, np.inf), method="lm"): #Doesn't work with current Python 3 setup (NEED TO FIX THIS), but use this to fit absorption + emission

    nlines = num_lines
    polyorder = poly

    if np.all(np.isinf(bounds)) == True and method == "lm": #Default options
        popt, pcov = curve_fit(lambda wavelengths, *a: fcosx_function(wavelengths, nlines, polyorder, ltp, a), wavelengths, fluxes, p0=a)

    else:
        popt, pcov = curve_fit(lambda wavelengths, *a: fcosx_function(wavelengths, nlines, polyorder, ltp, a), wavelengths, fluxes, p0=a, bounds=bounds, method=method)

    perr = np.sqrt(np.diag(pcov)) #Take diagonal entries of covariance matrix (variance of parameter estimates), then take square root to get one standard deviation errors

    f_model = fcosx_function(wavelengths, nlines, polyorder, ltp, popt)

    return f_model, popt, perr


def stis_lsf_cruzaguirre(wave,grat,slit,norm=False):
    data=pd.read_csv(f'LSF_{grat}_{wave}.txt',skiprows=1,sep='\\s{2,}',engine='python')
    #\\s{2,} treats spaces from 2 - infinity as a separator, this is a regex separator that needs the python engine to work
    lsfp=data['Rel pixel']
    lsfy=data[slit]
    if grat=='G140L':
        disp=0.6 #0.6 A/px for G140L, from handbook
    lsfx = lsfp*disp #relative angstron space
    if norm:
        return lsfx,lsfy/np.sum(lsfy) #return wavelength vs normalized psf
    else:
        return lsfx,lsfy #return wavelength vs raw psf