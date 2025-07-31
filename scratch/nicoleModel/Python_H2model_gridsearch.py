"""Script to do MCMC resampling of Python H2 models
Written by N. Arulanantham
Last Updated: November 2024
"""
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import math
import scipy.special as special
import time
import sys
import pickle
from scipy.io import readsav
import matplotlib.pyplot as plt
import gc
import os
from decimal import Decimal
from itertools import product
from scipy.optimize import minimize
from tqdm import tqdm
#Import custom modules
from ccm_unred import which_spec_range, ccm_unred
from linefittingfunctions import cos_lsf_arulanantham
from OOP_HoadleyH2model import H2Model

#Global constants
#Gravitational constant [cm^3 g^-1 s^-2]
G = 6.67259e-8
#Speed of light [cm/s]
c = 3.e10
#Planck constant [erg s]
h = 6.626e-27
#Boltzmann constant [erg/K]
kB = 1.38e-16
#Convert [AU -> cm]
AUtoCM = 1.496e13
#Convert [cm -> km]
CMtoKM = 1.e-5
#Convert [parsecs -> astronomical units]
PCtoAU = 206265.
#Convert [M_solar -> grams]
MSUNtoG = 1.99e33
#Convert [Angstroms -> cm]
ANGtoCM = 1.e-8
#Extinction dependent on line-of-sight
Rv = 3.1

np.seterr(divide='ignore')

def find_wave_loc(wave, ref_wave_array):
    """Match wavelengths in data arrays to wavelengths in spreadsheet with measured line properties
    :param wave: wavelength to find in spreadsheet (float)
    :param ref_wave_array: list of wavelengths in spreadsheet (array of floats)
    :return: wave_loc (float)
    """
    ref_diff = np.absolute(ref_wave_array - wave)
    wave_loc = np.argmin(ref_diff)
    return wave_loc

def get_prog_data(lines_to_fit, ref_wave_array, velocity_array, flux_array, fluxerr_array):
    """Concatenate emission lines for model fitting
    :param lines_to_fit: Wavelengths of emission lines (usually for a single progression)
    :param ref_wave_array: List of wavelengths, in order that the data is stored in Keri-style arrays
    :param velocity_array: Velocity space of all emission lines
    :param flux_array: Fluxes of all emission lines
    :param fluxerr_array: Flux uncertainties for all emission lines
    """
    #Get indices of emission lines to concatenate in data arrays
    prog_wave_keys = [str(round(wave, 2))[0:7] for wave in lines_to_fit]
    #Start concatenating!
    prog_velocity = np.concatenate([velocity_array[key] for key in prog_wave_keys])
    prog_flux = np.concatenate([flux_array[key] for key in prog_wave_keys])
    prog_fluxerr = np.concatenate([fluxerr_array[key] for key in prog_wave_keys])

    return prog_velocity, prog_flux, prog_fluxerr

def grid_from_radius(targ_d, param_dict):
    """Builds grids to describe disk based on stellar distance and disk inclination
    """
    rinout = 10.**(0.02 * np.arange(0, 201, 1) - 1.4)
    r_grid = rinout / targ_d
    phi_grid = np.arange(0., 360., 2.) * math.pi / 180.

    xgrid = np.zeros((len(r_grid), len(phi_grid)))
    ygrid = np.zeros((len(r_grid), len(phi_grid)))

    tan_phi = np.tan(phi_grid)

    for r_idx, r in enumerate(r_grid):
        for az_idx, az in enumerate(tan_phi):
            if az >= 90. * math.pi / 180. and az <= 270. * math.pi / 180.:
                x_dummy = -1. * math.sqrt(r ** 2. / (1. + az))
            else:
                x_dummy = math.sqrt(r ** 2. / (1. + az ** 2.))
            xgrid[r_idx, az_idx] = x_dummy / math.cos(param_dict["Inclination"])

            y_dummy = math.sqrt(r ** 2. / (1. + (1. / az ** 2.)))
            if phi_grid[az_idx] >= math.pi:
                ygrid[r_idx, az_idx] = -1. * y_dummy
            else:
                ygrid[r_idx, az_idx] = y_dummy

    r_H2 = np.zeros((len(r_grid), len(phi_grid)))
    phi_H2 = np.zeros((len(r_grid), len(phi_grid)))
    for r_idx, r in enumerate(r_grid):
        for az_idx, az in enumerate(tan_phi):
            r_dummy = math.sqrt(xgrid[r_idx, az_idx] ** 2. + ygrid[r_idx, az_idx] ** 2.) * targ_d
            r_H2[r_idx, az_idx] = r_dummy

            phi_dummy = math.atan(ygrid[r_idx, az_idx] / xgrid[r_idx, az_idx])
            if ygrid[r_idx, az_idx] <= 0.:
                if xgrid[r_idx, az_idx] <= 0.:
                    phi_H2[r_idx, az_idx] = phi_dummy + math.pi
                else:
                    phi_H2[r_idx, az_idx] = phi_dummy + 2.*math.pi
            else:
                if xgrid[r_idx, az_idx] <= 0.:
                    phi_H2[r_idx, az_idx] = phi_dummy + math.pi
                else:
                    phi_H2[r_idx, az_idx] = phi_dummy
    return {"xgrid": xgrid,
            "ygrid": ygrid,
            "rgrid": (r_grid * targ_d) - 0.03,
            "phigrid": phi_grid * 180. / math.pi,
            "rH2": r_H2,
            "phiH2": phi_H2 * 180. / math.pi}

def radial_velocity(r_grid, targ_M):
    """Done checking against IDL version!
    """
    v_k = np.sqrt(G * targ_M * MSUNtoG / (r_grid * AUtoCM)) * CMtoKM

    return v_k

def vrad_to_vobs(r_grid, phi_grid):
    """Done checking against IDL version!
    """
    v_k = radial_velocity(r_grid, targ_M)
    vobs = [v_k * math.sin(param_dict["Inclination"]) * math.sin(phi * math.pi / 180.) for phi in phi_grid]
    vobs = np.array(vobs).astype(np.float64)
    vobs = np.swapaxes(vobs, 0, 1)
    return vobs

def collapsed_line_profile(lineprofs, param_dict, inds):

    binmin = min(v_obs_bin)
    binmax = max(v_obs_bin)
    total_line_profile = np.zeros((len(v_obs_bin), len(lineprofs[:, 0, 0])), dtype=np.float64)
    total_line_profile.fill(1.e-75)
    for prog_idx in range(0, np.shape(lineprofs)[0]):
        flat_fluxes = lineprofs[prog_idx, :, :].flatten()
        test_df = pd.DataFrame({"Velocities": inds,
                                "Fluxes": flat_fluxes})

        g = test_df.groupby(["Velocities"])

        profile = np.array(g.sum()).astype(np.float64).flatten()

        pads = (600 - len(profile))

        if pads % 2 == 0:
            pad_left = pads / 2
            pad_right = pads / 2
        else:
            pad_left = math.floor(pads / 2)
            pad_right = math.ceil(pads / 2)

        final_profile = np.pad(profile, (pad_left, pad_right), 'edge')

        total_line_profile[:, prog_idx] = final_profile
        total_line_profile[300, prog_idx] = total_line_profile[301, prog_idx]

    return total_line_profile

def get_H2model(theta, theta_keys, param_dict, inds):
    """Concatenate model emission lines into single array for fitting to data
    :param params: List of model parameters (floats)
    :return: full_interp_model (array), array of model emission lines
    """
    start_time = time.time()
    for key_idx, key in enumerate(theta_keys):
        param_dict[key] = theta[key_idx]

    target = H2Model(targ_M, targ_AV, targ_d, param_dict)
    EBV = -1. * targ_AV / Rv

    #This is the optical depth - you can plot it if you want
    tau = target.H2_number_density()

    lineprof = target.total_intensity()

    rvals = target.grid_dict["rgrid"]
    Tprof = target.radial_T(rvals)

    #Uncomment this section to plot the temperature profile
    # plt.clf()
    # plt.close()
    # plt.xscale("log")
    # plt.title("Radial Temperature Distribution: q = "+str(param_dict["q"]))
    # plt.plot(rvals,
    #          Tprof,
    #          color="k", lw=2.)
    # plt.xlim(0.01, 15.)
    # plt.gca().minorticks_on()
    # plt.show(block=True)

    intensity_binned = collapsed_line_profile(lineprof, param_dict, inds)
    v_obs_bin = np.arange(-300., 300., 1.)

    intensity = np.zeros((len(np.arange(-300., 300., 1.)), len(target.lambda_props["Wavelength"])))

    for wave_idx, wave in enumerate(target.lambda_props["Wavelength"]):
        prog_idx = np.argmin(np.absolute(target.lambda_props["Jl"][wave_idx] - target.J_initial))
        intens_unred = intensity_binned[:, prog_idx] * target.lambda_props["Bul"][wave_idx]

        intens_red = ccm_unred(np.zeros(len(intens_unred)) + wave, intens_unred, EBV)
        intens_red[np.argwhere(np.isnan(intens_red) == True).flatten()] = 1.e-30

        #Smooth/convolve with LSF
        lsfx, lsfy = cos_lsf_arulanantham(wave, "LTP1", False)

        lsfy_norm = lsfy.flatten() / np.sum(lsfy.flatten()) #Normalize LSF

        intensity[:, wave_idx] = np.convolve(intens_red, lsfy_norm, mode="same")

    full_interp_model = []
    for prog_idx, prog in enumerate(lines_to_fit_dict.keys()):
        prog_wave_idxs = [find_wave_loc(wave, model_refwaves) for wave in lines_to_fit_dict[prog]]
        for wave_idx in prog_wave_idxs:
            smoothflux = intensity[:, wave_idx]
            key = str(round(model_refwaves[wave_idx], 2))[0:7]
            interp_smoothflux = griddata(v_obs_bin, smoothflux, target_velocity[key], method = "nearest")
            full_interp_model.append(interp_smoothflux)
    end_time = time.time()

    return np.concatenate(full_interp_model)

def lnlike(theta, velocity, flux, flux_err, theta_keys, param_dict, inds):
    """Log-likelihood function for H2 models (assuming Gaussian distribution, for lack of a better guess)
    :param theta: List of model parameters
    :param velocity: velocity space of all emission lines to be modeled
    :param flux: fluxes for all emission lines to be modeled
    :param flux_err: flux uncertainties for all emission lines to be modeled
    :return: value of log-likelihood function for the given model parameters
    """
    full_interp_model = get_H2model(theta, theta_keys, param_dict, inds)
    str_list = [str(val) for val in theta]
    save_fig_str = "_".join(str_list)

    loglike = np.sum(np.square(flux - full_interp_model) / flux_err) / (len(theta) - 1)
    #print (loglike, theta)

    return loglike

def lnprior(theta, theta_keys, prior_dict):
    """Define bounds of parameter space to probe
    :param theta: List of model parameters
    :return: 0 if all the parameters are in bounds (nothing added to log-likelihood function), -Infinity otherwise (invalidate log-likelihood function)
    """
    param_dict = {key:theta[key_idx] for key_idx, key in enumerate(theta_keys)}
    in_bounds = 0
    for param_idx, param in enumerate(param_dict.keys()):
        prior_bounds = prior_dict[param]
        param_val = param_dict[param]
        if param_val < prior_bounds[0] or param_val > prior_bounds[1]:
            in_bounds += 1
    if in_bounds == 0:
        return 0.
    else:
        return -np.inf

def lnprob(theta, velocity, flux, flux_err, theta_keys, param_dict, prior_dict, inds):
    """Calculate Bayesian posterior probability that these parameters are correct, given the data
    :param theta: List of model parameters
    :param velocity: velocity space of all emission lines to be modeled
    :param flux: fluxes for all emission lines to be modeled
    :param flux_err: flux uncertainties for all emission lines to be modeled
    :return: float, proportional to Bayesian posterior probability (prior*probability of observing data given model parameters)
    """
    start_time = time.time() #Want to time how long it takes for each model to run
    lp = lnprior(theta, theta_keys, prior_dict)
    if not np.isfinite(lp):
        return -np.inf
    else:
        posterior = lp + lnlike(theta, velocity, flux, flux_err, theta_keys, param_dict, inds)
        #print '%s seconds' % (time.time()-start_time)
        return posterior
#########################################################################################################################
#CONSTANTS
lines_to_fit_dict = {'[0,1]': [1460.17], '[1,7]': [1467.08]}

#READ IN DATA - ALL FILES WILL BE LOADED INTO DICTIONARIES
targ = "V4046Sgr"
#Disk Constants
targ_d = 72.3 #Gaia DR 2 (Bailer-Jones et al. 2018)
targ_AV = 0. #France et al. 2017 (HI corrected)
targ_M = 0.86+0.69 #France et al. 2017
targ_inclination = 33.

target_wavelengths = pd.read_pickle("../target_wavelengths_RULupi.pkl")
target_velocity = pd.read_pickle('../RULupi_H2velocity.pkl')
target_fluxes = pd.read_pickle('../RULupi_H2flux.pkl')
target_fluxerr = pd.read_pickle('../RULupi_H2fluxerr.pkl')

#Convert dictionaries to arrays
model_refwaves = np.array(target_wavelengths["labwave"]).astype(np.float64) #Wavelengths that models were created for

line_params = pd.read_csv("Hoadley2015_Table2_V4046Sgr.csv", sep=",") #Table with line properties, observed wavelengths/velocity shifts (measured with GUI)
ref_waves = np.array(line_params["Wavelength"]).astype(np.float64) #Rest wavelengths, from Hoadley 2015
obs_waves = np.array(line_params["Observed Wavelength"]).astype(np.float64) #Observed wavelengths
obs_vel = np.array(line_params["Velocity Shift"]).astype(np.float64) #Observed velocities
abs_vel = 150. #Velocity range of emission line to fit (+/-150 km/s)

#Need to shift lines according to measured velocity (values from GUI measurements) - this is done for all data, regardless of whether or not we'll actually fit the line
for wave_idx, wave in enumerate(model_refwaves):

    keri_loc = find_wave_loc(wave, ref_waves) #Find where difference between reference wavelengths in spreadsheet and model reference wavelengths is minimized (i.e. where the reference wavelength is equal to the model wavelength)
    v_shift = obs_vel[keri_loc] #Pull correct velocity shift from spreadsheet

    if np.isnan(v_shift) == False:
        key = str(round(wave, 2))[0:7]
        #Need to fit continuum/shift lines
        line_indices = np.argwhere(np.absolute(target_velocity[key]) <= abs_vel).flatten()
        linevelocity = np.array(target_velocity[key] - v_shift).astype(np.float64)

        cont_indices = np.argwhere(np.absolute(linevelocity) >= 150).flatten()

        flux_line = np.array(target_fluxes[key]).astype(np.float64)
        flux_err = np.array(target_fluxerr[key]).astype(np.float64)

        v_tofit = [linevelocity[idx] for idx in cont_indices]
        flux_tofit = [flux_line[idx] for idx in cont_indices]
        fluxerr_tofit = [flux_err[idx] for idx in cont_indices]

        p = np.polyfit(v_tofit, flux_tofit, 1) #1st order polynomial fit - this isn't great for some lines, but those lines are pretty noisy anyway

        continuum_fluxes = p[1] + (linevelocity*p[0])
        norm_fluxes = flux_line - continuum_fluxes

        #Uncomment this section to make sure the continuum-subtracted emission lines look good
        # plt.clf()
        # plt.close()
        # plt.title(key+" Continuum-Subtracted Flux")
        # plt.plot(linevelocity,
        #          norm_fluxes,
        #          color="k", lw=3.)
        # plt.plot(linevelocity,
        #          continuum_fluxes,
        #          color="r", lw=3.)
        # plt.show(block=True)

        target_fluxes[key] = norm_fluxes[line_indices]
        target_velocity[key] = linevelocity[line_indices]
        target_fluxerr[key] = flux_err[line_indices]

#Concatenate individual lines into single arrays for model fitting
prog_14_datavel, prog_14_dataflux, prog_14_dataerr = get_prog_data(lines_to_fit_dict['[1,7]'], model_refwaves, target_velocity, target_fluxes, target_fluxerr)

#Define set of model parameters for Python version of code
param_dict = {"Inclination": targ_inclination * math.pi / 180.,
              "flux_level": 0.,
              "lnf": 0.,
              "Molecule": "H2"
              }

#Read in reconstructed LyA profile
LyA_keys = ["LyA01flux", "LyA02flux", "LyA14flux", "LyA17flux"]
LyA_pumping = [1217.205, 1217.643, 1216.070, 1215.726]
LyA_file = "V4046Sgr_LyAprof_France2014.csv"
LyA_df = pd.read_csv(LyA_file)

for pumpwave_idx, pump_wave in enumerate(LyA_pumping):
    LyA_tokeep = np.array(LyA_df['ry_out'].loc[(LyA_df['lambda'] >= pump_wave - 0.01)
                                               & (LyA_df['lambda'] <= pump_wave + 0.01)])[0]
    param_dict[LyA_keys[pumpwave_idx]] = LyA_tokeep

prior_dict = {'z': [2., 7.],
              'gamma': [0., 1.99],
              'T': [500., 5000.],
              'q': [-2.5, 2.5],
              'rchar': [0., 50.],
              'MH2': [1.e-16, 1.e-1],
              'inclination': [0., 89.],
              'LyA01flux': [1.e-11, 3.e-11], #LyA bounds set for RY Lupi (see April 7-8 notes for the rest)
              'LyA02flux': [9.e-12, 2.e-11],
              'LyA14flux': [4.e-11, 1.e-10],
              'LyA17flux': [2.e-11, 7.e-11],
              'flux_level': [1.e-16, 1.e-9],
              'lnf': [-5., 5.]
              }

theta_keys = ['z', 'gamma', 'T', 'q', 'rchar', 'MH2']

grid_dict = grid_from_radius(targ_d, param_dict)
vobs = vrad_to_vobs(grid_dict["rgrid"], grid_dict["phigrid"])

vobs_flat = vobs.flatten() #same as np.ravel with order 'C'
bins, edges = np.histogram(vobs_flat, bins=599, range=(-300., 300.))
inds = np.digitize(vobs_flat, edges)
v_obs_bin = np.arange(-300., 300., 1.)

#Grid of parameters to search
z_vals = np.linspace(2,7,5)#[5., 6., 7.]
gamma_vals = np.linspace(0,2,5)#[1.75]
T_vals = [2200.]#[1500., 2200.]
q_vals = [0.5]#[-0.5, 0.5]
rchar_vals = [0.03]
MH2_vals = np.linspace(1e-5,1e-1,5)#[5.e-10]

#Loop through grid of parameters, save chi^2 values to file
for param_list in tqdm(list(product(z_vals, gamma_vals, T_vals, q_vals, rchar_vals, MH2_vals))):
    theta = list(param_list)

    log_like = lnprob(theta, prog_14_datavel, prog_14_dataflux, prog_14_dataerr, theta_keys, param_dict, prior_dict, inds)
    theta.insert(0, log_like)
    str_to_write = "/".join([str(val) for val in theta])
    outF = open("%s_chi2_H2.txt" % targ, "a")
    outF.write(str_to_write)
    outF.write("\n")
    outF.close()

#Read in file with chi^2 values saved, to pull out the best-fit model
data_file = "%s_chi2_H2.txt" % targ
data_table = pd.read_csv(data_file, sep="/", header=None)
data_table.columns = ["chi^2", "z", "gamma", "T", "q", "rchar", "MH2"]


sorted_data_table = data_table.assign(x=data_table['chi^2'].replace(-np.inf, np.nan)).sort_values("x", ascending=True, na_position='last').drop('x',axis=1)
index_sorted = sorted_data_table.index
best_fit_idx = index_sorted[0]
print (sorted_data_table.loc[best_fit_idx])

theta = [sorted_data_table.loc[best_fit_idx]["z"], sorted_data_table.loc[best_fit_idx]["gamma"], sorted_data_table.loc[best_fit_idx]["T"], sorted_data_table.loc[best_fit_idx]["q"], sorted_data_table.loc[best_fit_idx]["rchar"], sorted_data_table.loc[best_fit_idx]["MH2"]]
full_interp_model = get_H2model(theta, theta_keys, param_dict, inds)

plt.clf()
plt.close()
plt.title("V4046 Sgr: Best-Fit H2 Model Emission Lines")
plt.plot(np.arange(len(prog_14_dataflux)),
         prog_14_dataflux,
         color="k", lw=2.)
plt.plot(np.arange(len(prog_14_dataflux)),
         full_interp_model,
         color="teal", lw=1.)
plt.gca().minorticks_on()
plt.show(block=True)
