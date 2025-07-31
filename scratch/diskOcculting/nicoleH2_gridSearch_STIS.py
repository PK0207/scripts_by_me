# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 18:21:41 2025

@author: ecruzaguirre
"""

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from itertools import product
from scipy.interpolate import griddata
#Import custom modules
from ccm_unred import ccm_unred
from OOP_HoadleyH2model import H2Model
from linefittingfunctions import stis_lsf_cruzaguirre

def loadData(dataFile):
    raw=np.load(dataFile,allow_pickle=True)
    dat=raw['data'].item() #returns dictionary
    lab,obs,off,prg,lya,vel,flx,err=[],[],[],[],[],[],[],[]
    ltf={} #lines to fit, progressions as keys, with lab wavelengths
    for key in dat.keys():
        one=dat[key] #load the data list for each emission line
        lab.append(one[0])
        obs.append(one[1])
        off.append(one[2])
        prg.append(one[3])
        lya.append(one[4])
        vel.append(one[5])
        flx.append(one[6])
        err.append(one[7])
        if one[3] not in ltf.keys():
            ltf[one[3]]=[one[0]] #create new key
        else:
            ltf[one[3]].append(one[0]) #add to existing key
    return lab,obs,off,prg,lya,vel,flx,err,ltf

def grid_from_radius(targ_d, param_dict):
    """Builds grids to describe disk based on stellar distance and disk inclination
    """
    rinout = 10.0**(0.02 * np.arange(0, 201, 1) - 1.4)
    r_grid = rinout / targ_d
    phi_grid = np.arange(0.0, 360.0, 2.0) * (np.pi / 180.0)

    xgrid = np.zeros((len(r_grid), len(phi_grid)))
    ygrid = np.zeros((len(r_grid), len(phi_grid)))

    tan_phi = np.tan(phi_grid)

    for r_idx, r in enumerate(r_grid):
        for az_idx, az in enumerate(tan_phi):
            if az >= 90.0 * (np.pi / 180.0) and az <= 270.0 * (np.pi / 180.0):
                x_dummy = -1.0 * np.sqrt(r ** 2.0 / (1.0 + az))
            else:
                x_dummy = np.sqrt(r ** 2.0 / (1.0 + az ** 2.0))
            xgrid[r_idx, az_idx] = x_dummy / np.cos(param_dict["Inclination"])

            y_dummy = np.sqrt(r ** 2.0 / (1.0 + (1.0 / az ** 2.0)))
            if phi_grid[az_idx] >= np.pi:
                ygrid[r_idx, az_idx] = -1.0 * y_dummy
            else:
                ygrid[r_idx, az_idx] = y_dummy

    r_H2 = np.zeros((len(r_grid), len(phi_grid)))
    phi_H2 = np.zeros((len(r_grid), len(phi_grid)))
    for r_idx, r in enumerate(r_grid):
        for az_idx, az in enumerate(tan_phi):
            r_dummy = np.sqrt(xgrid[r_idx, az_idx] ** 2.0 + ygrid[r_idx, az_idx] ** 2.0) * targ_d
            r_H2[r_idx, az_idx] = r_dummy

            phi_dummy = np.arctan(ygrid[r_idx, az_idx] / xgrid[r_idx, az_idx])
            if ygrid[r_idx, az_idx] <= 0.0:
                if xgrid[r_idx, az_idx] <= 0.0:
                    phi_H2[r_idx, az_idx] = phi_dummy + np.pi
                else:
                    phi_H2[r_idx, az_idx] = phi_dummy + 2.0*np.pi
            else:
                if xgrid[r_idx, az_idx] <= 0.0:
                    phi_H2[r_idx, az_idx] = phi_dummy + np.pi
                else:
                    phi_H2[r_idx, az_idx] = phi_dummy
    return {"xgrid": xgrid,
            "ygrid": ygrid,
            "rgrid": (r_grid * targ_d) - 0.03,
            "phigrid": phi_grid * (180.0 / np.pi),
            "rH2": r_H2,
            "phiH2": phi_H2 * (180.0 / np.pi)}

def vrad_to_vobs(r_grid, phi_grid):
    v_k = np.sqrt(G * targ_M * MSUNtoG / (r_grid * AUtoCM)) * CMtoKM
    vobs = [v_k * np.sin(param_dict["Inclination"]) * np.sin(phi * (np.pi / 180.0)) for phi in phi_grid]
    vobs = np.array(vobs).astype(np.float64)
    vobs = np.swapaxes(vobs, 0, 1)
    return vobs

def collapsed_line_profile(lineprofs, param_dict, inds):

    total_line_profile = np.full((len(v_obs_bin), len(lineprofs[:, 0, 0])),1.0e-75,dtype=np.float64) #initialize array with dummy value, shape (RV bin, pump_prog)
    for prog_idx in range(0, np.shape(lineprofs)[0]):
        flat_fluxes = lineprofs[prog_idx, :, :].flatten()
        test_df = pd.DataFrame({"Velocities": inds,"Fluxes": flat_fluxes})
        g = test_df.groupby(["Velocities"])
        profile = np.array(g.sum()).astype(np.float64).flatten()
        pads = (4*abs_vel - len(profile))

        if pads % 2 == 0:
            pad_left = pads / 2
            pad_right = pads / 2
        else:
            pad_left = int(np.floor(pads / 2))
            pad_right = int(np.ceil(pads / 2))

        final_profile = np.pad(profile, (pad_left, pad_right), 'edge')

        total_line_profile[:, prog_idx] = final_profile #save the profile
        total_line_profile[2*abs_vel, prog_idx] = total_line_profile[2*abs_vel+1, prog_idx] #!!! is this to make the profile symmetric?

    return total_line_profile

def vradWave(rad,cen):
    c=2.998e5 #speed of light in km/s
    w=cen*(1+rad/c)
    return w

def prepLSF(orig_lsf_wave,orig_lsf,data_wave):
    '''
    prepare the lower resolution STIS LSF for use with the higher resolution model 
    '''
    data_wave_spacing = data_wave[1]-data_wave[0]
    data_wave_length = len(data_wave)
    lsf_lam_min = np.round(np.min(orig_lsf_wave)/data_wave_spacing) * data_wave_spacing
    lsf_lam_onesided = np.arange(lsf_lam_min,0,data_wave_spacing)  ### Make sure it's even and doesn't include zero
    if len(lsf_lam_onesided) % 2 != 0:
        lsf_lam_onesided = lsf_lam_onesided[1::] # get rid of the first element of the array
    
    lsf_lam_flipped = lsf_lam_onesided[::-1]
    lsf_lam_pluszero=np.append(lsf_lam_onesided,np.array([0]))
    lsf_lam=np.append(lsf_lam_pluszero,lsf_lam_flipped) # should be odd
    
    lsf_interp = np.interp(lsf_lam,orig_lsf_wave,orig_lsf/np.sum(orig_lsf))
    lsf_interp_norm = lsf_interp/np.sum(lsf_interp) 
    
    if data_wave_length < len(lsf_interp_norm): #if this is true, returned convolution has the length of the lsf, not the data, which is bad
        lsf_interp_norm = np.delete(lsf_interp_norm,np.where(lsf_interp_norm == 0)) #get rid of all zeros on sides
        lsf_interp_norm = np.insert(lsf_interp_norm,0,0) #add one 0.0 on each side
        lsf_interp_norm = np.append(lsf_interp_norm,0) 
        #this may not always be enough, but these cuts do not affect the convolution
    return lsf_interp_norm

def get_H2model(theta, theta_keys, param_dict, inds, plot=False, retT=False):
    """Concatenate model emission lines into single array for fitting to data
    :param params: List of model parameters (floats)
    :return: full_interp_model (array), array of model emission lines
    """
    
    for key_idx, key in enumerate(theta_keys):
        param_dict[key] = theta[key_idx] #add the current set of parameters to the param dictionary

    #pump the H2 at four LyA points, four progressions which share the same upper state
    target = H2Model(targ_M, targ_AV, targ_d, param_dict) #modeled disk
    num_den = target.H2_number_density() #number density (pump_prog, z, 201)
    lineprof = target.total_intensity() #shape (pump_prog, r, phi)
    rvals = target.grid_dict["rgrid"] #radial disk array
    Tprof = target.radial_T(rvals) #temperature profile

    if plot:
        fig,ax=plt.subplots()
        plt.xscale("log")
        plt.title("Radial Temperature Distribution: q = "+str(param_dict["q"]))
        plt.plot(rvals,Tprof,color="k",lw=2.0)
        plt.xlim(0.01, 15.0)
        plt.gca().minorticks_on()
        plt.show()
        
        fig,ax=plt.subplots()
        plt.xscale("log")
        plt.title("Number Density Distribution: z = 0")
        for p in range(len(num_den)):
            plt.plot(rvals,num_den[p,0,:],lw=2.0,label=LyA_prog[p])
        plt.xlim(0.01, 15.0)
        plt.legend()
        plt.gca().minorticks_on()
        plt.show()

    intensity_binned = collapsed_line_profile(lineprof, param_dict, inds) #intensity as a function of disk position, shape (RV bin, pump_prog)
    intensity = np.zeros((len(v_obs_bin), len(target.lambda_props["Wavelength"]))) #velocity distribution for each progression (RV bin, lambda_props lines)

    for wave_idx, wave in enumerate(target.lambda_props["Wavelength"]): 
        #grab the flux vs radial velocity data from the disk for each emission line in the H2 model
        prog_idx = np.argmin(np.absolute(target.lambda_props["Jl"][wave_idx] - target.J_initial)) #select which of the four progressions this line is part of
        intens_unred = intensity_binned[:, prog_idx] * target.lambda_props["Bul"][wave_idx] #get the intensity contribution of the total intensity from just this line
        intens_red = ccm_unred(np.zeros(len(intens_unred)) + wave, intens_unred, EBV) #redden the intensity profile
        intens_red[np.argwhere(np.isnan(intens_red) == True).flatten()] = 1.0e-30 #remove any nans

        #Smooth/convolve with LSF
        lsfx, lsfy = stis_lsf_cruzaguirre(1500,stisGrat,stisSlit,True) #get wave vs norm LSF
        lsfy_prep = prepLSF(lsfx,lsfy,vradWave(v_obs_bin,wave)) #prepare the LSF for the model wavegrid
        if len(lsfy_prep)>len(v_obs_bin):
            #if the LSF is too large, the resulting profile will not be the correct length
            diffLen=len(lsfy_prep)-len(v_obs_bin) #find the difference in length
            lsfy_trnc=lsfy_prep[int(np.ceil(diffLen/2.0)):-int(np.ceil(diffLen/2.0))] #truncate, at the cost of flux loss
            intensity[:, wave_idx] = np.convolve(intens_red, lsfy_trnc, mode="same")
        else:
            intensity[:, wave_idx] = np.convolve(intens_red, lsfy_prep, mode="same")
        
    full_interp_model = []
    for idx,wave in enumerate(model_refwaves):
        try: #of the emission lines in the model, grab the ones that we have data for
            mod_idx=np.where(wave==target.lambda_props["Wavelength"])[0][0] #grab the modeled emission line
            smoothflux = intensity[:, mod_idx]
            interp_smoothflux = griddata(v_obs_bin, smoothflux, target_velocityL[idx], method = "nearest")
            full_interp_model.append(interp_smoothflux)
        except(IndexError):
            #if the line is not in the model's list of emission lines, return zeros
            full_interp_model.append(np.zeros(len(target_velocityL[idx])))
    
    if retT:
        return np.concatenate(full_interp_model),target
    else:            
        return np.concatenate(full_interp_model) #return the emission lines in the same format as the imput emission lines

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

    lp = lnprior(theta, theta_keys, prior_dict) #check that the parameters are within the model bounds
    if not np.isfinite(lp):
        return -np.inf #if out of bounds, don't run the model
    else: #otherwise, run the model
        #Log-likelihood function for H2 models (assuming Gaussian distribution, for lack of a better guess)
        full_interp_model = get_H2model(theta, theta_keys, param_dict, inds)
        loglike = np.sum(np.square((flux - full_interp_model) / flux_err)) / (len(theta) - 1)
        posterior = lp + loglike
        return posterior
        



#Global constants
G = 6.67259e-8 #Gravitational constant [cm^3 g^-1 s^-2]
c = 3.0e10 #Speed of light [cm/s]
h = 6.626e-27 #Planck constant [erg s]
kB = 1.38e-16 #Boltzmann constant [erg/K]
AUtoCM = 1.496e13 #Convert [AU -> cm]
CMtoKM = 1.0e-5 #Convert [cm -> km]
PCtoAU = 206265.0 #Convert [parsecs -> astronomical units]
MSUNtoG = 1.99e33 #Convert [M_solar -> grams]
ANGtoCM = 1.0e-8 #Convert [Angstroms -> cm]
Rv = 3.1 #Extinction dependent on line-of-sight, 3.1 is the standard value for diffuse ISM targets

#error silencing
np.seterr(divide='ignore') #divide by 0 silencing
np.seterr(invalid='ignore') #correct_factor silencing

#Disk Constants
targ = "PDS 70"
targ_d = 112 
targ_AV = 0.05 
targ_M = 0.76 
targ_inclination = 52.0
EBV = -1.0 * targ_AV / Rv #reddening/attenuation

#read in the disk data
preppedData='PDS70_STIS_G140L_MAST_Data.npz'
ref_waves,obs_waves,obs_vel,progTest,lyawTest,target_velocity,target_fluxes,target_fluxerr,lines_to_fit_dict=loadData(preppedData)
model_refwaves=ref_waves #check if deepcopy is needed

#observation details
stisGrat='G140L'
stisSlit='52x0.2' #from the x1d file


#Need to shift lines according to measured velocity, done for all data regardless of whether or not we'll actually fit the line
abs_vel=500 #cutoff between line and continuum
target_velocityL,target_fluxesL,target_fluxerrL=[],[],[] #prepared line data
for v_idx,v_shift in enumerate(obs_vel):
    if np.isnan(v_shift) == False:
        #Shift the lines
        linevelocity = target_velocity[v_idx] - v_shift #shifted line velocities

        #isolate the emission line
        line_indices = np.argwhere(np.absolute(linevelocity) <= abs_vel).flatten() #indices covering the emisison line
        cont_indices = np.argwhere(np.absolute(linevelocity) >= abs_vel).flatten() #indices in the b/g continuum

        #grab the flux data
        flux_line = target_fluxes[v_idx]
        flux_err = target_fluxerr[v_idx]

        #fit the continuum
        v_tofit = [linevelocity[idx] for idx in cont_indices]
        flux_tofit = [flux_line[idx] for idx in cont_indices]
        p = np.polyfit(v_tofit, flux_tofit, 1) #1st order polynomial fit
        continuum_fluxes = p[1] + (linevelocity*p[0])
        # fig,ax=plt.subplots()
        # plt.title(f"{ref_waves[v_idx]} Continuum Fitting")
        # plt.plot(v_tofit,flux_tofit,color="k",lw=3)
        # plt.plot(linevelocity,continuum_fluxes,color='r',lw=3)
        # plt.show()
        norm_fluxes = flux_line - continuum_fluxes #subtract the background

        # #Uncomment this section to make sure the continuum-subtracted emission lines look good
        # fig,ax=plt.subplots()
        # plt.title(f"{ref_waves[v_idx]} Continuum-Subtracted Flux")
        # plt.plot(linevelocity,flux_line,color="k",lw=3)
        # plt.plot(linevelocity,norm_fluxes,color="g",lw=3)
        # plt.plot(linevelocity,continuum_fluxes,color="r",lw=3)
        # plt.show()

        target_velocityL.append(linevelocity[line_indices])
        target_fluxesL.append(norm_fluxes[line_indices])
        target_fluxerrL.append(flux_err[line_indices]) #only save the line data
        
#Concatenate individual lines into single arrays for model fitting
all_datavel=np.concatenate(target_velocityL)
all_dataflux=np.concatenate(target_fluxesL)
all_dataerr=np.concatenate(target_fluxerrL)

#Define set of model parameters for Python version of code
param_dict = {"Inclination": targ_inclination * (np.pi / 180.0),
              "flux_level": 0.,
              "lnf": 0.,
              "Molecule": "H2"}

#Read in reconstructed LyA profile
LyA_keys = ["LyA01flux", "LyA02flux", "LyA14flux", "LyA17flux"] #keys corresponding to pumping wavelengths for each progression
LyA_prog=['[0,1]','[0,2]','[1,4]','[1,7]']
LyA_pumping = [1217.205, 1217.643, 1216.070, 1215.726]
LyA_file = "V4046Sgr_LyAprof_France2014.csv" 
LyA_df = pd.read_csv(LyA_file) #profile containing the intrinsic emission and outflow absorbed emission

for pumpwave_idx, pump_wave in enumerate(LyA_pumping): #pull the corresponding flux from the model line profile (as seen from Earth)
    LyA_tokeep = np.array(LyA_df['ry_out'].loc[(LyA_df['lambda'] >= pump_wave - 0.01)
                                               & (LyA_df['lambda'] <= pump_wave + 0.01)])[0]
    param_dict[LyA_keys[pumpwave_idx]] = LyA_tokeep
    
#define which parameters will be varied within the model
#theta_keys = ['z', 'gamma', 'T', 'q', 'rchar', 'MH2']
theta_keys = ['z', 'gamma', 'T', 'q', 'rchar', 'MH2', 'LyA14flux', 'LyA17flux']
    
#Set bounds for each model parameter
# prior_dict = {'z': [2.0, 7.0],
#               'gamma': [0.0, 1.99],
#               'T': [500., 5000.],
#               'q': [-2.5, 2.5],
#               'rchar': [0.0, 50.0],
#               'MH2': [1.0e-16, 1.0e-1],
#               'inclination': [0.0, 89.0],
#               'LyA01flux': [1.0e-11, 3.0e-11], #LyA bounds set for RY Lupi
#               'LyA02flux': [9.0e-12, 2.0e-11],
#               'LyA14flux': [4.0e-11, 1.0e-10],
#               'LyA17flux': [2.0e-11, 7.0e-11],
#               'flux_level': [1.0e-16, 1.0e-9],
#               'lnf': [-5.0, 5.0]}
prior_dict = {'z': [2.0, 7.0],
              'gamma': [0.0, 1.99],
              'T': [500., 5000.],
              'q': [-2.5, 2.5],
              'rchar': [0.0, 50.0],
              'MH2': [1.0e-16, 1.0e-1],
              'inclination': [0.0, 89.0],
              'LyA01flux': [1.0e-11, 3.0e-11], #LyA bounds set for RY Lupi
              'LyA02flux': [9.0e-12, 2.0e-11],
              'LyA14flux': [1.0e-12, 1.0e-10],
              'LyA17flux': [1.0e-12, 1.0e-10],
              'flux_level': [1.0e-16, 1.0e-9],
              'lnf': [-5.0, 5.0]}

#initialize the PPD grid and associated properties
grid_dict = grid_from_radius(targ_d, param_dict) #contains the cartesian (x,y) and cylindrical (r,phi) coordinates for each disk grid point (r: 201, phi:180)
vobs = vrad_to_vobs(grid_dict["rgrid"], grid_dict["phigrid"]) #the radial velocity at each grid point
vobs_flat = vobs.flatten() #flatten disk RV into a 1D array of len 201*180 
bins, edges = np.histogram(vobs_flat, bins=599, range=(-2*abs_vel, 2*abs_vel)) #create a histogram of disk RV values
inds = np.digitize(vobs_flat, edges) #np.digitize creates an array of which bin each value falls within
v_obs_bin = np.arange(-2*abs_vel, 2*abs_vel, 1.0) #bin values (cleaner than edges --> whole numbers) 

#Grid of parameters to search
#on average on my computer, ~3.5it/s or ~0.286s/it
#2**6 --> 20 sec
#3**6 --> 4 min
#4**6 --> 20 min
#5**6 --> 75 min
#10**6 --> 80 hrs
res=2 #resolution for all parameters
z_vals = np.linspace(2,7,res)
gamma_vals = np.linspace(0,1.99,res)
T_vals = np.linspace(500,5000,res)
q_vals = np.linspace(-2.5,2.5,res)
rchar_vals = np.linspace(0.0,50.0,res)
MH2_vals = np.logspace(-16,-1,res) #logspace will space things ount nicely within the range of input powers of 10
LyA14flux = np.logspace(-12,-10,res)
LyA17flux = np.logspace(-12,-10,res)
# z_vals = [6]
# gamma_vals = [1.75]
# T_vals = [2500.0]
# q_vals = [0.5]
# rchar_vals = np.linspace(0.0,50.0,res)
# MH2_vals = np.logspace(-16,-1,res) #logspace will space things ount nicely within the range of input powers of 10

#Loop through grid of parameters, save chi^2 values to file
timeRaw=time.localtime()
timeFmt=time.strftime('%b-%d-%Y_%I-%M_%p',timeRaw) #generate one output file per model run
#test set
z_vals = [6.0]
gamma_vals = [1.75]
T_vals = [2200]
q_vals = [0.5]
rchar_vals = [0.03]
MH2_vals = [10**-9.82]
LyA14flux = [10**-11]
LyA17flux = [10**-11]
timeFmt='TEST'
for param_list in tqdm(list(product(z_vals, gamma_vals, T_vals, q_vals, rchar_vals, MH2_vals, LyA14flux, LyA17flux))):
    #create a list of each unique combination of parameters
    theta = list(param_list)
    
    #run the model and save the data
    log_like = lnprob(theta, all_datavel, all_dataflux, all_dataerr, theta_keys, param_dict, prior_dict, inds)
    theta.insert(0, log_like)
    str_to_write = "/".join([str(val) for val in theta])
    outF = open(f"PDS70_modeling_chi2_H2_{timeFmt}.txt",'a') #save model results to output file
    outF.write(str_to_write+'\n')
    outF.close()
    
#Read in file with chi^2 values saved, to pull out the best-fit model
data_file = f"PDS70_modeling_chi2_H2_{timeFmt}.txt"
#data_cols = ["chi^2", "z", "gamma", "T", "q", "rchar", "MH2"]
data_cols = ["chi^2", "z", "gamma", "T", "q", "rchar", "MH2", "LyA14flux", "LyA17flux"]
data_table = pd.read_csv(data_file, sep="/", header=None, names=data_cols)
closeChi=1.0 #pick how close to a specific chi^2 value you want to be
closeChi=200
best_fit_idx=np.argmin(np.abs(data_table['chi^2'].to_numpy() - closeChi)) #index of the best fit (reduced chi^2 closest to 1) [don't need to worry about -inf]
for i in range(len(data_table.loc[best_fit_idx])):
    print(f'{data_cols[i]: <5} : {data_table.loc[best_fit_idx].iloc[i]}') #pretty print
#print(data_table.loc[best_fit_idx])
# thetaB = [data_table.loc[best_fit_idx]["z"],data_table.loc[best_fit_idx]["gamma"],
#           data_table.loc[best_fit_idx]["T"],data_table.loc[best_fit_idx]["q"],
#           data_table.loc[best_fit_idx]["rchar"],data_table.loc[best_fit_idx]["MH2"]]
thetaB = [data_table.loc[best_fit_idx]["z"],data_table.loc[best_fit_idx]["gamma"],
          data_table.loc[best_fit_idx]["T"],data_table.loc[best_fit_idx]["q"],
          data_table.loc[best_fit_idx]["rchar"],data_table.loc[best_fit_idx]["MH2"],
          data_table.loc[best_fit_idx]["LyA14flux"],data_table.loc[best_fit_idx]["LyA17flux"]]
full_interp_modelB = get_H2model(thetaB, theta_keys, param_dict, inds)
loglikeB = np.sum(np.square((all_dataflux - full_interp_modelB) / all_dataerr)) / (len(thetaB) - 1)

fig,ax=plt.subplots()
plt.title("PDS 70: Best-Fit H2 Model Emission Lines")
plt.plot(np.arange(len(all_dataflux)),all_dataflux,color="k",lw=2.0,label='Data')
plt.fill_between(np.arange(len(all_dataflux)),all_dataflux+all_dataerr,all_dataflux-all_dataerr,alpha=0.35,color='k')
plt.plot(np.arange(len(all_dataflux)),full_interp_modelB,color="teal", lw=1.0,label='Model')
plt.legend()
plt.gca().minorticks_on()
plt.show()


'''
To Do
-----
* Introduce lmfit to find a best fit
* Use the MCMC capabilities of lmfit to run an mcmc about the best fit 
* Get Steps 3-7 running and then backtrack to fix things
'''