# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 14:19:59 2025

@author: ecruzaguirre
"""

import time
import corner #make fancy parameter plots
import pickle #to save the mcmc results so I don't have to do it over again
import arviz as az #convert pandas dataframe to a data type compatible with corner (corner+pandas is now deprecated)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from lmfit import Model
from scipy.interpolate import griddata
#Import custom modules
from ccm_unred import ccm_unred
from OOP_HoadleyH2model import H2Model
from linefittingfunctions import cos_lsf_arulanantham_ullyses

def loadData(dataFile):
    raw=np.load(dataFile,allow_pickle=True)
    dat=raw['data'].item() #returns dictionary
    lab,off,prg,lya,vel,flx,err=[],[],[],[],[],[],[]
    ltf={} #lines to fit, progressions as keys, with lab wavelengths
    for key in dat.keys():
        one=dat[key] #load the data list for each emission line
        lab.append(one[0])
        off.append(one[1])
        prg.append(one[2])
        lya.append(one[3])
        vel.append(one[4])
        flx.append(one[5])
        err.append(one[6])
        if one[2] not in ltf.keys():
            ltf[one[2]]=[one[0]] #create new key
        else:
            ltf[one[2]].append(one[0]) #add to existing key
    return lab,off,prg,lya,vel,flx,err,ltf

def grid_from_radius(targ_d, inc):
    """Builds grids to describe disk based on stellar distance and disk inclination
    """
    rinout = 10.0**(0.02 * np.arange(0, 201, 1) - 1.4)
    r_grid = rinout / targ_d
    phi_grid = np.arange(0.0, 360.0, 2.0) * (np.pi / 180.0)

    xgrid = np.zeros((len(r_grid), len(phi_grid)))
    ygrid = np.zeros((len(r_grid), len(phi_grid)))

    tan_phi = np.tan(phi_grid)

    for r_idx, r in enumerate(r_grid):
        for azi_idx, azi in enumerate(tan_phi):
            if azi >= 90.0 * (np.pi / 180.0) and azi <= 270.0 * (np.pi / 180.0):
                x_dummy = -1.0 * np.sqrt(r ** 2.0 / (1.0 + azi))
            else:
                x_dummy = np.sqrt(r ** 2.0 / (1.0 + azi ** 2.0))
            xgrid[r_idx, azi_idx] = x_dummy / np.cos(inc * (np.pi / 180))

            y_dummy = np.sqrt(r ** 2.0 / (1.0 + (1.0 / azi ** 2.0)))
            if phi_grid[azi_idx] >= np.pi:
                ygrid[r_idx, azi_idx] = -1.0 * y_dummy
            else:
                ygrid[r_idx, azi_idx] = y_dummy

    r_H2 = np.zeros((len(r_grid), len(phi_grid)))
    phi_H2 = np.zeros((len(r_grid), len(phi_grid)))
    for r_idx, r in enumerate(r_grid):
        for azi_idx, azi in enumerate(tan_phi):
            r_dummy = np.sqrt(xgrid[r_idx, azi_idx] ** 2.0 + ygrid[r_idx, azi_idx] ** 2.0) * targ_d
            r_H2[r_idx, azi_idx] = r_dummy

            phi_dummy = np.arctan(ygrid[r_idx, azi_idx] / xgrid[r_idx, azi_idx])
            if ygrid[r_idx, azi_idx] <= 0.0:
                if xgrid[r_idx, azi_idx] <= 0.0:
                    phi_H2[r_idx, azi_idx] = phi_dummy + np.pi
                else:
                    phi_H2[r_idx, azi_idx] = phi_dummy + 2.0*np.pi
            else:
                if xgrid[r_idx, azi_idx] <= 0.0:
                    phi_H2[r_idx, azi_idx] = phi_dummy + np.pi
                else:
                    phi_H2[r_idx, azi_idx] = phi_dummy
    return {"xgrid": xgrid,
            "ygrid": ygrid,
            "rgrid": (r_grid * targ_d) - 0.03,
            "phigrid": phi_grid * (180.0 / np.pi),
            "rH2": r_H2,
            "phiH2": phi_H2 * (180.0 / np.pi)}

def vrad_to_vobs(r_grid, phi_grid, inc):
    v_k = np.sqrt(G * targ_M * MSUNtoG / (r_grid * AUtoCM)) * CMtoKM
    vobs = [v_k * np.sin(inc * (np.pi / 180)) * np.sin(phi * (np.pi / 180.0)) for phi in phi_grid]
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
        pads = (2*abs_vel_mlt*abs_vel - len(profile))

        if pads > 0:
            if pads % 2 == 0:
                pad_left = pads / 2
                pad_right = pads / 2
            else:
                pad_left = int(np.floor(pads / 2))
                pad_right = int(np.ceil(pads / 2))
        else:
            #can't pad if it's negative
            pad_left=0
            pad_right=0
            

        final_profile = np.pad(profile, (pad_left, pad_right), 'edge')

        total_line_profile[:, prog_idx] = final_profile #save the profile
        total_line_profile[abs_vel_mlt*abs_vel, prog_idx] = total_line_profile[abs_vel_mlt*abs_vel+1, prog_idx] #!!! is this to make the profile symmetric?

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

def get_H2model_all(vradData,LyA01,LyA02,LyA14,LyA17,z,gamma,T,q,rchar,MH2):
    """Concatenate model emission lines into single array for fitting to data
    :param params: List of model parameters (floats)
    :return: full_interp_model (array), array of model emission lines
    """
    
    #setup a param_dict, to be compatible with Nicole's model framework
    param_dict = {"Inclination": targ_inclination * (np.pi / 180.0),
                  "flux_level": 0.,
                  "lnf": 0.,
                  "Molecule": "H2",
                  "LyA01flux": 10.0**LyA01,
                  "LyA02flux": 10.0**LyA02,
                  "LyA14flux": 10.0**LyA14,
                  "LyA17flux": 10.0**LyA17,
                  "z": z,
                  "gamma": gamma,
                  "T": T,
                  "q": q,
                  "rchar": rchar,
                  "MH2": 10.0**MH2}

    #pump the H2 at four LyA points, four progressions which share the same upper state
    target = H2Model(targ_M, targ_AV, targ_d, param_dict) #modeled disk
    lineprof = target.total_intensity() #shape (pump_prog, r, phi)

    intensity_binned = collapsed_line_profile(lineprof, param_dict, inds) #intensity as a function of disk position, shape (RV bin, pump_prog)
    intensity = np.zeros((len(v_obs_bin), len(target.lambda_props["Wavelength"]))) #velocity distribution for each progression (RV bin, lambda_props lines)

    lsfXall,lsfYall=[],[]
    for wave_idx, wave in enumerate(target.lambda_props["Wavelength"]): 
        #grab the flux vs radial velocity data from the disk for each emission line in the H2 model
        prog_idx = np.argmin(np.absolute(target.lambda_props["Jl"][wave_idx] - target.J_initial)) #select which of the four progressions this line is part of
        intens_unred = intensity_binned[:, prog_idx] * target.lambda_props["Bul"][wave_idx] #get the intensity contribution of the total intensity from just this line
        intens_red = ccm_unred(np.zeros(len(intens_unred)) + wave, intens_unred, EBV) #redden the intensity profile
        intens_red[np.argwhere(np.isnan(intens_red) == True).flatten()] = 1.0e-30 #remove any nans

        #Smooth/convolve with LSF
        # lsfx, lsfy = stis_lsf_cruzaguirre(1500,stisGrat,stisSlit,True) #get wave vs norm LSF
        # lsfy_prep = prepLSF(lsfx,lsfy,vradWave(v_obs_bin,wave)) #prepare the LSF for the model wavegrid
        # if len(lsfy_prep)>len(v_obs_bin):
        #     #if the LSF is too large, the resulting profile will not be the correct length
        #     diffLen=len(lsfy_prep)-len(v_obs_bin) #find the difference in length
        #     lsfy_trnc=lsfy_prep[int(np.ceil(diffLen/2.0)):-int(np.ceil(diffLen/2.0))] #truncate, at the cost of flux loss
        #     intensity[:, wave_idx] = np.convolve(intens_red, lsfy_trnc, mode="same")
        # else:
        #     intensity[:, wave_idx] = np.convolve(intens_red, lsfy_prep, mode="same")
        lsfx, lsfy = cos_lsf_arulanantham_ullyses(wave, "LTP1", False) #!!!replace with my version of the LSF prep, need LP1 for G130M and G160M 
        lsfy_norm = lsfy.flatten() / np.sum(lsfy.flatten()) #Normalize LSF
        intensity[:, wave_idx] = np.convolve(intens_red, lsfy_norm, mode="same")
        lsfXall.append(lsfx)
        lsfYall.append(lsfy_norm)
        
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
    
    return np.concatenate(full_interp_model),full_interp_model,lineprof,lsfXall,lsfYall,param_dict,target

def get_H2model(vradData,LyA01,LyA02,LyA14,LyA17,z,gamma,T,q,rchar,MH2):
    """Concatenate model emission lines into single array for fitting to data
    :param params: List of model parameters (floats)
    :return: full_interp_model (array), array of model emission lines
    """
    
    #setup a param_dict, to be compatible with Nicole's model framework
    param_dict = {"Inclination": targ_inclination * (np.pi / 180.0),
                  "flux_level": 0.,
                  "lnf": 0.,
                  "Molecule": "H2",
                  "LyA01flux": 10.0**LyA01,
                  "LyA02flux": 10.0**LyA02,
                  "LyA14flux": 10.0**LyA14,
                  "LyA17flux": 10.0**LyA17,
                  "z": z,
                  "gamma": gamma,
                  "T": T,
                  "q": q,
                  "rchar": rchar,
                  "MH2": 10.0**MH2}

    #pump the H2 at four LyA points, four progressions which share the same upper state
    target = H2Model(targ_M, targ_AV, targ_d, param_dict) #modeled disk
    lineprof = target.total_intensity() #shape (pump_prog, r, phi)

    intensity_binned = collapsed_line_profile(lineprof, param_dict, inds) #intensity as a function of disk position, shape (RV bin, pump_prog)
    intensity = np.zeros((len(v_obs_bin), len(target.lambda_props["Wavelength"]))) #velocity distribution for each progression (RV bin, lambda_props lines)

    for wave_idx, wave in enumerate(target.lambda_props["Wavelength"]): 
        #grab the flux vs radial velocity data from the disk for each emission line in the H2 model
        prog_idx = np.argmin(np.absolute(target.lambda_props["Jl"][wave_idx] - target.J_initial)) #select which of the four progressions this line is part of
        intens_unred = intensity_binned[:, prog_idx] * target.lambda_props["Bul"][wave_idx] #get the intensity contribution of the total intensity from just this line
        intens_red = ccm_unred(np.zeros(len(intens_unred)) + wave, intens_unred, EBV) #redden the intensity profile
        intens_red[np.argwhere(np.isnan(intens_red) == True).flatten()] = 1.0e-30 #remove any nans

        #Smooth/convolve with LSF
        # lsfx, lsfy = stis_lsf_cruzaguirre(1500,stisGrat,stisSlit,True) #get wave vs norm LSF
        # lsfy_prep = prepLSF(lsfx,lsfy,vradWave(v_obs_bin,wave)) #prepare the LSF for the model wavegrid
        # if len(lsfy_prep)>len(v_obs_bin):
        #     #if the LSF is too large, the resulting profile will not be the correct length
        #     diffLen=len(lsfy_prep)-len(v_obs_bin) #find the difference in length
        #     lsfy_trnc=lsfy_prep[int(np.ceil(diffLen/2.0)):-int(np.ceil(diffLen/2.0))] #truncate, at the cost of flux loss
        #     intensity[:, wave_idx] = np.convolve(intens_red, lsfy_trnc, mode="same")
        # else:
        #     intensity[:, wave_idx] = np.convolve(intens_red, lsfy_prep, mode="same")
        lsfx, lsfy = cos_lsf_arulanantham_ullyses(wave, "LTP1", False) #!!!replace with my version of the LSF prep, need LP1 for G130M and G160M 
        lsfy_norm = lsfy.flatten() / np.sum(lsfy.flatten()) #Normalize LSF
        intensity[:, wave_idx] = np.convolve(intens_red, lsfy_norm, mode="same")
        
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
    
    return np.concatenate(full_interp_model) #return the emission lines in the same format as the input emission lines

def fit_H2Model(vradD,fluxD,errrD,ig,pb):
    '''
    use lmfit to find a best fit model to the data
    '''
    
    #define the model
    H2M=Model(get_H2model) 
    
    H2P=H2M.make_params(LyA01=dict(value=ig[0],vary=False),
                        LyA02=dict(value=ig[1],vary=False),
                        LyA14=dict(value=ig[2],vary=False),
                        LyA17=dict(value=ig[3],vary=False),
                        z=dict(value=ig[4],min=pb[4][0],max=pb[4][1]),
                        gamma=dict(value=ig[5],min=pb[5][0],max=pb[5][1]),
                        T=dict(value=ig[6],min=pb[6][0],max=pb[6][1]),
                        q=dict(value=ig[7],min=pb[7][0],max=pb[7][1]),
                        rchar=dict(value=ig[8],min=pb[8][0],max=pb[8][1]),
                        MH2=dict(value=ig[9],min=pb[9][0],max=pb[9][1]))
    
    #fitting results and other useful things
    totResult=H2M.fit(fluxD,H2P,vradData=vradD,weights=1.0/errrD,scale_covar=False)
    totBestfit=totResult.best_fit
    totReport=totResult.fit_report()
    totSqr=totResult.chisqr
    totRed=totResult.redchi
    totChi=[totSqr,totRed]
    totBic=totResult.bic
    
    return totBestfit,totResult.params,totReport,totChi,totBic,H2M

def extractParam(paramList,paramName,paramPrnt=False):
    '''
    extract, and optionally print, all parameters
    '''
    allVal,allErr=[],[]
    for name in paramName:
        val=paramList[name].value
        err=paramList[name].stderr
        
        if paramPrnt:
            print(f'{name: <5} : {val:8.3f} +/- {err:.3f}')
            
        allVal.append(val)
        allErr.append(err)

    return allVal,allErr

def plotPower(array):
    pltpow=np.floor(np.log10(np.nanmax(array)))
    pltdiv=10**(pltpow)
    return int(pltpow),pltdiv

def makeFlatchain(combArr,keys):
    arrShape=combArr.shape
    dim1=arrShape[0]*arrShape[1]
    dim2=arrShape[2]
    flat=combArr.reshape(dim1,dim2)
    flatDF=pd.DataFrame(flat,columns=keys)
    return flatDF #return the flatchain in the same format as created by emcee

def combineChain(chains,addls,keys):
    fullSteps=0
    if len(chains)==1:
        print('Only 1 chain provided!')
        return chains[0],makeFlatchain(chains[0],keys),addls[0]['steps']
    for i in range(1,len(chains)):
        if addls[i-1]['thin']==addls[i]['thin']: #thinning parameters must match
            if i==1:
                #concatenate the steps (axis=0), keeping walkers and variables constant
                fullChain=np.concatenate([chains[i-1],chains[i]],axis=0)
                fullSteps+=addls[i-1]['steps']
                fullSteps+=addls[i]['steps']
            else:
                fullChain=np.concatenate([fullChain,chains[i]],axis=0)
                fullSteps+=addls[i]['steps']
            flatChain=makeFlatchain(fullChain,keys)
        else:
            print('Thinning Parameters Do Not Match!')
            return None,None,None
    return fullChain,flatChain,fullSteps

def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i

def autocorr_func_1d(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf

def auto_window(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1

def autocorr_new(y, c=5.0):
    f = np.zeros(y.shape[1])
    for yy in y:
        f += autocorr_func_1d(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]



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
targ = "EP Cha"
targ_d = 98.816 
targ_AV = 0.0 
targ_M = 0.8 
targ_inclination = 70.0
EBV = -1.0 * targ_AV / Rv #reddening/attenuation

#read in the disk data
preppedData='EP_Cha_COS_ULLYSES_H2_Data.npz'
ref_waves,obs_vel,progVals,lyawVals,target_velocity,target_fluxes,target_fluxerr,lines_to_fit_dict=loadData(preppedData)
model_refwaves=ref_waves #check if deepcopy is needed

#observation details (from ULLYSES)
grtng=[] #The ULLYSES data switches to G160M data at 1470A
cwave=[] #G130M combined 1327 and 1291 data, G160M combined 1577, 1600, and 1623 data
lifep=[] #G130M was all LP1, G160M was also all LP1
#[1,4]: G130M, G160M, G130M, G160M
#[1,7]: G130M, G160M, G160M, G160M
#[0,1]: G130M, G130M


#Need to shift lines according to measured velocity, done for all data regardless of whether or not we'll actually fit the line
abs_vel=100 #cutoff between line and continuum
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

#make an initial guess for each model parameter, and define parameter bounds for the model
paramName=['LyA01',      'LyA02',      'LyA14',      'LyA17',      'z',      'gamma',   'T',           'q',       'rchar',   'MH2'] #parameter names (for printing mainly)
# initGuess=[-11.23,       -11.50,       -10.97,        -11.3,       6.0,      0.5,       1000.0,        -1.25,     25.0,      -4.0] #initial parameter guesses
# paramBnds=[[-14.0,-10.0],[-14.0,-10.0],[-14.0,-10.0],[-14.0,-10.0],[2.0,7.0],[0.0,1.99],[500.0,5000.0],[-2.5,2.5],[0.0,50.0],[-5.0,-1.0]] #parameter bounds
initGuess=[-11.23,       -11.50,       -10.97,        -11.3,       6.0,      0.5,       1000.0,        -1.25,     25.0,      -10.0] #initial parameter guesses
paramBnds=[[-14.0,-10.0],[-14.0,-10.0],[-14.0,-10.0],[-14.0,-10.0],[2.0,7.0],[0.0,1.99],[500.0,5000.0],[-2.5,2.5],[0.0,50.0],[-16.0,-9.0]] #parameter bounds

#initialize the PPD grid and associated properties
grid_dict = grid_from_radius(targ_d, targ_inclination) #contains the cartesian (x,y) and cylindrical (r,phi) coordinates for each disk grid point (r: 201, phi:180)
vobs = vrad_to_vobs(grid_dict["rgrid"], grid_dict["phigrid"], targ_inclination) #the radial velocity at each grid point
vobs_flat = vobs.flatten() #flatten disk RV into a 1D array of len 201*180 
abs_vel_mlt=3 #multiplier for the final profile length (more lines = more profile length)
bins, edges = np.histogram(vobs_flat, bins=599, range=(-abs_vel_mlt*abs_vel, abs_vel_mlt*abs_vel)) #create a histogram of disk RV values
inds = np.digitize(vobs_flat, edges) #np.digitize creates an array of which bin each value falls within
v_obs_bin = np.arange(-abs_vel_mlt*abs_vel, abs_vel_mlt*abs_vel, 1.0) #bin values (cleaner than edges --> whole numbers) 

#run the model
print('Finding the Best Fit with LMFIT')
start=time.time()
bestFit,bestPar,mReport,bestChi,bestBic,modlObj=fit_H2Model(all_datavel,all_dataflux,all_dataerr,initGuess,paramBnds)
ended=time.time()
print(f'Time to Best Fit: {ended-start:.1f} seconds')

parVal,parErr=extractParam(bestPar,paramName,True)

ppow,pscl=plotPower(all_dataflux)
fig,ax=plt.subplots(figsize=(11.5,8),layout='constrained')
plt.title(f"EP Cha: Best Fit ($\\chi^2_R$ = {bestChi[1]:.3f})",fontsize=20)
plt.plot(np.arange(len(all_dataflux)),all_dataflux/pscl,color="k",lw=2,label='Data')
plt.fill_between(np.arange(len(all_dataflux)),(all_dataflux+all_dataerr)/pscl,(all_dataflux-all_dataerr)/pscl,alpha=0.35,color='k')
plt.plot(np.arange(len(all_dataflux)),bestFit/pscl,'-o',color="teal", lw=2,label='Model')
plt.xlabel('Data Points (To Be Changed)',fontsize=18)
plt.ylabel('Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',fontsize=18)
plt.legend(fontsize=18)
plt.gca().minorticks_on()
plt.grid(alpha=0.35)
plt.tick_params(labelsize=18)
plt.show()

combFlux,indiFlux,lineProf,lsfX,lsfY,paramD,diskModel=get_H2model_all(all_datavel,*parVal)

fig,ax=plt.subplots()
for x,y in zip(lsfX,lsfY):
    plt.plot(x,y)
plt.show()

preConvLines=diskModel.all_emission_line_intensities()
fig,ax=plt.subplots()
for i in range(preConvLines.shape[1]):
    plt.plot(preConvLines[:,i])
plt.show()

# fig,ax=plt.subplots(figsize=(8,6),layout='constrained')
# plt.semilogx(diskModel.grid_dict["rgrid"],np.sum(lineProf,axis=(0,2)),color='teal',linewidth=2) #!!!why are going out to 100,000 AU???
# plt.xlabel('Disk Radius (AU)',fontsize=18)
# plt.ylabel('F$_{H_2}$(r) / F$_{H_2}$(max)',fontsize=18)
# plt.title('EP Cha Radial Flux Distribution',fontsize=18)
# plt.tick_params(labelsize=18)
# plt.grid(alpha=0.35)
# plt.xlim(1e-1,1e1)
# plt.show()
fig,ax=plt.subplots(figsize=(8,6),layout='constrained')
plt.semilogx(diskModel.grid_dict["rgrid"],np.sum(lineProf,axis=(0,2))/np.max(np.sum(lineProf,axis=(0,2))),color='teal',linewidth=2) 
plt.xlabel('Disk Radius (AU)',fontsize=18)
plt.ylabel('F$_{H_2}$(r) / F$_{H_2}$(max)',fontsize=18)
plt.title('EP Cha Radial Flux Distribution',fontsize=18)
plt.tick_params(labelsize=18)
plt.grid(alpha=0.4)
plt.grid(which='minor',axis='x',alpha=0.2)
plt.xlim(1e-1,1e1)
plt.show()

#following code written by Rachel Pauline, slightly modified by me
rpColors=['cornflowerblue', 'orange', 'green', 'firebrick', 'mediumorchid', 'orangered', 'pink', 'gray', 'olive', 'cyan', 'deeppink', 'navy', 'khaki', 'lightsalmon', 'lime'] #colors from Rachel

ppow,pscl=plotPower(all_dataflux)
fig,ax=plt.subplots(figsize=(11.5,8),layout='constrained')
for i in range(len(target_velocity)):
    ax.plot(target_velocityL[i],target_fluxesL[i]/pscl,'-',markersize=2,color=rpColors[i],label=f'{ref_waves[i]} $\\AA$ - {progVals[i]}')
plt.xlabel('Velocity ($km$ $s^{-1}$)',fontsize=18)
plt.ylabel('Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',fontsize=18)
plt.title('H2 Emission Lines (Corrected for Stellar RV)',fontsize=18)
plt.tick_params(labelsize=18)
plt.legend(loc='upper right',fontsize=16)
plt.grid(alpha=0.35)
plt.show()

line_lengths=[0] #lengths of individual emission lines
for f in target_fluxesL:
    line_lengths.append(line_lengths[-1]+len(f))

ppow,pscl=plotPower(all_dataflux)
fig,ax=plt.subplots(4,3,figsize=(14,10),layout='constrained') #!!!change subplot size/shape as needed
axs=ax.flatten()
for i in range(len(target_velocityL)):
    axs[i].plot(target_velocityL[i],target_fluxesL[i]/pscl,color='k',linewidth=2)
    axs[i].plot(target_velocityL[i],bestFit[line_lengths[i]:line_lengths[i+1]]/pscl,color='teal',linewidth=2)
    axs[i].set_title(f'{ref_waves[i]} $\\AA$ - {progVals[i]}',fontsize=18)
    axs[i].tick_params(labelsize=18)
    axs[i].grid(alpha=0.35)
for j in range(i+1,len(axs)):
    #turn off remaining subplots
    axs[j].set_visible(False)
fig.supxlabel('Velocity ($km$ $s^{-1}$)',fontsize=18)
fig.supylabel('Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',fontsize=18)
axs[-1].plot([],[], color='k', label='Data',linewidth=2)
axs[-1].plot([],[], color='teal',label='LMFIT',linewidth=2)
fig.legend(fontsize=18,loc='upper left',bbox_to_anchor=(0.8,0.2)) #anchor the upper left point
plt.show()

#setup mcmc parameters
nWalkers=80 #how many walkers will roam the parameter space
nSteps=250 #how many samples are drawn by each walker, how many "steps" they take
nBurn=20 #how many of the first samples/steps are discarded for each walker (set to 0 if continuing a previous chain)
nThin=10 #thin tells it to accept 1 in every x samples (reduces data output)
#can't have less than 2*dimension walkers
#try to go for an integer value of (nSteps-nBurn)/nThin 
#nW/nS/nB/nT predicted to take X time units, actually took Y time units

#run the mcmc starting at the best fit location
emcee_kws = dict(nwalkers=nWalkers,steps=nSteps,burn=nBurn,thin=nThin,is_weighted=True,progress=True)
emcee_params = bestPar.copy() #copies the final parametrs from the previous lmfit result
#!!!copypasta begins here
newChain=False
#newChain=True 
if newChain:
    print(f'Starting MCMC with {nWalkers} Walkers, {nSteps} Steps, a Burn of {nBurn}, and a Thinning Factor of {nThin}')
    startM=time.time()
    result_emcee=modlObj.fit(all_dataflux,emcee_params,vradData=all_datavel,weights=1.0/all_dataerr,method='emcee',fit_kws=emcee_kws)
    endedM=time.time()
    print(f'MCMC Runtime: {(endedM-startM)/60:.1f} minutes')
    print(result_emcee.fit_report()) 
    #the results are the median parameter values +/- the 1sigma values about the median
    #+/-1sigma is half the difference of the 15.87 and the 84.13 percentiles
    #the chain attribute contains the saved samples, has the following shape:
    #((steps-burn)//thin, nwalkers, nvarys) where nvarys=the number of free variables
    #the flatchain attribute flattens the chain for each free variable
    #the lnprob attribute has the log(probability) for each chain sample, highest prob = maximum likelihood estimate
    #the acor attribute has the computed autocorrelation times if it could be computed from the chain
    #the acceptance_fraction attribute is an array of the acceptance fractions of accepted steps for each walker
    
    #plot the acceptance fraction
    fig,ax=plt.subplots()
    plt.plot(result_emcee.acceptance_fraction)
    plt.axhline(0.2,linestyle='--',color='r')
    plt.axhline(0.5,linestyle='--',color='r')
    plt.xlabel('walker')
    plt.ylabel('acceptance fraction')
    plt.show()
    #acceptance fraction should be between ~0.2 and ~0.5 after convergence
    #this is the fraction of steps where the walker did NOT move back to its previous step
    
    #make a corner plot of the parameter distributions
    varMask=[True if emcee_params[key].vary else False for key in emcee_params.keys()] #determine which are variables
    paramVar=np.array(paramName)[varMask] #get names of variable parameters
    emcee_corner=corner.corner(az.convert_to_inference_data(result_emcee.chain),
                               labels=result_emcee.var_names,
                               truths=np.array(list(result_emcee.params.valuesdict().values()))[varMask]) #remove the first two fixed params
    
    ppow,pscl=plotPower(all_dataflux)
    fig,ax=plt.subplots(2,1,sharex=True)
    plt.suptitle('EP Cha MCMC Results',fontsize=24)
    fig.text(0.04,0.5,'Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',va='center',rotation='vertical',fontsize=18)  
    plt.sca(ax[0])
    plt.plot(np.arange(len(all_dataflux)),all_dataflux/pscl,label='Data',color='k')
    plt.fill_between(np.arange(len(all_dataflux)),(all_dataflux+all_dataerr)/pscl,(all_dataflux-all_dataerr)/pscl,alpha=0.35,color='k')
    plt.plot(np.arange(len(all_dataflux)),bestFit/pscl,label='LMFIT',color="teal")
    plt.plot(np.arange(len(all_dataflux)),result_emcee.best_fit/pscl,label='MCMC',linestyle='--',color='darkgreen')
    plt.legend(fontsize=18)
    plt.tick_params(labelsize=18)
    OoM=ax[0].yaxis.get_offset_text()
    OoM.set_size(18)
    plt.sca(ax[1])
    plt.plot(np.arange(len(all_dataflux)),(bestFit-result_emcee.best_fit)/pscl,label='LMFIT - MCMC',color='indigo')
    plt.axhline(0.0,color='k',linestyle='--')
    plt.xlabel('Data Points (To Be Changed)',fontsize=18)
    plt.legend(fontsize=18)
    plt.tick_params(labelsize=18)
    OoM=ax[1].yaxis.get_offset_text()
    OoM.set_size(18)
    plt.show()
    
    #find and print the maximum likelihood estimates compared to the median estimates
    highest_prob=np.argmax(result_emcee.lnprob) #find the highest probability for each variable
    hp_loc=np.unravel_index(highest_prob,result_emcee.lnprob.shape) #find the locaitons of the highest probs
    mle_soln=result_emcee.chain[hp_loc] #get the maximum likelihood variable values
    mle_params=bestPar.copy() #copy the parameter structure
    
    track=0 #track how many variables have been added 
    for i, par in enumerate(mle_params):
        if par in paramVar:
            mle_params[par].value = mle_soln[track] #save the mle params
            track+=1 #!!!now that I have varMAsk, maybe I can do this loop more efficiently?
    
    print('\nMaximum Likelihood Estimation from emcee       ')
    print('-------------------------------------------------')
    print('Parameter  MLE Value   Median Value   Uncertainty')
    fmt = '  {:5s}  {:11.5f} {:11.5f}   {:11.5f}'.format
    for name, param in mle_params.items():
        if name in paramVar:
            print(fmt(name, param.value, result_emcee.params[name].value,
                      result_emcee.params[name].stderr))
            
    #get the actual 1sigma and 2sigma values, not just the percentile difference cut in half
    print('\nError estimates from emcee:')
    print('------------------------------------------------------')
    print('Parameter    -2sigma    -1sigma     median    +1sigma    +2sigma')
    for name in result_emcee.params.keys():
        if name in paramVar:
            quantiles = np.percentile(result_emcee.flatchain[name],
                                      [2.275, 15.865, 50, 84.135, 97.275])
            median = quantiles[2]
            err_m2 = quantiles[0] - median
            err_m1 = quantiles[1] - median
            err_p1 = quantiles[3] - median
            err_p2 = quantiles[4] - median
            fmt = '  {:5s}   {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}'.format
            print(fmt(name, err_m2, err_m1, median, err_p1, err_p2))
            
    #print the autocorrelation times if they were able to be computed
    if hasattr(result_emcee, "acor"):
        print("\nAutocorrelation time for the parameters:")
        print("----------------------------------------")
        for i in range(len(paramVar)):
            print(f'{paramVar[i]: <5}: {result_emcee.acor[i]:.2f} ({int(result_emcee.acor[i]*50)} steps')
    else:
        print("\nUnable to determine autocorrelation time for the parameters from emcee. Estimating manually from the chain.")
        for i in range(len(paramVar)):
            estimate=nThin*autocorr_new(np.transpose(result_emcee.chain[:,:,i]))
            print(f'{paramVar[i]: <5}: {estimate:.2f} (Estimate, {int(estimate*50)} steps)')
    #autocorrelation time is not really a time, it is a convergence rate
    #AT is ~the number of markov chain transitions equivalent to a single independent draw of the distribution
    #mcmc samples aren't independent, a number of samples needs to be drawn to effectively get 1 independent sample
    #AT is the number of samples where the mcmc "forgets" where it started
    #N/AT is the effective number of samples
    #should at minimum do 50*autocorrelation number of effective samples
    
    #save=True
    save=False
    if save:
        try:
            important=[result_emcee.fit_report(),result_emcee.acceptance_fraction,result_emcee.flatchain,
                       result_emcee.var_names,result_emcee.params,result_emcee.best_fit,
                       result_emcee.lnprob,result_emcee.chain,emcee_kws,result_emcee.acor]
        except(AttributeError):
            important=[result_emcee.fit_report(),result_emcee.acceptance_fraction,result_emcee.flatchain,
                       result_emcee.var_names,result_emcee.params,result_emcee.best_fit,
                       result_emcee.lnprob,result_emcee.chain,emcee_kws] #if acor attribute doesn't exist
        file_mcmc=open('EP_Cha_MCMC_Results_and_Params.obj','wb')
        pickle.dump(important,file_mcmc)
        file_mcmc.close()
        #save the MCMC result data so that it can be loaded in later
        
contMCMC=False
#contMCMC=True #!!!make sure you are loading the right files in
if contMCMC:
    #to continue the mcmc from where it left off, open the data file which has the last saved positions of the walkers
    runFiles=glob('EP_Cha_MCMC*.obj') 
    runDat,runAux,runKey=[],[],[]
    for f in runFiles:
        with open(f,'rb') as file:
            loadData=pickle.load(file) #reloads the list of 9 (or 10) objects in the 'important' list defined above
            runDat.append(loadData[7]) #grab and save the relevant data for continuing the chain
            runAux.append(loadData[8])
            runKey.append(loadData[2].keys())
    #combine the chains and estimate autocorrelation time (more steps = more accuracy)
    combined,combflat,allsteps=combineChain(runDat,runAux,runKey)     
    chainEst=0
    if isinstance(combined,np.ndarray):
        # Compute the estimators for autocorrelation time
        T=runAux[0]['thin'] #the thinning parameter (already checked that this is the same for both)
        new=np.zeros(combined.shape[2])
        for i in range(0,combined.shape[2]):
            new[i]=T*autocorr_new(np.transpose(combined[:,:,i]))
            print(f'ACT Estimates for {paramVar[i]}: {new[i]:.2f} ({int(new[i]*50)} Steps Needed)')
            if allsteps<new[i]*50:
                print(f'Chain Length of {allsteps} is Shorter than 50*ACT, Run a Longer Chain!\n')
                if new[i]*50>chainEst:
                    chainEst=new[i]*50 #save the largest estimate for total chain length
    if allsteps<chainEst:
        print(f'Run the Chain for {int(chainEst-allsteps)} More Steps')
        
    #make combined chain output plots
    cmbtruth=[np.percentile(combflat[name],50) for name in paramVar]
    emcee_corner=corner.corner(az.convert_to_inference_data(combined),
                               labels=paramVar,truths=cmbtruth)
    cmbmodel=get_H2model(all_datavel,initGuess[0],initGuess[1],*cmbtruth)
    ppow,pscl=plotPower(all_dataflux)
    fig,ax=plt.subplots(2,1,sharex=True)
    plt.suptitle('EP Cha MCMC Results',fontsize=24)
    fig.text(0.04,0.5,'Flux Density ($\\times$10$^{'+str(ppow)+'}$ $erg$ $cm^{-2}$ $s^{-1}$ $\\AA^{-1}$)',va='center',rotation='vertical',fontsize=18)  
    plt.sca(ax[0])
    plt.plot(np.arange(len(all_dataflux)),all_dataflux/pscl,label='Data',color='k')
    plt.fill_between(np.arange(len(all_dataflux)),(all_dataflux+all_dataerr)/pscl,(all_dataflux-all_dataerr)/pscl,alpha=0.35,color='k')
    plt.plot(np.arange(len(all_dataflux)),bestFit/pscl,label='LMFIT',color="teal")
    plt.plot(np.arange(len(all_dataflux)),cmbmodel/pscl,label='MCMC',linestyle='--',color='darkgreen')
    plt.legend(fontsize=18)
    plt.tick_params(labelsize=18)
    OoM=ax[0].yaxis.get_offset_text()
    OoM.set_size(18)
    plt.sca(ax[1])
    plt.plot(np.arange(len(all_dataflux)),(bestFit-cmbmodel)/pscl,label='LMFIT - MCMC',color='indigo')
    plt.axhline(0.0,color='k',linestyle='--')
    plt.xlabel('Wavelength ($\\AA$)',fontsize=18)
    plt.legend(fontsize=18)
    plt.tick_params(labelsize=18)
    OoM=ax[1].yaxis.get_offset_text()
    OoM.set_size(18)
    plt.show()
    
    
    #print useful combined chain outputs
    print('\nError estimates from emcee:')
    print('------------------------------------------------------')
    print('Parameter    -2sigma    -1sigma     median    +1sigma    +2sigma')
    for name in paramName:
        if name in paramVar:
            quantiles = np.percentile(combflat[name],[2.275, 15.865, 50, 84.135, 97.275])
            median = quantiles[2]
            err_m2 = quantiles[0] - median
            err_m1 = quantiles[1] - median
            err_p1 = quantiles[3] - median
            err_p2 = quantiles[4] - median
            fmt = '  {:5s}   {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}'.format
            print(fmt(name, err_m2, err_m1, median, err_p1, err_p2))
    
    #set the parameters and and allow the walkers to keep roaming
    nSteps=2000 #how many additional steps to go for
    nBurn=0 #no need to burn in when continuing the chain
    emcee_kws=dict(nwalkers=nWalkers,steps=nSteps,burn=nBurn,thin=nThin,pos=combined[-1],is_weighted=True,progress=True) #continue from most recent walker positions
    #copypasta code above to let it rip
    