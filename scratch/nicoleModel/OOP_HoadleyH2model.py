"""Python version of Keri Hoadley's H2 fluorescence modeling code
Converted from IDL by Nicole Arulanantham
Last updated April 2019
"""
import numpy as np
import pandas as pd
import math
import itertools
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow
#Import script with Python version of IDL's ccm_unred
from ccm_unred import which_spec_range, ccm_unred
from itertools import product
#Import script that generates wavelength-dependent line-spread functions for COS
from linefittingfunctions import cos_lsf_arulanantham

#Define global constants
#Gravitational constant [cm^3 g^-1 s^-2]
G = 6.67259e-8
#Speed of light [cm/s]
c = 3.e10
#Planck constant [erg s]
h = 6.626e-27
#Doppler parameter [km/s]
b_km = 5.
#Electron charge [esu]
e = 4.803206814e-10
#Electron mass [g]
m_e = 9.109389754e-28
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

class CircumstellarDisk(object):
    """Initialize a disk around a star, based on stellar properties and viewing angle
    Should be able to do radiative transfer modeling (i.e. of H2 or CO) on this physical structure
    """
    #Doppler velocity of the target [km/s]
    doppler = 0.0
    #Turbulent velocity of gas in disk [km/s]
    turbulence = 0.5
    #Radii limits in disk (r_in -> r_mid will have higher resolution)
    #r_mid -> r_out will have lower resolution
    r_in = 0.01
    r_mid = 1.
    r_out = 50.
    #Keri kept r_LyA constant - this drives where LyA photons originate in the disk/accretion front
    r_LyA = 0.01

    def __init__(self, Mstar, Av, d, param_dict):
        #Stellar mass [Msun]
        self.Mstar = Mstar
        #Visual extinction [mag]
        self.Av = Av
        #Distance [pc]
        self.d = d
        #Variable disk parameters
        #These can be fit in modeling procedure
        self.param_dict = param_dict

    def calc_phi_H2(self, xval, yval):
        """Helper method for phi_H2 calculation in grid_from_radius()
        """
        if xval <= 0.:
            return math.atan(yval / xval) + math.pi
        elif yval <= 0. and xval > 0.:
            return math.atan(yval / xval) + (2. * math.pi)
        else:
            return math.atan(yval / xval)

    def grid_from_radius(self):
        """Builds grids to describe disk based on stellar distance and disk inclination
        """
        start_time = time.time()

        rinout = 10.**(0.032 * np.arange(0., 201., 1.) - 1.4)
        r_grid = rinout / self.d

        phi_grid = np.arange(0., 360., 2.) * math.pi / 180.
        tan_phi = np.tan(phi_grid)

        xgrid = np.array([-1. * np.sqrt(np.square(r_grid) / (1. + az ** 2.)) / math.cos(self.param_dict["Inclination"]) if math.pi / 2. <= az <= 3. * math.pi / 2.
                          else np.sqrt(r_grid ** 2. / (1. + az ** 2.)) / math.cos(self.param_dict["Inclination"]) for az in tan_phi]).astype(np.float64)
        xgrid = np.swapaxes(xgrid, 0, 1)

        ygrid = np.array([-1. * np.sqrt(r_grid ** 2. / (1. + (1. / az ** 2.))) if phi_grid[az_idx] >= math.pi
                          else np.sqrt(r_grid ** 2. / (1. + (1. / az ** 2.))) for az_idx, az in enumerate(tan_phi)]).astype(np.float64)
        ygrid = np.swapaxes(ygrid, 0, 1)

        r_H2 = np.array([np.sqrt(np.square(xgrid[:, az_idx]) + np.square(ygrid[:, az_idx])) * self.d for az_idx in range(len(tan_phi))]).astype(np.float64)
        r_H2 = np.swapaxes(r_H2, 0, 1)

        phi_H2 = np.array([[self.calc_phi_H2(xgrid[r_idx, az_idx], ygrid[r_idx, az_idx]) for r_idx, r in enumerate(r_grid)] for az_idx, az in enumerate(tan_phi)]).astype(np.float64)
        phi_H2 = np.swapaxes(phi_H2, 0, 1)

        return {"xgrid": xgrid,
                "ygrid": ygrid,
                "rgrid": rinout, #np.array([0.1, 0.2, 0.3, 0.4, 0.5, 1., 2., 3., 5., 10.]).astype(np.float64)
                "phigrid": phi_grid * 180. / math.pi,
                "rH2": r_H2,
                "phiH2": phi_H2 * 180. / math.pi} #"rgrid": r_grid * self.d

    def sightline_x(self, r_array, phi):
        """Helper method for sightline_r_phi_z()
        """
        return (r_array / self.d) * math.cos(phi * math.pi / 180.)

    def sightline_y(self, r_array, phi):
        """Helper method for sightline_r_phi_z()
        """
        return -1. * (r_array / self.d) * math.sin(phi * math.pi / 180.) - self.d * math.sin(self.param_dict["Inclination"])

    def sightline_r_phi_z(self, r_grid, phi_grid):
        """Done checking against IDL version!
        """
        start_time = time.time()
        sz = -1. * self.d * math.cos(self.param_dict["Inclination"])

        sightline = np.array([np.sqrt(np.square(sz) + np.square(self.sightline_x(r_grid, phi)) + np.square(self.sightline_y(r_grid, phi)))
                      for phi in phi_grid]).astype(np.float64)
        sightline = np.swapaxes(sightline, 0, 1)

        return sightline

    def radial_velocity(self, r_grid):
        """Done checking against IDL version!
        """
        v_k = np.sqrt(G * self.Mstar * MSUNtoG / (r_grid * AUtoCM)) * CMtoKM

        return v_k

    def vrad_to_vobs(self, r_grid, phi_grid):
        """Done checking against IDL version!
        """
        v_k = self.radial_velocity(r_grid)
        vobs = np.array([v * math.sin(self.param_dict["Inclination"]) * np.sin(phi_grid * math.pi / 180.) for v in v_k]).astype(np.float64)

        return vobs

    def radial_T(self, r_grid):
        """Done checking against IDL version!
        """
        if self.param_dict["Molecule"] == "H2":
            T_fixed_loc = 1. #[au]
            T_grid = self.param_dict["T"]*np.power(r_grid / T_fixed_loc, -1.*self.param_dict["q"])
            T_H2 = [T if T <= 5000. and T >= 500. else 0. for T in T_grid]
        elif self.param_dict["Molecule"] == "CO":
            T_fixed_loc = 10. #[au]
            T_grid = self.param_dict["T"]*np.power(r_grid / T_fixed_loc, -1.*self.param_dict["q"])
            T_H2 = [T if T <= 1800. and T >= 20. else 0. for T in T_grid]

        return np.array(T_H2).astype(np.float64)

    def scale_height(self):
        """Done checking against IDL version!
        """
        r_grid = self.grid_dict["rgrid"]
        T_grid = self.T

        if self.param_dict["Molecule"] == "H2":
            mmw = 2.33 * 1.67e-24 #mean molecular weight for circumstellar material (Ruden & Pollack 1991)
        elif self.param_dict["Molecule"] == "CO":
            mmw = 28. * 1.66054e-24 #[g]

        thermal = kB * T_grid / mmw
        kinetic = (np.power(r_grid * AUtoCM, 3.) / (G * self.Mstar * MSUNtoG))

        if len(np.shape(T_grid)) > 1:
            H_z = [np.sqrt(thermal[:, z_idx] * kinetic) * r_grid for z_idx in range(np.shape(T_grid)[1])]
            H_z = np.array(H_z).astype(np.float64)
        else:
            H_z = np.sqrt(thermal * kinetic) / AUtoCM #[au]

        return H_z

    def exp_disk_density(self):
        """Done checking against IDL version!
        """
        z_array = np.arange(0., 11.)
        z_step = 11.
        # z_array = np.arange(0., 11.)
        # z_step = 11.
        H_z = self.scale_height()

        #rho_disk = np.zeros((len(H_z), len(z_array)))
        rho_disk = np.zeros(len(z_array))

        for z_idx, z in enumerate(z_array):
            #rho_disk[:, z_idx] = np.exp(-0.5 * np.square(self.param_dict["z"] - z) / np.square(H_z)) #
            rho_disk[z_idx] = np.exp(-0.5 * np.square(self.param_dict["z"] - (z / z_step)))

        return rho_disk

    def surface_density(self):
        """Done checking against IDL version!
        """
        r_grid = self.grid_dict["rgrid"]
        power_law = np.power(r_grid / self.param_dict["rchar"], -1. * self.param_dict["gamma"])
        exp_law = np.exp(-1. * np.power(r_grid / self.param_dict["rchar"], 2. - self.param_dict["gamma"]))

        return power_law * exp_law

    def angular_grid(self, dist_grid):

        rolled_in_r = np.roll(dist_grid, 1)
        ang_grid = np.absolute(dist_grid - rolled_in_r)
        # ang_grid[0] = ang_grid[1]

        return ang_grid

    def collapsed_line_profile(self, lineprof):
        """Integrate disk fluxes in velocity space to produce emission line profiles
        """
        vobs = self.vrad_to_vobs(self.grid_dict["rgrid"], self.grid_dict["phigrid"])
        vobs_flat = vobs.flatten() #same as np.ravel with order 'C'
        bins, edges = np.histogram(vobs_flat, bins=599, range=(-300., 300.))
        inds = np.digitize(vobs_flat, edges)
        v_obs_bin = np.arange(-300., 300., 1.)
        binmin = min(v_obs_bin)
        binmax = max(v_obs_bin)

        total_line_profile = np.zeros((len(v_obs_bin), len(lineprof[:, 0, 0])), dtype=np.float64)
        total_line_profile.fill(1.e-75)

        for prog_idx in range(0, len(lineprof[:, 0, 0])):
            flat_fluxes = lineprof[prog_idx, :, :].flatten()
            test_df = pd.DataFrame({"Velocities": inds,
                                    "Fluxes": flat_fluxes})

            g = test_df.groupby(["Velocities"])
            profile = np.array(g.sum()).astype(np.float64).flatten()

            pads = int((600 - len(profile)) / 2)
            if pads < 600:
                pad_left = pads+1
                pad_right = pads
            elif pads > 600:
                pad_left = pads-1
                pad_right = pads

            final_profile = np.pad(profile, (pad_left, pad_right), 'edge')

            total_line_profile[:, prog_idx] = final_profile
            total_line_profile[300, prog_idx] = total_line_profile[301, prog_idx]

        return total_line_profile

    def __del__(self):
        # print "Deleted a disk!"
        del self

class H2Model(CircumstellarDisk):
    #States to model
    J_initial = np.array([0, 1, 5, 6]).astype(np.float64) #J_ROT_STATE
    g_initial = np.array([2 * J + 1 if J % 2 == 0 else 3 * (2 * J + 1) for J in J_initial]).astype(np.float64)
    J_progression = np.array([1, 2, 4, 7]).astype(np.float64)
    g_progression = np.array([2 * J + 1 if J % 2 == 0 else 3 * (2 * J + 1) for J in J_progression]).astype(np.float64)
    lower_energy_eV = np.array([1.00, 1.02, 1.20, 1.27]).astype(np.float64) #From Table 2 of Herczeg et al. 2004
    #energy of the lower state of the transition, in ergs
    lower_energy = np.array([val*1.60217657e-12 for val in lower_energy_eV]).astype(np.float64)
    # [0,1], [0,2], [1,4], [1,7]
    LyA_pumping = np.array([1217.205, 1217.643, 1216.070, 1215.726]).astype(np.float64)
    # total Aul for each progression (from Abgrall et al. 1993)
    Aul_total = np.array([1.8e9, 1.8e9, 1.7e9, 1.7e9]).astype(np.float64)
    #ground vibration state, before LyA pumping
    v1 = 2
    #Disk filling fraction
    eta = 0.25
    #Mean molecular weight of H2
    mmw_H2 = 2.33 * 1.67e-24
    lambda_props = pd.DataFrame({"Wavelength": np.array([1338.56, 1342.26, 1393.96, 1398.95, 1402.65, 1446.12, 1460.17, 1463.83, 1467.08, 1489.57, 1500.45, 1504.76, 1521.59, 1524.65, 1580.67]).astype(np.float64),
                                 "Aul": np.array([3.1E8, 2.8E8, 1.6E8, 2.6E8, 2.3E8, 1.4E8, 1.5E8, 1.4E8, 1.3E8, 1.6E8, 1.7E8, 2.0E8, 6.0E8, 1.9E8, 1.1E8]).astype(np.float64),
                                 "Jl": np.array([0, 1, 1, 0, 1, 5, 0, 1, 6, 5, 6, 5, 0, 6, 6]).astype(np.float64),
                                 "Ju": np.array([1, 2, 2, 1, 2, 4, 1, 2, 7, 4, 7, 4, 1, 7, 7]).astype(np.float64),
                                 "Jf": np.array([2, 3, 1, 2, 3, 5, 2, 3, 8, 3, 6, 5, 2, 8, 8]).astype(np.float64),
                                 "Aul_total": np.array([1.8E9, 1.8E9, 1.8E9, 1.8E9, 1.8E9, 1.7E9, 1.8E9, 1.8E9, 1.7E9, 1.7E9, 1.7E9, 1.7E9, 1.8E9, 1.7E9, 1.7E9]).astype(np.float64)
                                 })
    lambda_props["Bul"] = lambda_props["Aul"] / lambda_props["Aul_total"]
    lambda_props["gl"] = np.array([2. * Jl + 1 if Jl % 2 == 0 else 3. * (2. * Jl + 1.) for Jl in lambda_props["Jl"]]).astype(np.float64)
    lambda_props["glp1"] = np.array([2. * Ju + 1 if Ju % 2 == 0 else 3. * (2. * Ju + 1.) for Ju in lambda_props["Ju"]]).astype(np.float64)

    Z_df = pd.read_csv("H2_PPD_modeling_code/disk_partition_values.csv")
    Z_TH2 = []
    for key in list(Z_df):
        Z_TH2.append(np.array(Z_df[key]).astype(np.float64))
    Z_TH2 = np.array(Z_TH2).astype(np.float64)

    column_tau_df = pd.read_csv("H2_PPD_modeling_code/column_tau.csv")
    column_tau = []
    for key in list(column_tau_df):
        column_tau.append(np.array(column_tau_df[key]).astype(np.float64))
    column_tau = np.array(column_tau).astype(np.float64)

    T_tau_df = pd.read_csv("H2_PPD_modeling_code/T_tau.csv")
    T_tau = []
    for key in list(T_tau_df):
        T_tau.append(np.array(T_tau_df[key]).astype(np.float64))
    T_tau = np.array(T_tau).astype(np.float64)

    eff_tau_df = pd.read_csv("H2_PPD_modeling_code/eff_tau.csv")
    eff_tau = []
    for key in list(eff_tau_df):
        eff_tau.append(np.array(eff_tau_df[key]).astype(np.float64))
    eff_tau = np.array(eff_tau).astype(np.float64)

    def __init__(self, Mstar, Av, d, param_dict,
                 DALI_rgrid=np.array([99.999]),
                 DALI_Tgrid=np.array([99.999]),
                 DALI_ncol=np.array([[[99.999]]]),
                 DALI_zarray = np.array([[[99.999]]]),
                 DALI_zdiffarray=np.array([[[99.999]]])):
        super(H2Model, self).__init__(Mstar, Av, d, param_dict)
        self.grid_dict = self.grid_from_radius()
        if DALI_rgrid[0] != 99.999:
            self.grid_dict["rgrid"] = DALI_rgrid

        if len(np.shape(DALI_Tgrid)) > 1:
            self.T = DALI_Tgrid
        else:
            self.T = self.radial_T(self.grid_dict["rgrid"])

        # self.DALI_ncol = DALI_ncol
        # self.z_diff_array = DALI_zdiffarray
        # self.z_array = DALI_zarray

    def unred_flux_to_H2(self):

        #Define color excess E(B-V)
        EBV = self.Av / Rv

        LyA_fluxes = np.array([self.param_dict['LyA01flux'], self.param_dict['LyA02flux'],
                               self.param_dict['LyA14flux'], self.param_dict['LyA17flux']])
        self.LyA_fluxes = LyA_fluxes
        flux_unreddened = ccm_unred(self.LyA_pumping, LyA_fluxes, EBV)

        flux_LyA_target = flux_unreddened * (self.d * PCtoAU / self.r_LyA) ** 2.

        rgrid = self.grid_dict["rgrid"]

        flux_at_H2 = []
        for flux in flux_LyA_target:
            prog_flux_at_H2 = [flux * (self.r_LyA / r) ** 2. for r in rgrid]
            flux_at_H2.append(prog_flux_at_H2)

        return np.array(flux_at_H2).astype(np.float64)

    def trans_temp(self):

        exp_trans = np.zeros((len(self.T), len(self.lower_energy)))
        for E_idx, E in enumerate(self.lower_energy):
            exp_trans[:, E_idx] = np.exp(-1. * E / (kB * self.T))

        return exp_trans

    def sigma_0(self):

        sig0 = (np.power(self.LyA_pumping * ANGtoCM, 3.) / (8. * math.pi * c)) * (self.g_progression / self.g_initial) * self.Aul_total

        return sig0

    def match_grid_to_file(self, prop, grid_val):
        """Helper method

        """
        return np.argmin(np.absolute(grid_val - prop))

    def match_2D_grid_to_file(self, prop1, prop2, grid_val1, grid_val2):
        """Helper method
        """
        tot_diff = np.absolute(grid_val1 - prop1) + np.absolute(grid_val2 - prop2)
        return np.argmin(tot_diff)

    def partition_function(self):
        """Do this by progression (instead of grid point)
        """
        T_H2_range = range(0, np.shape(self.Z_TH2)[1]) #Partition function defined for T in increments of 1 Kelvin
        part_fnct = np.zeros((len(self.T), len(self.J_initial)))
        for prog_idx in range(0, np.shape(self.Z_TH2)[0]):
            T_array = [self.Z_TH2[prog_idx, self.match_grid_to_file(T_H2_range, T)] if T != 0. else 0. for T in self.T]
            part_fnct[:, prog_idx] = T_array

        return part_fnct

    def H2_number_density(self):

        num_density_constant = ((self.param_dict["MH2"] * MSUNtoG * (2. - self.param_dict["gamma"]))
                                / (self.mmw_H2 * np.sqrt(8. * math.pi ** 3.)
                                   * (self.param_dict["rchar"] * AUtoCM) ** 2.))

        rho_surface = self.surface_density()
        rho_disk = self.exp_disk_density()
        thermal_distr = self.trans_temp()
        part_fnct = self.partition_function()
        Hz = self.scale_height()

        n_H2 = np.zeros((len(self.J_initial), len(rho_disk), len(rho_surface)))

        for prog_idx in range(0, len(self.J_initial)):

            part_0_loc = np.where(part_fnct[:, prog_idx] == 0.)[0]
            correct_factor = np.array([0. if r_idx in part_0_loc else (num_density_constant
                                                              * rho_surface[r_idx] * thermal_distr[r_idx, prog_idx]
                                                              / (part_fnct[r_idx, prog_idx] * Hz[r_idx] * AUtoCM) * self.g_initial[prog_idx]) for r_idx in range(0, len(rho_surface))]).astype(np.float64)
            for z_idx, z in enumerate(rho_disk):
                n_H2[prog_idx, z_idx, :] = z * correct_factor

        return n_H2

    def opacity(self, n_H2):

        sigma = self.sigma_0()
        kappa_H2 = []
        for prog_idx in range(0, len(self.J_initial)):
            kappa_H2.append(sigma[prog_idx] * n_H2[prog_idx, :, :])

        return np.array(kappa_H2).astype(np.float64)

    def effective_tau(self, NH2, prog_idx, z_idx):

        H2_col_density = np.array([np.log10(col) if col != 0. else 0. for col in NH2]).astype(np.float64)

        step_T = 100
        step_N = 0.2
        H2_col_density = np.rint(H2_col_density / step_N) * step_N

        tau_eff = np.zeros(len(H2_col_density))

        if len(np.shape(self.T)) == 1:
            tau_indices = [self.match_2D_grid_to_file(self.column_tau, self.T_tau, H2_col_density[r_idx], Tr) for r_idx, Tr in enumerate(self.T)]
            tau_eff = [self.eff_tau[prog_idx, best_idx] if self.eff_tau[prog_idx, best_idx] != 0. else 1. for best_idx in tau_indices]
        else:
            tau_eff = []
            tau_indices = []
            for r_idx, H2_col in enumerate(H2_col_density):
                if self.T[r_idx, z_idx] >= 1000. and self.T[r_idx, z_idx] <= 5000. and H2_col >= 10. and H2_col <= 24.:
                    best_idx = self.match_2D_grid_to_file(self.column_tau, self.T_tau, H2_col, self.T[r_idx, z_idx])
                    tau_eff.append(self.eff_tau[prog_idx, best_idx], )
                    tau_indices.append(best_idx)
                else:
                    tau_eff.append(0.)
                    tau_indices.append(0)

        return np.array(tau_eff).astype(np.float64), tau_indices

    def optical_depth(self):

        n_H2 = self.H2_number_density()

        Hp = self.scale_height()
        kappa_H2 = self.opacity(n_H2)
        tau = np.zeros(np.shape(n_H2))
        NH2 = np.zeros(np.shape(n_H2))
        t_eff = np.zeros(np.shape(n_H2))

        for prog_idx in range(0, len(kappa_H2[:, 0, 0])):
            for z_idx in range(0, len(kappa_H2[0, :, 0])):
                NH2[prog_idx, z_idx, :] = (float(z_idx + 1) / len(kappa_H2[0, :, 0])) * n_H2[prog_idx, z_idx, :] * Hp * AUtoCM
                tau[prog_idx, z_idx, :] = kappa_H2[prog_idx, z_idx, :] * Hp * AUtoCM * (float(z_idx + 1) / len(kappa_H2[0, :, 0]))
                t_eff[prog_idx, z_idx, :], tau_indices = self.effective_tau(NH2[prog_idx, z_idx, :], prog_idx, z_idx) #z_idx is ignored in this case

        tau = np.cumsum(tau[::-1], axis=2)[::-1] #Cumulative sum in z, increasing toward midplane
        tau = tau * t_eff

        return tau

    def total_intensity(self):

        EBV = -1. * self.Av / Rv
        LyA_flux_unred = self.unred_flux_to_H2()

        ang_grid = self.angular_grid(self.grid_dict["rgrid"])
        tau = self.optical_depth()
        sight_line = self.sightline_r_phi_z(self.grid_dict["rgrid"], self.grid_dict["phigrid"])

        intens_H2 = []
        for prog_idx in range(0, len(self.J_initial)):
            flux_H2_disk = np.zeros(len(tau[prog_idx, 0, :]))
            for z_idx in range(0, len(tau[prog_idx, :, 0])):
                for r_idx in range(0, len(tau[prog_idx, 0, :])):
                    LyA = LyA_flux_unred[prog_idx, r_idx]
                    flux_emission = LyA * self.eta * (1. - np.exp(-1. * tau[prog_idx, z_idx, r_idx]))
                    flux_H2_disk[r_idx] += flux_emission
            flux_H2_obs = [flux_H2_disk * np.square(ang_grid / (sight_line[:, phi_idx] * PCtoAU)) * np.square(math.cos(self.param_dict["Inclination"])) for phi_idx, phi in enumerate(self.grid_dict["phigrid"])]
            intens_H2.append(flux_H2_obs)

        intens_H2 = 2. * np.array(intens_H2).astype(np.float64)
        intens_H2 = np.swapaxes(intens_H2, 1, 2)

        return intens_H2

    def all_emission_line_intensities(self):

        start_time = time.time()
        EBV = -1. * self.Av / Rv

        intensity = np.zeros((len(np.arange(-300., 300., 1.)), len(self.lambda_props["Wavelength"])))
        lineprof = self.total_intensity()

        intensity_binned = self.collapsed_line_profile(lineprof)

        for wave_idx, wave in enumerate(self.lambda_props["Wavelength"]):
            prog_idx = np.argmin(np.absolute(self.lambda_props["Jl"][wave_idx] - self.J_initial))
            intens_unred = intensity_binned[:, prog_idx] * self.lambda_props["Bul"][wave_idx]

            intens_red = ccm_unred(np.zeros(len(intens_unred)) + wave, intens_unred, EBV)

            #Smooth/convolve with LSF
            lsfx, lsfy = cos_lsf_arulanantham(wave, "LTP1", False)

            lsfy_norm = lsfy.flatten() / np.sum(lsfy.flatten()) #Normalize LSF

            intensity[:, wave_idx] = np.convolve(intens_red, lsfy_norm, mode="same")

        end_time = time.time()
        print ("Python finishes in %s seconds!" % (end_time - start_time))
        return intensity

    def __del__(self):
        # print "Deleted a child!"
        super(H2Model, self).__del__()

class COModel(CircumstellarDisk):
    #States to model
    pump_props = pd.read_csv("H2_PPD_modeling_code/UVCO_pump_props_143_PQR.csv")
    #energy of the lower state of the transition, in ergs
    pump_props["lower_energy"] = np.array([val*1.60217657E-12 for val in pump_props["lower_energy_eV"]]).astype(np.float64)
    #ground vibration state, before LyA pumping
    v1 = 0
    #Disk filling fraction
    eta = 0.25
    #Mean molecular weight of CO
    mmw_CO = 28. * 1.66054e-24 #[g]
    lambda_props = pd.read_csv("H2_PPD_modeling_code/UVCO_lambda_props_143.csv")
    lambda_props["Bul"] = lambda_props["Aul"] / lambda_props["Aul_total"]
    lambda_props["gl"] = np.array([2. * Jl + 1 if Jl % 2 == 0 else 3. * (2. * Jl + 1.) for Jl in lambda_props["Jl"]]).astype(np.float64)
    lambda_props["glp1"] = np.array([2. * Ju + 1 if Ju % 2 == 0 else 3. * (2. * Ju + 1.) for Ju in lambda_props["Ju"]]).astype(np.float64)

    Z_df = pd.read_csv("H2_PPD_modeling_code/disk_partition_values_CO_143_PQR.csv") #disk_partition_values.csv - need to update these for CO (using H2 as proxy for now)
    Z_TCO = []
    for key in list(Z_df):
        Z_TCO.append(np.array(Z_df[key]).astype(np.float64))
    Z_TCO = np.array(Z_TCO).astype(np.float64)

    column_tau_df = pd.read_csv("H2_PPD_modeling_code/column_tau_CO_143_PQR.csv") #column_tau.csv
    column_tau = []
    for key in list(column_tau_df):
        column_tau.append(np.array(column_tau_df[key]).astype(np.float64))
    column_tau = np.array(column_tau).astype(np.float64)

    T_tau_df = pd.read_csv("H2_PPD_modeling_code/T_tau_CO_143_PQR.csv") #T_tau.csv
    T_tau = []
    for key in list(T_tau_df):
        T_tau.append(np.array(T_tau_df[key]).astype(np.float64))
    T_tau = np.array(T_tau).astype(np.float64)

    eff_tau_df = pd.read_csv("H2_PPD_modeling_code/eff_tau_CO_143_PQR.csv") #eff_tau.csv - need to update for CO (using H2 as proxy for now)
    eff_tau = []
    for key in list(eff_tau_df):
        eff_tau.append(np.array(eff_tau_df[key]).astype(np.float64))
    eff_tau = np.array(eff_tau).astype(np.float64)

    def __init__(self, Mstar, Av, d, param_dict):
        super(COModel, self).__init__(Mstar, Av, d, param_dict)
        self.grid_dict = self.grid_from_radius()
        self.T = self.radial_T(self.grid_dict["rgrid"])

    def unred_flux_to_CO(self):

        #Define color excess E(B-V)
        EBV = self.Av / Rv

        self.pump_props["LyA_fluxes_scaled"] = np.zeros(len(self.pump_props["LyA_pumping"])).astype(np.float64)

        for pump_wave_idx, pump_wave in enumerate(self.pump_props["LyA_pumping"]):
            if pump_wave in self.param_dict["LyA_pumping"]:
                min_diff_loc = np.argmin(np.absolute(pump_wave - self.param_dict["LyA_pumping"]))
                self.pump_props.loc[pump_wave_idx, "LyA_fluxes_scaled"] = self.param_dict["LyA_fluxes"][min_diff_loc]

        prog_flux_unreddened = ccm_unred(self.pump_props["LyA_pumping"], self.pump_props["LyA_fluxes_scaled"], EBV)

        flux_LyA_target = prog_flux_unreddened * (self.d * PCtoAU / self.r_LyA) ** 2.

        flux_at_CO = []
        for flux in flux_LyA_target:
            prog_flux_at_CO = [flux * (self.r_LyA / r) ** 2. for r in self.grid_dict["rgrid"]]
            flux_at_CO.append(prog_flux_at_CO)

        flux_at_CO = np.array(flux_at_CO).astype(np.float64)

        return flux_at_CO

    def trans_temp(self):

        exp_trans = np.zeros((len(self.T), len(self.pump_props["lower_energy"])))
        for E_idx, E in enumerate(self.pump_props["lower_energy"]):
            exp_trans[:, E_idx] = [np.exp(-1. * E / (kB * T)) if T != 0. else 0. for T in self.T]

        return exp_trans

    def sigma_0(self):

        sig0 = (np.power(self.pump_props["LyA_pumping"] * ANGtoCM, 3.)
                / (8. * math.pi * c)) * (self.pump_props["g_progression"] / self.pump_props["g_initial"]) * self.pump_props["Aul_total"]

        return sig0

    def match_grid_to_file(self, prop, grid_val):
        """Helper method
        """
        return np.argmin(np.absolute(grid_val - prop))

    def match_2D_grid_to_file(self, prop1, prop2, grid_val1, grid_val2):
        """Helper method
        """
        prop1_diff = np.absolute(grid_val1 - prop1)
        prop2_diff = np.absolute(grid_val2 - prop2)
        tot_diff = prop1_diff + prop2_diff
        min_loc = np.argmin(tot_diff)

        return np.argmin(tot_diff)

    def partition_function(self):
        """Do this by progression (instead of grid point)
        """
        part_fnct = np.zeros((len(self.T), np.shape(self.Z_TCO)[0]))
        for prog_idx in range(np.shape(part_fnct)[1]):
            for T_idx, T in enumerate(self.T):
                if T != 0.:
                    match_loc = self.match_grid_to_file(self.T_tau[prog_idx, :], T)
                    part_fnct[T_idx, prog_idx] = self.Z_TCO[prog_idx, match_loc]
                else:
                    part_fnct[T_idx, prog_idx] = np.nan
        return part_fnct

    def CO_mass_density(self):

        num_density_constant = ((self.param_dict["MCO"] * MSUNtoG * (2. - self.param_dict["gamma"]))
                                / (self.mmw_CO * np.sqrt(8. * math.pi ** 3.)
                                   * (self.param_dict["rchar"] * AUtoCM) ** 2))
        XCO = 1.e-4

        rho_surface = num_density_constant * XCO * self.surface_density() #len = 201
        rho_disk = self.exp_disk_density() #shape = (11)
        Hz = self.scale_height()

        mass_density = np.zeros((len(rho_surface), len(rho_disk)))
        for z_idx in range(len(rho_disk)):
            mass_density[:, z_idx] = (rho_surface / (Hz * AUtoCM)) * rho_disk[z_idx] #rho_disk[:, z_idx] #np.exp(-0.5 * np.square(self.param_dict["z"] * z_val / Hz))

        return mass_density

    def CO_number_density(self):

        mass_density = self.CO_mass_density()
        part_fnct = self.partition_function()

        n_CO = np.zeros((len(self.pump_props["J_progression"]), np.shape(mass_density)[1], np.shape(mass_density)[0]))
        for prog_idx in range(len(self.pump_props["J_progression"])):
            for z_idx in range(np.shape(mass_density)[1]):
                n_CO[prog_idx, z_idx, :] = (mass_density[:, z_idx] / self.mmw_CO) * part_fnct[:, prog_idx] * self.pump_props.iloc[prog_idx]["g_initial"]

        return n_CO

    def opacity(self, n_CO):

        sigma = self.sigma_0()
        pump_props_sorted = self.pump_props.groupby(["J_progression"])
        kappa_CO_full = np.array([sigma[prog_idx] * n_CO[prog_idx, :, :] for prog_idx in range(len(self.pump_props["J_progression"]))]).astype(np.float64)

        return kappa_CO_full

    def effective_tau(self, NCO, prog_idx):

        CO_col_density = np.log10(NCO)

        step_N = 0.2
        CO_col_density = np.rint(CO_col_density / step_N) * step_N

        tau_eff = np.zeros(len(CO_col_density))

        tau_indices = [self.match_2D_grid_to_file(self.column_tau[prog_idx, :], self.T_tau[prog_idx, :], CO_col_density[r_idx], Tr) for r_idx, Tr in enumerate(self.T)]

        tau_eff = [self.eff_tau[prog_idx, best_idx] if best_idx != 0. else 1. for best_idx in tau_indices]

        return tau_eff

    def optical_depth(self):

        n_CO = self.CO_number_density()
        kappa_CO = self.opacity(n_CO)
        Hp = self.scale_height()

        tau = np.zeros(np.shape(n_CO))
        NCO = np.zeros(np.shape(n_CO))
        t_eff = np.zeros(np.shape(n_CO))

        for prog_idx in range(0, len(kappa_CO[:, 0, 0])):
            for z_idx in range(0, len(kappa_CO[0, :, 0])):
                NCO[prog_idx, z_idx, :] = n_CO[prog_idx, z_idx, :] * Hp * AUtoCM * (float(z_idx + 1) / 11.)
                tau[prog_idx, z_idx, :] = kappa_CO[prog_idx, z_idx, :] * Hp * AUtoCM * (float(z_idx + 1) / 11.)
                t_eff[prog_idx, z_idx, :] = self.effective_tau(NCO[prog_idx, z_idx, :], prog_idx)

        tau = np.cumsum(tau, axis=1)
        tau = tau * t_eff

        return tau

    def total_intensity(self):

        EBV = -1. * self.Av / Rv
        LyA_flux_unred = self.unred_flux_to_CO()

        ang_grid = self.angular_grid(self.grid_dict["rgrid"])
        tau = self.optical_depth()

        sight_line = self.sightline_r_phi_z(self.grid_dict["rgrid"], self.grid_dict["phigrid"])
        vobs = self.vrad_to_vobs(self.grid_dict["rgrid"], self.grid_dict["phigrid"])

        flux_CO_disk = np.array([LyA_flux_unred[prog_idx, :] * self.eta * (1. - np.exp(-1. * tau[prog_idx, :, :])) for prog_idx in range(len(self.pump_props["J_progression"]))]).astype(np.float64)
        flux_CO_disk = np.sum(flux_CO_disk, axis=1) #Sum fluxes in z

        intens_CO = np.zeros((np.shape(flux_CO_disk)[0], np.shape(sight_line)[0], np.shape(sight_line)[1]))

        for phi_idx, phi in enumerate(self.grid_dict["phigrid"]):
            flux_CO_obs = np.array([2. * flux_CO_disk[prog_idx, :] * np.square(ang_grid * math.cos(self.param_dict["Inclination"])) / np.square(sight_line[:, phi_idx] * PCtoAU)
                                    for prog_idx in range(len(self.pump_props["J_progression"]))]).astype(np.float64)
            intens_CO[:, :, phi_idx] = flux_CO_obs

        return intens_CO

    def all_emission_line_intensities(self):

        EBV = -1. * self.Av / Rv

        intensity = np.zeros((len(np.arange(-300., 300., 1.)), len(self.lambda_props["Wavelength"])))
        lineprof = self.total_intensity()

        intensity_binned = self.collapsed_line_profile(lineprof)

        for wave_idx, wave in enumerate(self.lambda_props["Wavelength"]):
            prog_idx = np.argmin(np.absolute(self.lambda_props["Jl"][wave_idx] - self.pump_props["J_initial"]))
            intens_unred = intensity_binned[:, prog_idx] * self.lambda_props["Bul"][wave_idx]

            intens_red = ccm_unred(np.zeros(len(intens_unred)) + wave, intens_unred, EBV)

            #Smooth/convolve with LSF
            lsfx, lsfy = cos_lsf_arulanantham(wave, "LTP1", False)

            lsfy_norm = lsfy.flatten() / np.sum(lsfy.flatten()) #Normalize LSF

            intensity[:, wave_idx] = np.convolve(intens_red, lsfy_norm, mode="same")

        return intensity
