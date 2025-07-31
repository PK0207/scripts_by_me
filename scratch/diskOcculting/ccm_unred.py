"""Function to deredden spectrum (translated to Python from IDL ccm_unred.pro)
Written by N. Arulanantham
Last Edited: July 2018
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io.idl import readsav

def which_spec_range(wave):

    #Infrared
    if wave > 0 and wave < 1.1:
        a = 0.574 * np.power(wave, 1.61)
        b = -0.527 * np.power(wave, 1.61)
        return a, b
    #Optical/NIR, with new constants from O'Donnell (1994)
    elif wave >= 1.1 and wave < 3.3:
        #Numpy takes these polynomial coefficients in reverse order!
        c1 = [1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
        c2 = [0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]

        y = wave - 1.82
        a = np.polyval(c1[::-1], y)
        b = np.polyval(c2[::-1], y)
        return a, b
    #Mid-UV, range 1
    elif wave >= 3.3 and wave <= 5.9:
        a = 1.752 - 0.316 * wave - (0.104 / ((wave - 4.67) ** 2. + 0.341))
        b = -3.090 + 1.825 * wave + (1.206 / ((wave - 4.62) ** 2. + 0.263))
        return a, b
    #Mid-UV, range 2
    elif wave > 5.9 and wave < 8:
        y = wave - 5.9
        Fa = -0.04473 * y ** 2. - 0.009779 * y ** 3.
        Fb = 0.2130 * y ** 2. + 0.1207 * y ** 3.

        a = 1.752 - 0.316 * wave - (0.104 / ((wave - 4.67) ** 2. + 0.341)) + Fa
        b = -3.090 + 1.825 * wave + (1.206 / ((wave - 4.62) ** 2. + 0.263)) + Fb
        return a, b
    #Far-UV
    elif wave >= 8. and wave <= 11.:
        y = wave - 8.
        #IDL takes polynomial coefficients in the opposite direction as Python!
        c1 = [-1.073, -0.628, 0.137, -0.070]
        c2 = [13.670, 4.257, -0.420, 0.374]

        a = np.polyval(c1[::-1], y)
        b = np.polyval(c2[::-1], y)
        return a, b
    return np.nan, np.nan

def ccm_unred(wave_array, flux_array, ebv, Rv=3.1):
  """Deredden fluxes from a spectrum
  :param wave_array: wavelength array [Angstroms]
  :param flux_array: flux array, same length as wave
  :param ebv: scalar, color excess
  :param Rv: scalar (optional), ratio of total to selective extinction (Av/ebv)
  :return: unreddened flux array
  """

  wave_nums = 1.e4 / np.array(wave_array).astype(np.float64) #[inverse microns]
  npts = len(wave_nums)

  ab_vals = [list(which_spec_range(wave)) for wave in wave_nums]

  ab_lists = [list(ab) for ab in zip(*ab_vals)]

  a = np.array(ab_lists[0]).astype(np.float64)
  b = np.array(ab_lists[1]).astype(np.float64)

  Av = Rv * ebv
  A_lamda = Av * (a + b / Rv)

  return flux_array * np.power(10., 0.4 * A_lamda)

def main():
    # import pidly
    #
    # idl = pidly.IDL('/Applications/exelis/idl85/bin/idl') #Call to IDL functions
    #
    # idl('.compile /H2_PPD_modeling_code/ccm_unred.pro')
    data_table = readsav("/Users/Nicole/Documents/CU/Research/RY_Lup/Data/Comparison/DETau/DETAU_codd.sav")
    wavelength = np.array(data_table["wave"]).astype(np.float64) #[Angstroms] (probably, for UV data)
    flux_original = np.array(data_table["flux"]).astype(np.float64)
    flux_err = np.array(data_table["err"]).astype(np.float64)

    # path = "/Users/Nicole/Documents/CU/Research/RY_Lup/Scripts/SELFiE_for_Zac/Targets/V4046Sgr/"
    # data_frame = pd.read_csv("".join([path, "V4046SGR_071910kf.csv"]))

    flux_dereddened = ccm_unred(wavelength,
                                flux_original,
                                0.29)

    data_frame = pd.DataFrame({"Wavelength": wavelength,
                               "Flux": flux_dereddened})

    data_frame.to_csv("DETau_forETC.dat",
                      index=False,
                      sep=" ",
                      columns=["Wavelength", "Flux"])

    # IDL_unred = []
    # for flux_idx, flux in enumerate(data_frame["col2"]):
    #     idl('ccm_unred, ' + str(data_frame["col1"].iloc[flux_idx]) + ', '
    #         + str(flux) + ', ' + str(1. / 3.1) + ', IDL_unred')
    #     IDL_unred.append(idl.IDL_unred[0])

    plt.clf()
    plt.close()
    plt.plot(data_frame["Wavelength"],
             data_frame["Flux"],
             color="k", lw=1.5,
             label = "Python Dereddened")
    plt.plot(wavelength,
             flux_original,
             color=(0.8, 0., 0.2), lw=1.,
             label = "Original")
    # plt.plot(data_frame["col1"],
    #          IDL_unred,
    #          color=(0.2, 0., 0.8), lw=0.8,
    #          label = "IDL Dereddened")
    plt.legend(loc="best")
    plt.xlabel("Wavelength [Angstroms]")
    plt.ylabel("Flux [cgs]")
    plt.show(block=True)
    #plt.savefig("dereddening_comp.pdf", format="pdf")

if __name__ == "__main__":
    main()
