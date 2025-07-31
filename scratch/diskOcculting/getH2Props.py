# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 13:40:31 2025

@author: ecruzaguirre
"""

import numpy as np
from scipy.io import readsav

def getLineProps(wave,tolr):
    h2data=readsav("C:/Users/ecruzaguirre/Documents/HUSL/H2Spec/h2SyntheticSpectra-master/h2SyntheticSpectra-master/fluormod_trans_python.idlsav")
    h2wave=h2data.h2_trans.wave #wavelength array
    h2diff=abs(h2wave-wave) #get separation from the target wavelength
    h2indx=np.where(h2diff<=tolr)[0] #find all emission lines within the wavelength tolerance
    for i in h2indx:
        idWave=h2wave[i]
        idAulL=h2data.h2_trans.avalue[i]
        idJlLn=h2data.h2_trans.jl[i]
        idJuLn=h2data.h2_trans.ju[i]
        idAulT=h2data.h2_trans.atotal[i]
        print(f'Wavelength: {idWave} Ang')
        print(f'Aul Line  : {idAulL:.1e} Hz')
        print(f'Jl Line   : {idJlLn}')
        print(f'Ju Line   : {idJuLn}')
        print(f'Aul Total : {idAulT:.1e} Hz\n')

getLineProps(1431.010,0.01)
getLineProps(1556.870,0.01)