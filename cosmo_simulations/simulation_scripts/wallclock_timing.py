from unyt import ms, hour
import matplotlib.pyplot as plt
import numpy as np
import sys

def get_timesteps(timesteps_file):
    f = np.genfromtxt(timesteps_file, usecols=(2,12), invalid_raise = False)
    
    a = f[:,0] # Scale-factor
    t = f[:,1] * ms # Wallclock time
    t = np.cumsum(t).to(hour)
    return a, t

filename = sys.argv[1]
a, t = get_timesteps(filename)
np.savetxt('wallclock.txt', np.asarray([a, t]), delimiter = ',')
