import sys
import sys
sys.path.append('../Bayesflow/')

import numpy as np
import numba
from numba import njit
import keras
import matplotlib.pyplot as plt
import seaborn as sns
import notebook

RNG = np.random.default_rng(2024)

def meta():
	n_obs = np.random.randint(15, 21)
	return dict(n_obs = n_obs)

def prior():
	v  = RNG.gamma(2, 1/5)    
	a  = RNG.gamma(6, 1/12)
	t0 = RNG.gamma(4, 1/15) 
	return dict(v=v, a=a, t0=t0)
 

@njit
def ddm(v, a, t0, n_obs, z=0.0, dt=0.001, s=1, max_iter=20000):
    c = np.sqrt(dt * s)
    out = np.zeros((n_obs, 2))
    
    for n in range(n_obs):
        n_iter = 0
        x = a * z

        while x <= a and n_iter < max_iter:
            x += v*dt + c*np.random.randn()
            n_iter += 1
            rt = n_iter * dt
            if n_iter == max_iter and x < a:
                out[n, 0] = 0
                out[n, 1] = 1 # is it timeout?/ Yes
            else:
                out[n, 0] = rt + t0
                out[n, 1] = 0 # is it timeout?/ No
    return out


def likelihood(v, a, t0, n_obs):
    out = ddm(v, a, t0, n_obs)
    out1 = dict(rts=out[:, 0], timeouts=out[:,1])
    return out1