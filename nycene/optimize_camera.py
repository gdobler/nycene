#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import minimize
from utils import *


def loss(params, xyz, uv):
    """
    Define the loss.
    """

    # -- unpack params
    omega, phi, kappa, xs, ys, zs, ff = params

    # -- set limits
    if (omega < 0.0) or (omega >= 2.0 * np.pi):
        return 1e9 + omega**4
    elif (phi < -0.5 * np.pi) or (phi >= 0.5 * np.pi):
        return 1e9 + phi**4
    elif (kappa < 0.0) or (kappa >= 2.0 * np.pi):
        return 1e9 + kappa**4
    elif zs < 0.0:
        return 1e9 + zs**4
    elif ff < 0.0:
        return 1e9 + ff**4
    elif xs > 1.0e10:
        return 1e9 + xs**4
    elif ys > 1.0e10:
        return 1e9 + ys**4

    return ((colin(params, xyz.astype(float)) - uv)**2).sum()



def optimize_camera(guess, xyz, uv, imsize, niter=100, method="Powell", 
                    ftol=1e-6, verbose=False):
    """
    Optimize the camera parameters.
    """

    # -- utilities
    params = guess.copy()
    score  = 1e10

    # -- center image coordinates
    cuv = uv * np.array([-1, 1]) + \
        np.array([0.5 * imsize[0], -0.5 * imsize[1]])

    # -- minimize loss
    for ii in range(niter):
        res = minimize(loss, params, args=(xyz, cuv), method=method, 
                       options={"ftol":ftol})
        if res.fun < score:
            score  = res.fun
            params = res.x
            if verbose:
                print("params, score : {0}, {1}".format(params, score))

    return params, score


if __name__ == "__main__":

    # -- make a guess for the camera parameters
    lat    = 40.689872
    lon    = -73.988305
    xx, yy = latlon_to_ny(lat, lon)
    zz     = 400.
    kappa  = 0.5 * np.pi
    phi    = 0.0
    omega  = 0.0
    ff     = 800.
    guess  = np.array([kappa, phi, omega, xx, yy, zz, ff])

    # -- get the fiducials and set the image size
    xyz, uv = get_fiducials("day_ref.csv")
    imsize  = 2120, 4056

    # -- solve for the best fit solution
    params, score = optimize_camera(guess, xyz, uv, imsize)

    print params, score
