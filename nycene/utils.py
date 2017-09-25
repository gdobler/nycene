#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pyproj
import numpy as np
import pandas as pd


def colin(params, xyz):

    # -- Unwrap params
    kappa, phi, omega, xs, ys, zs, ff = params

    omega = float(omega)
    phi   = float(phi) + 0.5*np.pi
    kappa = float(kappa)
    xs    = float(xs)
    ys    = float(ys)
    zs    = float(zs)
    ff    = float(ff)

    # -- utilities
    co = np.cos(omega)
    so = np.sin(omega)
    cp = np.cos(phi)
    sp = np.sin(phi)
    ck = np.cos(kappa)
    sk = np.sin(kappa)

    a1 =  cp * ck + sp * so * sk
    b1 =  cp * sk + sp * so * ck
    c1 =  sp * co
    a2 = -co * sk
    b2 =  co * ck
    c2 =  so
    a3 =  sp * ck + cp * so * sk
    b3 =  sp * sk - cp * so * ck
    c3 =  cp * co

    # -- project to image plane
    unum  = a1 * (xyz[:,0] - xs) + b1 * (xyz[:,1] - ys) + c1 * (xyz[:,2] - zs)
    vnum  = a2 * (xyz[:,0] - xs) + b2 * (xyz[:,1] - ys) + c2 * (xyz[:,2] - zs)
    denom = a3 * (xyz[:,0] - xs) + b3 * (xyz[:,1] - ys) + c3 * (xyz[:,2] - zs)

    uu =  ff * unum / denom
    vv = -ff * vnum / denom

    return np.vstack([uu, vv]).T



def get_fiducials(fname, df=False):
    """ Get the fiducials """

    # -- set the path
    fname = os.path.join("..", "data", fname)

    # -- read the fiducials file
    fids = pd.read_csv(fname)

    # -- return full dataframe if desired, else parse the file
    if df:
        return fids
    else:
        return fids[["x", "y", "z"]].values, fids[["u", "v"]].values



def latlon_to_ny(lat, lon):
    """ lat/lon to NYS """

    proj = pyproj.Proj(init="epsg:2263", preserve_units=True)
    result = proj(lon, lat)

    return result



def ny_to_latlon(ny_lat, ny_lon):
    """ NYS to lat/lon """

    proj   = pyproj.Proj(init="epsg:2263", preserve_units=True)
    result = proj(ny_lat, ny_lon, inverse=True)

    return result[1], result[0]
