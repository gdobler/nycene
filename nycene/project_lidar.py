#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
from utils import *
from optimize_camera import *


def project_lidar(params, raster, imsize, rad=4000, buff=200, behind=500):
    """
    ADD DOCS!!!

    params : camera parameters
    raster : name of the raster file
    imsize : the size of the image
    rad    : max distance in front of camera plane to project
    buff   : min distance in front of the camera plane to project
    behind : distance north of the camera plane to include
    """

    # -- check raster
    if type(raster) is str:
        rast_full, hdr = read_raster(raster)
        xvec_full      = hdr.cmin + np.arange(hdr.ncol) # NYS
        yvec_full      = hdr.rmin + np.arange(hdr.nrow) # NYS
    else:
        print("PROJECT_LIDAR: input raster must be string for now.")
        return None

    # -- subselect the rast that is in front of the cmaera
    # GGD: FOR NOW I KNOW WHICH WAY THE CAMERA IS FACING!!!
    camr = int(params[4] - hdr.rmin)
    camc = int(params[3] - hdr.cmin)
    rast = rast_full[:camr + behind] # go 500 ft north camera may be tilted
    xvec = xvec_full
    yvec = yvec_full[:camr + behind]

    # -- sub-select raster points
    # GGD: OBJECTS BEHIND THE CAMERA ARE PROJECTING INTO THE IMAGE!!!
    t0 = time.time()
    rast_sub   = rast[-(rad + buff + behind) : -(buff + behind)]
    xvec_sub   = xvec
    yvec_sub   = yvec[-(rad + buff + behind) : -(buff + behind)]
    nrow, ncol = rast_sub.shape
    xvec_rast  = xvec_sub[(np.arange(nrow * ncol) % ncol) \
                              .reshape(nrow, ncol).flatten()]
    yvec_rast  = yvec_sub[(np.arange(nrow * ncol) // ncol) \
                              .reshape(nrow, ncol).flatten()]
    xyz        = np.vstack([xvec_rast, yvec_rast, rast_sub.flatten()]).T
    print("time to initialize points  : {0}s".format(time.time() - t0))

    # -- project points and recenter
    t0 = time.time()
    cuv  = colin(params, xyz)
    print("time to project to uv      : {0}s".format(time.time() - t0))
    t0   = time.time()
    uv   = ((cuv - np.array([0.5 * imsize[0], -0.5 * imsize[1]])) /
            np.array([-1, 1])).round().astype(int)
    print("time to recenter           : {0}s".format(time.time() - t0))

    # -- find points within the image
    t0 = time.time()
    gind = (uv[:, 0] >= 0) & (uv[:, 0] < imsize[0]) & (uv[:, 1] >= 0) & \
        (uv[:, 1] < imsize[1])

    # -- pull off parameters from these points
    us, vs     = uv[gind].T
    xs, ys, zs = xyz[gind].T
    dist2      = (xs - params[3])**2 + (ys - params[4])**2

    # -- fill image
    d2img = np.ones(imsize) * 1e20
    ximg  = np.zeros(imsize)
    yimg  = np.zeros(imsize)

    for ii, (uu, vv) in enumerate(zip(us, vs)):
        if d2img[uu, vv] > dist2[ii]:
            d2img[uu, vv] = dist2[ii]
            ximg[uu, vv]  = xs[ii]
            yimg[uu, vv]  = ys[ii]
    print("time to map to image       : {0}s".format(time.time() - t0))

    # -- remove overhangs
    t0 = time.time()
    for ii in range(1, d2img.shape[0]):
        for jj in range(d2img.shape[1]):
            if d2img[ii - 1, jj] < d2img[ii, jj]:
                d2img[ii, jj] = d2img[ii - 1, jj]
                ximg[ii, jj]  = ximg[ii - 1, jj]
                yimg[ii, jj]  = yimg[ii - 1, jj]

    # -- remove narrow gaps
    for ii in range(d2img.shape[0]):
        for jj in range(1, d2img.shape[1] - 1):
            if (d2img[ii, jj - 1] < d2img[ii, jj]) & \
                    (d2img[ii, jj + 1] < d2img[ii, jj]):
                d2img[ii, jj] = d2img[ii, jj - 1]
                ximg[ii, jj]  = ximg[ii, jj - 1]
                yimg[ii, jj]  = yimg[ii, jj - 1]

    # -- remove single points in the sky
    for ii in range(d2img.shape[0]):
        for jj in range(1, d2img.shape[1] - 1):
            if (d2img[ii, jj - 1] == d2img[ii, jj + 1]):
                d2img[ii, jj] = d2img[ii, jj - 1]
                ximg[ii, jj] = ximg[ii, jj - 1]
                yimg[ii, jj] = yimg[ii, jj - 1]
    print("time to post-process image : {0}s".format(time.time() - t0))

    return ximg, yimg, np.sqrt(d2img)
