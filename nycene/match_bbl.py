#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import geopandas as gp
from shapely.geometry import Point


def match_bbl(ximg, yimg, mappluto):
    """
    ADD DOCS!!!

    xgrid : an "image" of x coordinates in NYS
    ygrid : an "image" of y coordinates in NYS
    mappluto : the location of the mappluto shape file to use or geodataframe
    """

    # -- read the mappluto data
    if type(mappluto) == str:
        print("reading MapPLUTO data...")
        mp = gp.read_file(mappluto)
    else:
        mp = mappluto

    # -- extract centroids, geometries, and BBLs
    mp_cenx = np.array([i.x for i in mp['geometry'].centroid])
    mp_ceny = np.array([i.y for i in mp['geometry'].centroid])
    mp_geo  = np.array(mp['geometry'])
    mp_bbl  = np.array(mp['BBL'])


    # print("making pnts geodataframe...")
    # pnts   = [Point(i, j) for i, j in zip(ximg.flatten(), yimg.flatten())]
    # dfpnts = gp.GeoDataFrame({"pix" : np.arange(ximg.size), "geometry" : pnts})

    # print("performing spatial join...")
    # return gp.sjoin(mp, dfpnts, how="inner", op="contains")

    # -- initialize the BBL image
    nr, nc = ximg.shape
    bimg   = np.zeros((nr, nc))

    # -- loop through pixels
    rad2 = 500.**2
    for ii in range(nr):
        print "row {0} of {1}\r".format(ii + 1, nr),
        sys.stdout.flush()
        for jj in range(nc):
            xpos = ximg[ii, jj]
            ypos = yimg[ii, jj]
            pnt  = Point(xpos, ypos)
            ind  = ((mp_cenx - xpos)**2 + (mp_ceny - ypos)**2) < rad2
            for geo, bbl in zip(mp_geo[ind], mp_bbl[ind]):
                if geo.contains(pnt):
                    bimg[ii, jj] = bbl
                    continue

    return bimg
