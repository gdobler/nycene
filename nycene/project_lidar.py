#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
from utils import *
from optimize_camera import *


t00 = time.time()

# -- read in the raster
t0 = time.time()
rfile          = os.path.join(os.environ["REBOUND_WRITE"], "rasters",
                              "MN_raster.bin")
rast_full, hdr = read_raster(rfile)
xvec_full = hdr.cmin + np.arange(hdr.ncol) # NYS
yvec_full = hdr.rmin + np.arange(hdr.nrow) # NYS
print("time to read raster       : {0}s".format(time.time() - t0))


# -- get the params
t0 = time.time()
lat    = 40.7556029
lon    = -73.98712
xx, yy = latlon_to_ny(lat, lon)
zz     = 400.
kappa  = 1.5 * np.pi
phi    = 0.0
omega  = 0.0
ff     = 800.
guess  = np.array([kappa, phi, omega, xx, yy, zz, ff])
xyz, uv = get_fiducials("audubon_d6.csv")
imsize  = 2560, 1918
params  = optimize_camera(guess, xyz, uv, imsize)[0]
print("time to solve camera      : {0}s".format(time.time() - t0))


# -- subselect the rast that is in front of the cmaera
# GGD: FOR NOW I KNOW WHICH WAY THE CAMERA IS FACING!!!
camr = int(params[4] - hdr.rmin)
camc = int(params[3] - hdr.cmin)
#rast = rast_full[camr - 500:] # go 500 ft north in case camera is tilted
rast = rast_full[:camr + 500] # go 500 ft north in case camera is tilted
xvec = xvec_full
yvec = yvec_full[:camr + 500]


# -- project the nearest 1000 feet
# GGD: OBJECTS BEHIND THE CAMERA ARE PROJECTING INTO THE IMAGE!!!
#rast_sub = rast[:1000]
t0 = time.time()
rast_sub = rast[-6700:-700]
xvec_sub = xvec
yvec_sub = yvec[-6700:-700]

nrow, ncol = rast_sub.shape
xvec_rast  = xvec_sub[(np.arange(nrow * ncol) % ncol) \
                          .reshape(nrow, ncol).flatten()]
yvec_rast  = yvec_sub[(np.arange(nrow * ncol) // ncol) \
                          .reshape(nrow, ncol).flatten()]

xyz  = np.vstack([xvec_rast, yvec_rast, rast_sub.flatten()]).T
print("time to initialize points : {0}s".format(time.time() - t0))

t0 = time.time()
cuv  = colin(params, xyz)
print("time to project to uv     : {0}s".format(time.time() - t0))

t0 = time.time()
uv   = ((cuv - np.array([0.5 * imsize[0], -0.5 * imsize[1]])) /
            np.array([-1, 1])).round().astype(int)
print("time to recenter          : {0}s".format(time.time() - t0))


# -- put into image grid
t0 = time.time()
gind = (uv[:, 0] >= 0) & (uv[:, 0] < imsize[0]) & (uv[:, 1] >= 0) & \
    (uv[:, 1] < imsize[1])

us, vs = uv[gind].T
xs, ys = xvec_rast[gind], yvec_rast[gind]
d2s    = (xs - params[3])**2 + (ys - params[4])**2

dgrid = np.ones(imsize) * 1e20
xgrid = np.zeros(imsize)
ygrid = np.zeros(imsize)

for ii, (uu, vv) in enumerate(zip(us, vs)):
    if dgrid[uu, vv] > d2s[ii]:
        dgrid[uu, vv] = d2s[ii]
        xgrid[uu, vv] = xs[ii]
        ygrid[uu, vv] = ys[ii]
print("time to map to grid       : {0}s".format(time.time() - t0))


# -- process grid
t0 = time.time()
for ii in range(1, dgrid.shape[0]):
    for jj in range(dgrid.shape[1]):
        if dgrid[ii - 1, jj] < dgrid[ii, jj]:
            dgrid[ii, jj] = dgrid[ii - 1, jj]

for ii in range(dgrid.shape[0]):
    for jj in range(1, dgrid.shape[1] - 1):
        if (dgrid[ii, jj - 1] < dgrid[ii, jj]) & \
                (dgrid[ii, jj + 1] < dgrid[ii, jj]):
            dgrid[ii, jj] = dgrid[ii, jj - 1]

for ii in range(dgrid.shape[0]):
    for jj in range(1, dgrid.shape[1] - 1):
        if (dgrid[ii, jj - 1] == dgrid[ii, jj + 1]):
            dgrid[ii, jj] = dgrid[ii, jj - 1]
print("time to post-process grid : {0}s".format(time.time() - t0))


print("----------------------------------------------")
print("total elapsed time        : {0}s".format(time.time() - t00))


foo = dgrid.copy()
foo[foo == 1e20] = 0.0


# import numpy as np
# import sys
# import glob
# from math import sin,cos,sqrt,pi,acos
# import multiprocessing as mp
# from fiducials import image_dimensions
# from colin import colin



# def return_params(im_name):

#     ''' Final/Correct camera parameters are stored here '''

#     if im_name == "day_ref":
#         params = np. array([1.55288380e+00, -1.08786155e-01, -2.53207783e-02,
#                 9.87892749e+05, 1.91748294e+05, 4.00406708e+02,
#                 -1.63284930e+04])

#     elif im_name == "aps":
#         params = np.array([1.56926535e+00, -1.20789690e-01, \
#                 -3.05255789e-03,9.87920425e+05, 1.91912958e+05, \
#                 3.85333237e+02,-1.10001068e+04])

#     else:
#         print "Not a valid image"
#         return

#     omega, phi, kappa, xs, ys, zs, f = params

#     return omega, phi, kappa, xs, ys, zs, f



# def orientation(omega,phi,f):

#     r31 = sin(phi)
#     r32 = -1.0*sin(omega)*cos(phi)
#     r33 = cos(phi)*cos(omega)

#     orientation = np.array([r31,r32,r33])
#     orientation = np.sign(f)*orientation*(sqrt(orientation[0]**2 +
#                            orientation[1]**2 +
#                            orientation[2]**2)**(-1))

#     return orientation



# def distance(xs,ys,zs,x,y,z):

#     dist = 1.0*((x - xs)**2 + (y - ys)**2 + (z - zs)**2)**0.5

#     return dist



# def in_view(omega,phi,f, xs,ys,zs, x,y,z, view_angle=(pi/2)):
#     # In other words, check if a given point is within the vision cone

#     vector = np.array([(x - xs), (y - ys), (z - zs)])
#     normed_vector = vector / distance(xs,ys,zs, x,y,z)

#     orient = orientation(omega,phi,f)
#     dot = normed_vector[0]*orient[0] + normed_vector[1]*orient[1] + \
#     normed_vector[2]*orient[2]

#     return np.arccos(dot) < view_angle



# def in_picture(x,y,image_dimensions):
#     # Check if point gets mapped to a pixel within the specified x and y
#     # sizes of the image

#     is_in_picture = (x < image_dimensions[0])*(x > 0)*(y > 0)*\
#     (y < image_dimensions[1])

#     return is_in_picture



# def project(filename):

#     # Finds the desired projection

#     if globparams == None:
#         params = return_params(globname)
#     else:
#         params = globparams

#     omega, phi, kappa, xs, ys, zs, f = params
#     image_dims = image_dimensions(globname)
#     image_dims_reversed = np.array([image_dims[1], \
#         image_dims[0]])

#     # Rearrange
#     print "working on: ", filename
#     dat = np.load(filename).T.copy()

#     # Multiply by -1 because it apears as inverse; use orient?
#     pixel_xy = 1.0*colin(params, dat)

#     # un-center pixel (x,y)
#     x = image_dims[0]/2 - pixel_xy[:,0].astype(int)
#     y = image_dims[1]/2 + pixel_xy[:,1].astype(int)

#     is_in_picture = in_picture(x,y,image_dims)

#     index = np.arange(is_in_picture.size)[is_in_picture>0]

#     print "npix = ", index.size

#     distgrid = np.ones(image_dims_reversed)*(100000.0)
#     xgrid =  -1.*np.ones(image_dims_reversed)
#     ygrid = -1.*np.ones(image_dims_reversed)

#     if index.size==0:
#         print "no points, returning..."
#         return [distgrid, xgrid, ygrid]

#     n   = distance(xs,ys,zs, dat[index,0],dat[index,1],dat[index,2])
#     x   = x[index]
#     y   = y[index]
#     dat = dat[index]

#     # Add each point to the arrays, given it is visibile (vis[i] == 1)
#     # And it is closer to the camera than the current value stored in
#     # the corresponding pixel of the distance array

#     nx = distgrid.shape[1]-1
#     ny = distgrid.shape[0]-1

#     for ii in range(index.size):
#         if n[ii]<distgrid[ny-y[ii],nx-x[ii]] and n[ii]>500:
#             distgrid[ny-y[ii],nx-x[ii]] = n[ii]
#             xgrid[ny-y[ii],nx-x[ii]] = dat[ii,0]
#             ygrid[ny-y[ii],nx-x[ii]] = dat[ii,1]

#     print "Done with: ",filename
#     return [distgrid, xgrid, ygrid]



# def merge(final, new):

#     out = [0, 0, 0]
#     replace = np.greater(final[0], new[0])
#     out[0] = final[0]*np.logical_not(replace) + new[0]*replace
#     out[1] = final[1]*np.logical_not(replace) + new[1]*replace
#     out[2] = final[2]*np.logical_not(replace) + new[2]*replace
#     return out


# # Set up the parallel operation
# def project_parallel(file_list,nworkers):

#     print "Initializing workers..."
#     pool = mp.Pool(nworkers)
#     pic_slices = pool.map(project, file_list)
#     pool.close()
#     pool.join()
#     return pic_slices


# def get_file_names(directory):


#     file_list = []

#     all_files = glob.glob1(directory,"*.npy")

#     for file_name in all_files:
#         if len(file_name) == 12:
#             try:
#                 int(file_name[0:8])
#                 file_list.append(file_name)
#             except ValueError:
#                 pass

#     return file_list



# def run_project(name="new",directory="/scratch/share/gdobler/npydat/rawdat/",nworkers=10,params=None):

#     global globname
#     global globparams
#     globname = name
#     globparams = params

#     #name = "day_ref_emp_state"
#     image_dims = image_dimensions(name)
#     image_dims_reversed = np.array([image_dims[1], \
#         image_dims[0]])

#     if params == None:
#         omega, phi, kappa, xs, ys, zs, f = return_params(name)
#     else:
#         omega, phi, kappa, xs, ys, zs, f = params

#     final_grids = [np.ones(image_dims_reversed)*(10**8), \
#         -1*np.ones(image_dims_reversed), \
#         -1*np.ones(image_dims_reversed)]

#     file_list = []

#     #directory = "/home/erc399/data/lidar/"
#     lidar_file_names = get_file_names(directory)

#     num_files = len(lidar_file_names)
#     print 'Total files: ', num_files
#     print "Selecting files..."

#     file_list = [directory+i for i in lidar_file_names]

#     print 'Done selecting files...'
#     print "Files selected: ", len(file_list)


#     # Begin parallel computation of each file's projection
#     print "Begin Multiprocessing..."
#     projected_tiles = project_parallel(file_list,nworkers)


#     # Merge together all outputs into a single file
#     print "Multiprocessing complete..."
#     print "Merging files..."
#     for tile in projected_tiles:
#         final_grids = merge(final_grids, tile)

#     '''
#     # Cascade to fill holes by covering a given point if there is
#     # a closer point above it
#     print "Smoothing vertically..."
#     for i in range(0, len(final_grids[0])-1):
#         for j in range(0, len(final_grids[0][0])):
#             if final_grids[0][i][j] < final_grids[0][i + 1][j]:
#                 final_grids[0][i+1][j] = final_grids[0][i][j]
#                 final_grids[1][i+1][j] = final_grids[1][i][j]
#                 final_grids[2][i+1][j] = final_grids[2][i][j]

#     # Smooth pixels to get rid of vertical bars
#     print "Smoothing horizontally..."
#     for i in range(0, len(final_grids[0])-1):
#         for j in range(1, len(final_grids[0][0])-2):
#             if final_grids[0][i][j] != final_grids[0][i][j-1] and \
#                 final_grids[0][i][j+1] == final_grids[0][i][j+1]:

#                 final_grids[0][i][j] = final_grids[0][i][j+1]
#                 final_grids[1][i][j] = final_grids[1][i][j+1]
#                 final_grids[2][i][j] = final_grids[2][i][j+1]

#     # Smooth pixels two wide
#     print "Smoothing horizontally [2 pixels]..."
#     for i in range(0, len(final_grids[0])-1):
#         for j in range(1, len(final_grids[0][0])-4):
#             if final_grids[0][i][j] != final_grids[0][i][j-1] and \
#                 final_grids[0][i][j] == final_grids[0][i][j+1] and \
#                 final_grids[0][i][j-1] == final_grids[0][i][j+2]:

#                 final_grids[0][i][j] = final_grids[0][i][j+2]
#                 final_grids[1][i][j] = final_grids[1][i][j+2]
#                 final_grids[2][i][j] = final_grids[2][i][j+2]

#                 final_grids[0][i][j+1] = final_grids[0][i][j+2]
#                 final_grids[1][i][j+1] = final_grids[1][i][j+2]
#                 final_grids[2][i][j+1] = final_grids[2][i][j+2]
#     '''

#     # Save data
#     print "Saving data..."
#     np.save("../output/" + name + '_distgrid_raw', final_grids[0])
#     np.save("../output/" + name + '_xgrid_raw', final_grids[1])
#     np.save("../output/" + name + '_ygrid_raw', final_grids[2])

#     print "Done!"
#     return


