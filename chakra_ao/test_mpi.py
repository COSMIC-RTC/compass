import os

import cProfile
import pstats as ps

import sys
import numpy as np
import chakra as ch
import chakra_ao as ao
import time

rank=int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK'])
c=ch.chakra_context()
c.set_activeDevice(rank%c.get_ndevice())

#Delay import because of cuda_aware
#mpi_init called during the import
import mpi4py
from mpi4py import MPI



comm=MPI.COMM_WORLD
comm_size=comm.Get_size()
rank=comm.Get_rank()



if(len(sys.argv)!=2):
    error= 'command line should be:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
    raise StandardError(error)

#get parameters from file
param_file=sys.argv[1]
execfile(param_file)

#initialisation:
#   wfs
if(rank==0):print "->wfs"
wfs=ao.wfs_init(p_wfss,p_atmos,p_tel,p_geom,p_target,p_loop,comm_size,rank,p_dms)

#   atmos
if(rank==0):print "->atmos"
atm=p_atmos.atmos_init(c,p_tel,p_geom,p_loop,rank=rank)

#   dm 
if(rank==0):print "->dm"
dms=ao.dm_init(p_dms,p_wfs0,p_geom,p_tel)

#   target
if(rank==0):print "->target"
tar=p_target.target_init(c,p_atmos,p_geom,p_tel,p_wfss,wfs,p_dms)

#   rtc
if(rank==0):print "->rtc"
rtc=ao.rtc_init(wfs,p_wfss,dms,p_dms,p_geom,p_rtc,p_atmos,atm,p_loop,p_target)
comm.Barrier()


if(rank==0):
    print atm
    print wfs
    print dms
    print tar
    print rtc
comm.Barrier()

start = time.time()
#generate turbulence, raytrace through the atmos
# and display:
#   the turbulence
#   the target's image 
ao.see_atmos_target_disp_mpi(10,atm,tar,wfs,comm,alt=0,n_tar=0,f=0.15)

#for profiling purpose
#filename="Profile"+str(rank)+".prof"
#cProfile.runctx("ao.see_atmos_target_mpi(100,atm,tar,wfs,comm,alt=0,n_tar=0,f=0.15)", globals(), locals(),filename)

end = time.time()

for i in range(comm_size):
    if(rank==i):
        #for profiling purpose
        #s = ps.Stats(filename)
        #s.strip_dirs().sort_stats("time").print_stats()
        print end - start
    comm.Barrier()


#use to compare result with non mpi exec
# load numpy array with numpy.load("wfs_binimg_mpi.npy")
# for wfs_binimg_mpi.npy, data for process >0 are not send to process 0
# => compare only the data from process 0
"""
if(rank==0):
    shak=wfs.get_binimg(0)
    screen=atm.get_screen(0)
    img=tar.get_image(0,"se")

    np.save("wfs_binimg_mpi",shak)
    np.save("target_image_mpi",img)
    np.save("atmos_screen_mpi",screen)
"""
