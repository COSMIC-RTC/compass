import os

import cProfile
import pstats as ps

import sys
import numpy as np
import chakra as ch
import chakra_ao as ao
import time

rank=int(os.environ['OMPI_COMM_WORLD_RANK'])
c=ch.chakra_context()
c.set_activeDevice(rank%c.get_ndevice())

#Delay import because of cuda_aware
#mpi_init called during the import
import mpi4py
from mpi4py import MPI



comm=MPI.COMM_WORLD
comm_size=comm.Get_size()
rank=comm.Get_rank()

print "TEST CHAKRA_AO\n closed loop with MPI"

if(len(sys.argv)!=2):
    error= 'command line should be:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
    raise StandardError(error)

#get parameters from file
param_file=sys.argv[1]
execfile(param_file)

#initialisation:
#    wfs
print "->wfs"
wfs=ao.wfs_init(p_wfss,p_atmos,p_tel,p_geom,p_target,p_loop, comm_size,rank,p_dms)

#   atmos
print "->atmos"
atm=p_atmos.atmos_init(c,p_tel,p_geom,p_loop,rank=rank)

#   dm 
print "->dm"
dms=ao.dm_init(p_dms,p_wfs0,p_geom,p_tel)

#   target
print "->target"
tar=p_target.target_init(c,p_atmos,p_geom,p_tel,p_wfss,wfs,p_dms)

#   rtc
print "->rtc"
rtc=ao.rtc_init(wfs,p_wfss,dms,p_dms,p_geom,p_rtc,p_atmos,atm,p_loop,p_target,simul_name="simName")

comm.Barrier()
if(rank==0):
    print "===================="
    print "init done"
    print "===================="
    print "objects initialzed on GPU:"
    print "--------------------------------------------------------"
    print atm
    print wfs
    print dms
    print tar
    print rtc

    print "----------------------------------------------------";
    print "iter# | S.E. SR | L.E. SR | Est. Rem. | framerate";
    print "----------------------------------------------------";
comm.Barrier()

mimg = 0.# initializing average image

#import matplotlib.pyplot as pl


def loop( n):
    #if(rank==0):
    #    fig,((turbu,image),(shak,defMir))=pl.subplots(2,2, figsize=(15,15))
    #    pl.ion()
    #    pl.show()

    t0=time.time()
    for i in range(n):
        if(rank==0):
            atm.move_atmos()

            for t in range(p_target.ntargets):
                tar.atmos_trace(t,atm)
                tar.dmtrace(t,dms)
            for w in range(len(p_wfss)):
                wfs.sensors_trace(w,"all",atm,dms)

        wfs.Bcast_dscreen()
        for w in range(len(p_wfss)):
            wfs.sensors_compimg(w)
            wfs.gather_bincube(comm,w)
        if(rank==0):
            rtc.docentroids(0)
            rtc.docontrol(0)
            rtc.applycontrol(0,dms)

        if((i+1)%50==0):
            #s=rtc.getCentroids(0)
            if(rank==0):
                """ FOR DEBUG PURPOSE
                #turbu.clear()
                #image.clear()
                #shak.clear()
                #defMir.clear()

                screen=atm.get_screen(0.)

                im=tar.get_image(0,"se")
                im=np.roll(im,im.shape[0]/2,axis=0)
                im=np.roll(im,im.shape[1]/2,axis=1)

                sh=wfs.get_binimg(0)
                
                dm=dms.get_dm("pzt",0.)

                #f1=turbu.matshow(screen,cmap='Blues_r')
                #f2=image.matshow(im,cmap='Blues_r')
                #f3=shak.matshow(sh,cmap='Blues_r')
                #f4=defMir.matshow(dm)
                #pl.draw()

                c=rtc.getCom(0)
                v=rtc.getVoltage(0)

                sh_file="dbg/shak_"+str(i)+"_np_"+str(comm.Get_size())+".npy"
                im_file="dbg/imag_"+str(i)+"_np_"+str(comm.Get_size())+".npy"
                dm_file="dbg/DM_"+str(i)+"_np_"+str(comm.Get_size())+".npy"
                s_file="dbg/cent_"+str(i)+"_np_"+str(comm.Get_size())+".npy"
                c_file="dbg/comm_"+str(i)+"_np_"+str(comm.Get_size())+".npy"
                v_file="dbg/volt_"+str(i)+"_np_"+str(comm.Get_size())+".npy"

                np.save(sh_file,sh)
                np.save(im_file,im)
                np.save(dm_file,dm)
                np.save(s_file,s)
                np.save(c_file,c)
                np.save(v_file,v)
                """

               
                strehltmp = tar.get_strehl(0)
                print "%5d"%(i+1),"  %1.5f"%strehltmp[0],"  %1.5f"%strehltmp[1]

        t1=time.time()
        print rank, "| loop execution time:",t1-t0,"  (",n,"iterations), ",(t1-t0)/n,"(mean)


loop(p_loop.niter)
