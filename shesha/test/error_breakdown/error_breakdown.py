# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 11:00:14 2016

@author: fferreira
"""

import cProfile
import pstats as ps

import sys,os
import numpy as np
import naga as ch
import shesha as ao
import time
import matplotlib.pyplot as pl
import hdf5_utils as h5u
import pandas

if(len(sys.argv) < 2):
    error= 'command line should be at least:"python -i test.py parameters_filename"\n with "parameters_filename" the path to the parameters file'
    raise StandardError(error)


#get parameters from file
param_file=sys.argv[1]
if(param_file.split('.')[-1] == "py"):
    filename=param_file.split('/')[-1]
    param_path=param_file.split(filename)[0]
    sys.path.insert(0,param_path)
    exec("import %s as config" % filename.split(".py")[0])
    sys.path.remove(param_path)
elif(param_file.split('.')[-1] == "h5"):
    sys.path.insert(0,os.environ["SHESHA_ROOT"]+"/data/par/par4bench/")
    import scao_16x16_8pix as config
    sys.path.remove(os.environ["SHESHA_ROOT"]+"/data/par/par4bench/")
    h5u.configFromH5(param_file,config)
else:
    raise ValueError("Parameter file extension must be .py or .h5")

print "param_file is",param_file

if(len(sys.argv) > 2):
    device=int(sys.argv[2])
else:
    device = 0

############################################################################
#  _       _ _       
# (_)_ __ (_) |_ ___ 
# | | '_ \| | __/ __|
# | | | | | | |_\__ \
# |_|_| |_|_|\__|___/
############################################################################                    

if(hasattr(config,"simul_name")):
    if(config.simul_name is None):
        simul_name=""
    else:
        simul_name=config.simul_name
else:
    simul_name=""
print "simul name is",simul_name

matricesToLoad={}
if(simul_name==""):
    clean=1
else:
    clean=0
    param_dict = h5u.params_dictionary(config)
    matricesToLoad = h5u.checkMatricesDataBase(os.environ["SHESHA_ROOT"]+"/data/",config,param_dict)
#initialisation:
#   context
c=ch.naga_context()
c.set_activeDevice(device)

#    wfs
print "->wfs"
wfs,tel=ao.wfs_init(config.p_wfss,config.p_atmos,config.p_tel,config.p_geom,config.p_target,config.p_loop, 1,0,config.p_dms)

#   atmos
print "->atmos"
atm=ao.atmos_init(c,config.p_atmos,config.p_tel,config.p_geom,config.p_loop,config.p_wfss,config.p_target,rank=0, clean=clean, load=matricesToLoad)

#   dm
print "->dm"
dms=ao.dm_init(config.p_dms,config.p_wfss,config.p_geom,config.p_tel)

#   target
print "->target"
tar=ao.target_init(c,tel,config.p_target,config.p_atmos,config.p_geom,config.p_tel,config.p_wfss,wfs,config.p_dms)

print "->rtc"
#   rtc
rtc=ao.rtc_init(tel,wfs,config.p_wfss,dms,config.p_dms,config.p_geom,config.p_rtc,config.p_atmos,atm,config.p_tel,config.p_loop,tar,config.p_target,clean=clean,simul_name=simul_name, load=matricesToLoad)

if not clean:
    h5u.validDataBase(os.environ["SHESHA_ROOT"]+"/data/",matricesToLoad)

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

error_flag = True in [w.error_budget for w in config.p_wfss]

##################################################################################################### 
#   __                  _   _                 
#  / _|_   _ _ __   ___| |_(_) ___  _ __  ___ 
# | |_| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# |  _| |_| | | | | (__| |_| | (_) | | | \__ \
# |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
###################################################################################################
def error_breakdown(com,noise_com,alias_com,tomo_com,H_com,trunc_com,bp_com,wf_com,i):
    """
    Compute the error breakdown of the AO simulation. Returns the error commands of 
    each contributors. Suppose no delay (for now) and only 2 controllers : the main one, controller #0, (specified on the parameter file) 
    and the geometric one, controller #1 (automatically added if error_budget is asked in the parameter file)
    Commands are computed by applying the loop filter on various kind of slopes/commands : (see schema_simulation_budget_erreur_v2) 
    
        - A : Aliasing contribution 
            Obtained by computing commands from DM orthogonal phase (projection + slopes_geom)
            
        - B : Wavefront + Filtered modes
            Obtained as the commmands output of the geometric controller
        
        - C : Wavefront 
            Obtained by computing commands from DM parallel phase (projection + slopes_geom)
        
        - D : Wavefont + aliasing + ech/truncature + noise + tomo
            Obtained as the classical AO loop commands (just do a get_commands at the end of an iteration)
        
        - E : Wavefront + aliasing + ech/trunc + tomo
            Obtained by performing the AO loop iteration without noise on the WFS
        
        - F : Wavefront + aliasing + tomo
            Obtained by performing the AO loop iteration without noise on the WFS and using phase deriving slopes
    
        - G : tomo
    
    Note : rtc.getErr returns to -CMAT.slopes
        
    :parameters:
        noise_com : np.array((nactu,niter)) : Noise contribution
            Computed with D-E
            
        alias_com : np.array((nactu,niter)) : Aliasing contribution
            Computed with A
            
        tomo_com : np.array((nactu,niter)) : Tomographic error contribution
            Computed with C-B
            
        H_com : np.array((nactu,niter)) : Filtered modes error
            Computed with B
        
        trunc_com : np.array((nactu,niter)) : sampling/truncature error contribution
            Computed with E-F
            
        bp_com : np.array((nactu,niter)) : Bandwidth error
        
        wf_com : np.array((nactu,niter)) : Reconstructed wavefront 
        
        i : (int) : current iteration number
            
    """
    g = config.p_controllers[0].gain
#    g=1.
    D = rtc.getCentroids(0) # Classical slopes with all contributions
    Dcom = rtc.getCom(0)
    com[:,i] = Dcom
    tarphase = tar.get_phase(0)
    ###########################################################################
    ## Noise contribution
    ###########################################################################
    ideal_bincube = wfs.get_bincubeNotNoisy(0)
    bincube = wfs.get_bincube(0)
    if(config.p_centroiders[0].type_centro == "tcog"): # Select the same pixels with or without noise
        invalidpix = np.where(bincube <= config.p_centroiders[0].thresh)
        ideal_bincube[invalidpix] = 0
        rtc.setthresh(0,-1e16)
    wfs.set_bincube(0,ideal_bincube)
    rtc.docentroids(0)
    if(config.p_centroiders[0].type_centro == "tcog"):
        rtc.setthresh(0,config.p_centroiders[0].thresh)
    
    E = rtc.getCentroids(0) #Slopes without noise
    
    rtc.setCentroids(0,D-E) # D-E = noise on slopes
    rtc.docontrol(0)
    # Apply loop filter to get contribution of noise on commands
    noise_com[:,i] = noise_com[:,i-1]*(1-g) + g*rtc.getErr(0)
    
    ###########################################################################
    ## Sampling/truncature contribution
    ###########################################################################
    rtc.docentroids_geom(0)
    F = rtc.getCentroids(0) # Slopes without noise and sampling/truncature error
    rtc.setCentroids(0,E-F) # E-F = sampling/truncature error on slopes
    rtc.docontrol(0)
    # Apply loop filter to get contribution of sampling/truncature on commands
    trunc_com[:,i] = trunc_com[:,i-1]*(1-g) + g*rtc.getErr(0)

    ###########################################################################
    ## Wavefront + filtered modes reconstruction
    ###########################################################################
    tar.atmos_trace(0,atm,tel)
    rtc.docontrol_geo(1, dms, tar, 0)
    B = rtc.getCom(1)

    ###########################################################################
    ## Wavefront
    ###########################################################################  
    rtc.applycontrol(1,dms)
    for w in range(len(config.p_wfss)):
        wfs.sensors_trace(w,"dm",tel,atm,dms,rst=1)
        wfs.sensors_compimg(w)
    ideal_bincube = wfs.get_bincubeNotNoisy(0)
    wfs.set_bincube(0,ideal_bincube)
    rtc.docentroids_geom(0)
    rtc.docontrol(0)
    C = (-1.) * rtc.getErr(0)
    
    ###########################################################################
    ## Filtered modes error
    ###########################################################################    
    H_com[:,i] = B - C
    
    ###########################################################################
    ## Bandwidth error
    ###########################################################################
    wf_com[:,i] = wf_com[:,i-1]*(1-g) + g*C
    bp_com[:,i] = C - wf_com[:,i-1]
    
    ###########################################################################
    ## Aliasing error
    ###########################################################################       
#    tar.dmtrace(0,dms,do_phase_var=0)
    for w in range(len(config.p_wfss)):
        wfs.sensors_trace(w,"all",tel,atm,dms)
        wfs.sensors_compimg(w)
    ideal_bincube = wfs.get_bincubeNotNoisy(0)
    wfs.set_bincube(0,ideal_bincube)
    rtc.docentroids_geom(0)
    rtc.docontrol(0)
    A = rtc.getErr(0)
    alias_com[:,i] = alias_com[:,i-1]*(1-g) + g*A
    
    ###########################################################################
    ## Tomographic error
    ###########################################################################
    rtc.setCentroids(0,F)
    rtc.docontrol(0)
    G = rtc.getErr(0) - (C + A - com[:,i-1])
#    Fc = rtc.getErr(0)
#    rtc.setCom(0,com[:,i-1].astype(np.float32))
#    rtc.applycontrol(0,dms)
#    
#
#    for w in range(len(config.p_wfss)):
#        wfs.sensors_trace(w,"dm",tel,atm,dms,rst=1)
#        wfs.sensors_compimg(w)
#    rtc.docentroids_geom(0)
#    rtc.docontrol(0)
#    
#    G = Fc - (A+C+rtc.getErr(0))
    
    tomo_com[:,i] = tomo_com[:,i-1]*(1-g) + g*G

#    if(i>0):
#        raise ValueError
   
    # Without anyone noticing...
    tar.set_phase(0,tarphase)
    rtc.setCom(0,Dcom) 

    
def loop(n):
    """
    Performs the main AO loop for n interations
    
    :param n: (int) : number of iterations
    """
    if(error_flag):
    # Initialize buffers for error breakdown
        nactu = rtc.getCom(0).size
        com = np.zeros((nactu,n))
        noise_com = np.zeros((nactu,n))
        alias_com = np.copy(noise_com)
        wf_com = np.copy(noise_com)
        tomo_com = np.copy(noise_com)
        trunc_com = np.copy(noise_com)
        H_com = np.copy(noise_com)
        bp_com = np.copy(noise_com)
        
    t0=time.time()
    for i in range(n):
        atm.move_atmos()

        if(config.p_controllers[0].type_control == "geo"):
            for t in range(config.p_target.ntargets):
                tar.atmos_trace(t,atm,tel)
                rtc.docontrol_geo(0, dms, tar, 0)
                rtc.applycontrol(0,dms)
                tar.dmtrace(0,dms)
        else:
            for t in range(config.p_target.ntargets):
                tar.atmos_trace(t,atm,tel)
                tar.dmtrace(t,dms)
            for w in range(len(config.p_wfss)):
                wfs.sensors_trace(w,"all",tel,atm,dms)
                wfs.sensors_compimg(w)
            rtc.docentroids(0)
            rtc.docontrol(0)
            
            if(error_flag):
            #compute the error breakdown for this iteration
                error_breakdown(com,noise_com,alias_com,tomo_com,H_com,trunc_com,bp_com,wf_com,i)
            
            rtc.applycontrol(0,dms)
                
        if((i+1)%100==0):
            strehltmp = tar.get_strehl(0)
            print i+1,"\t",strehltmp[0],"\t",strehltmp[1]
    t1=time.time()
    print " loop execution time:",t1-t0,"  (",n,"iterations), ",(t1-t0)/n,"(mean)  ", n/(t1-t0),"Hz"
    if(error_flag):
    #Returns the error breakdown
        return com,noise_com,alias_com,tomo_com,H_com,trunc_com,bp_com,wf_com

###############################################################################################
#  _            _       
# | |_ ___  ___| |_ ___ 
# | __/ _ \/ __| __/ __|
# | ||  __/\__ \ |_\__ \
#  \__\___||___/\__|___/
###############################################################################################                       
G = rtc.get_imat(0)*0.9679
rtc.set_imat(0,G)
ao.cmat_init(0,rtc,config.p_rtc,config.p_wfss,config.p_atmos,config.p_tel, config.p_dms, clean=1)
com,noise_com,alias_com,tomo_com,H_com,trunc_com,bp_com,wf_com = loop(1000)
filename = "/home/fferreira/Data/breakdown_imat_hack.h5"
h5u.writeHdf5SingleDataset(filename,com,datasetName="com")
h5u.save_hdf5(filename,"noise_com",noise_com)
h5u.save_hdf5(filename,"alias_com",alias_com)
h5u.save_hdf5(filename,"tomo_com",tomo_com)
h5u.save_hdf5(filename,"H_com",H_com)
h5u.save_hdf5(filename,"trunc_com",trunc_com)
h5u.save_hdf5(filename,"bp_com",bp_com)
h5u.save_hdf5(filename,"wf_com",wf_com)

pl.ion()
pl.plot(np.var(noise_com,axis=1))
pl.plot(np.var(alias_com,axis=1))
pl.plot(np.var(tomo_com,axis=1))
pl.plot(np.var(H_com,axis=1))
pl.plot(np.var(trunc_com,axis=1))
pl.plot(np.var(bp_com,axis=1))
pl.plot(np.var(wf_com,axis=1))
pl.plot(np.var(com,axis=1))
pl.legend(["noise","aliasing","tomo","H","trunc","BP","wf","com"])
                                             