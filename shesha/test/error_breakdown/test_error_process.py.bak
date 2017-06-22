# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 08:53:38 2015

@author: fferreira
"""

import cProfile
import pstats as ps

import sys, os
import numpy as np
import naga as ch
import shesha as ao
import time
import matplotlib.pyplot as pl

import tools2 as tools
import hdf5_utils as h5u
import pandas

#get parameters from file
param_file= os.environ["SHESHA_ROOT"]+"/data/par/noise_error_process_80x80_linearity.py"
filename=param_file.split('/')[-1]
param_path=param_file.split(filename)[0]
sys.path.insert(0,param_path)
exec("import %s as config" % filename.split(".py")[0])

print "param_file is",param_file

if(hasattr(config,"simul_name")):
    if(config.simul_name is None):
        simul_name=""
    else:
        simul_name=config.simul_name
else:
    simul_name=""

if(simul_name==""):
    clean=1
else:
    clean=0

    
#clean = 1

#   context
c=ch.naga_context()

    
    
def Init(config,device=0):
    #initialisation:
    c.set_activeDevice(device)
    matricesToLoad={}
    if(clean == 0):
        param_dict = h5u.params_dictionary(config)
        matricesToLoad = h5u.checkMatricesDataBase(os.environ["SHESHA_ROOT"]+"/data/",config,param_dict)
    
    #    wfs
    print "->wfs"
    wfs, tel=ao.wfs_init(config.p_wfss,config.p_atmos,config.p_tel,config.p_geom,config.p_target,config.p_loop, 1,0,config.p_dms)
    
    #   atmos
    print "->atmos"
    atm=ao.atmos_init(c,config.p_atmos,config.p_tel,config.p_geom,config.p_loop,config.p_wfss,wfs,config.p_target,clean=clean,load=matricesToLoad,rank=0)
    
    #   dm 
    print "->dm"
    dms=ao.dm_init(config.p_dms,config.p_wfss,config.p_geom,config.p_tel)
    
    #   target
    print "->target"
    tar=ao.target_init(c,tel,config.p_target,config.p_atmos,config.p_geom,config.p_tel,config.p_dms)
    
    print "->rtc"
    #   rtc
    rtc=ao.rtc_init(tel,wfs,config.p_wfss,dms,config.p_dms,config.p_geom,config.p_rtc,config.p_atmos,atm,config.p_tel,config.p_loop,tar,config.p_target,clean=clean,simul_name=simul_name,load=matricesToLoad)
    
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
    
    return atm, wfs, dms, tar, rtc, tel

#############################################################################################################
################################################# MAIN ANALYSIS #############################################
#############################################################################################################
#Allocate buffers
def check_param(config,niter,save=False,device=5):
    
    param = [40,50,60,70,77]
    #param = [0.13,0.14,0.15,0.16]
    brvar = []
    nvar = []
    dsp_med = []
    f_med = []
    dsp_Kmed = []
    f_Kmed = []
    if(save):
        store = pandas.HDFStore("/home/fferreira/Data/resErrorProcessMicado2.h5")
        df=pandas.DataFrame(columns=["bisimus","monosimu","dspm","fm","dspKmed","fKmed","config_mono","config_bi","nm_rms_noise","SRnoise","SR"],dtype=object)
        store.put("results",df)
        store.close()
    for i in range(len(param)):
        bv, nv, dspm, fm, dspKmed, fKmed=process_err(config,niter,param[i],save=save,linear=True,device=device)
        if(save):
            store = pandas.HDFStore("/home/fferreira/Data/resErrorProcessMicado2.h5")
            df=store["results"]
            df.loc[param[i],"bisimus"]=bv
            df.loc[param[i],"monosimu"]=nv
            df.loc[param[i],"dspm"]=dspm
            df.loc[param[i],"fm"]=fm
            df.loc[param[i],"dspKmed"]=dspKmed
            df.loc[param[i],"fKmed"]=fKmed
            store.put("results",df)
            store.close()

        brvar.append(bv)
        nvar.append(nv)
        dsp_med.append(dspm)
        f_med.append(fm)
        dsp_Kmed.append(dspKmed)
        f_Kmed.append(fKmed)
        
    #store.close()
        
    
    return param,brvar,nvar,dsp_med,f_med,dsp_Kmed,f_Kmed
    
    
def process_err(config,niter,param,save=False,linear=False,device=0):
    
    config = reload(config)
    config.p_wfss[0].set_nxsub(param)
    config.p_dms[0].set_nact(param+1)
    config.p_wfss[0].set_noise(-1)
    if(linear):
        FoV = 2.
        RASC = 180/np.pi * 3600.
        d = config.p_tel.diam/config.p_wfss[0].nxsub
        pixsize = config.p_wfss[0].Lambda*1e-6/d * RASC / 2.5 #bit better than Shannon
        config.p_wfss[0].set_pixsize(pixsize)
        config.p_wfss[0].set_npix(long(FoV/pixsize))
    if(save):
        store=pandas.HDFStore("/home/fferreira/Data/resErrorProcessMicado2.h5")
        #config.simul_name = "errorProcess_27Jan2016"
        pdict = h5u.params_dictionary(config)
        #config.simul_name = ""
        
        df=store["results"]
        df.loc[param,"config_mono"]=[pdict]
        
    atm, wfs, dms, tar, rtc, tel =Init(config,device=device)
    nactu = rtc.getCom(0).size
    nslope = rtc.getCentroids(0).size
    idcom = np.zeros((nactu,niter),dtype=np.float32)
    full_com = np.zeros((nactu,niter),dtype=np.float32)
    noise_com = np.zeros((nactu,niter),dtype=np.float32)
    noise_mes = np.zeros((nslope,niter),dtype=np.float32)
    
    #cmat = np.load("cmat.npy")
    #rtc.set_cmat(0,cmat)
    loop(niter,idcom,noise_com,noise_mes,atm, wfs, dms, tar, rtc, tel, config)
    if(save):
        df.loc[param,"SR"] = tar.get_strehl(0)[1]


    #config.p_wfss[0].set_noise(noise)
    config = reload(config) 
    config.p_wfss[0].set_nxsub(param)
    config.p_dms[0].set_nact(param+1)
    S = np.pi/4.*(1-config.p_tel.cobs**2)*config.p_tel.diam**2
    Nph = 300.*(40./param)
    m = Nph / (config.p_wfss[0].zerop * config.p_wfss[0].optthroughput * config.p_loop.ittime * (config.p_tel.diam/config.p_wfss[0].nxsub)**2 ) * S
    m = np.log10(m)/0.4 * (-1)
    #config.p_wfss[0].set_gsmag(m)
    if(linear):
        FoV = 2.
        RASC = 180/np.pi * 3600.
        d = config.p_tel.diam/config.p_wfss[0].nxsub
        pixsize = config.p_wfss[0].Lambda*1e-6/d * RASC / 2.5 #bit better than Shannon
        config.p_wfss[0].set_pixsize(pixsize)
        config.p_wfss[0].set_npix(long(FoV/pixsize))

    if(save):
        #config.simul_name = "errorProcess_Micado"
        pdict = h5u.params_dictionary(config)
        #config.simul_name = ""
        df.loc[param,"config_bi"]=[pdict]


    

    
    atm, wfs, dms, tar, rtc, tel =Init(config,device=device)
    nactu = rtc.getCom(0).size
    nslope = rtc.getCentroids(0).size   
    full_com = np.zeros((nactu,niter),dtype=np.float32)
    noise_com = np.zeros((nactu,niter),dtype=np.float32)
    noise_mes = np.zeros((nslope,niter),dtype=np.float32)
    loop(niter,full_com,noise_com,noise_mes,atm, wfs, dms, tar, rtc, tel, config)
    if(save):
        df.loc[param,"SRnoise"] = tar.get_strehl(0)[1]
    
    
    RASC = 180/np.pi * 3600.
    sspdiam = config.p_tel.diam / config.p_wfss[0].nxsub
    
    noise_mes = noise_mes * sspdiam / RASC * 1e9
    
    print " "
    print "Mean noise : ",np.mean(np.std(noise_mes,axis=1)), "nm rms"
    
    if(save):
        df.loc[param,"nm_rms_noise"]=np.mean(np.std(noise_mes,axis=1)) 
        store.put("results",df)
        store.close()
    
    #idcom = np.load("idcom_ol.npy")
    brcom = full_com - idcom
    brvar = np.var(brcom,axis=1)
    nvar = np.var(noise_com,axis=1)
    
    err = brcom - noise_com
    #dsp,f = tools.DSP(err,1.,1./config.p_loop.ittime,1)
    dsp_med,f_med = tools.DSP(err,1.,1./config.p_loop.ittime,1,med=True)
    dsp_Kmed, f_Kmed = tools.DSP(err,1.,1./config.p_loop.ittime,1,K=128,med=True)
    #dsp_K, f_K = tools.DSP(err,1.,1./config.p_loop.ittime,1,K=128)
    
#    pl.ion()
#    pl.plot(brvar,color="blue")
#    pl.plot(nvar,color="red")
#    
#    pl.figure(3)
#    pl.plot(f_med,dsp_med)
#    pl.plot(f_Kmed,dsp_Kmed)
    
    #return idcom, full_com, noise_com, dsp_med, f_med
    return brvar, nvar, dsp_med, f_med, dsp_Kmed, f_Kmed


def loop(n,full_com,noise_com,noise_mes,atm, wfs, dms, tar, rtc, tel,config):
    """ Here, we want to get centroids obtain from a SH wfs seeing an atmos screen
    and compare it to the centroids obtain from a wfs seeing the parallel component
    and a wfs seeing the orthogonal component
    For this, 1 wfs is seeing the atmos in open loop, an other sees atmos+dm after
    it shape has been computed with a geo controller, and the last one just sees
    the dm
    """
    reset_com = rtc.getCom(0)*0.
    for i in range(n):
        atm.move_atmos()
        tar.atmos_trace(0,atm,tel)
        tar.dmtrace(0,dms)
        
        wfs.sensors_trace(0,"all",tel,atm,dms)
        wfs.sensors_compimg(0)
        ideal_bincube = wfs.get_bincubeNotNoisy(0)
        bincube = wfs.get_bincube(0)
        if(config.p_centroiders[0].type_centro == "tcog"):
            invalidpix = np.where(bincube <= config.p_centroiders[0].thresh)
        #noise_bincube = bincube-ideal_bincube
        
        #compute command normally
        rtc.docentroids(0)
        mesb = rtc.getCentroids(0)
        rtc.docontrol(0)
        full_com[:,i] = rtc.getCom(0)
        rtc.applycontrol(0,dms)
        current_com = rtc.getCom(0)
        
        #compute noise command
        rtc.setCom(0,reset_com)
        if(config.p_centroiders[0].type_centro == "tcog"):
            ideal_bincube[invalidpix] = 0
            rtc.setthresh(0,-1e16)
        wfs.set_bincube(0,ideal_bincube)
        rtc.docentroids(0)
        rtc.setthresh(0,config.p_centroiders[0].thresh)
        mes = rtc.getCentroids(0)
        noise = mesb - mes
        noise_mes[:,i] = noise
        rtc.setCentroids(0,noise)
        rtc.docontrol(0)
        noise_com[:,i] = noise_com[:,i-1]*(1-config.p_controllers[0].gain) + rtc.getCom(0) #rtc.getErr(0)*config.p_controllers[0].gain
        rtc.setCom(0,current_com)
        
        if((i+1)%1000==0):
            strehltmp = tar.get_strehl(0)
            print i+1,"\t",strehltmp[0],"\t",strehltmp[1]
        
        print "\r Recording... %d%%"%(i*100/n),
        
def openloop(n,full_com,noise_com,noise_mes, atm, wfs, dms, tar, rtc, tel,config):
    """ Here, we want to get centroids obtain from a SH wfs seeing an atmos screen
    and compare it to the centroids obtain from a wfs seeing the parallel component
    and a wfs seeing the orthogonal component
    For this, 1 wfs is seeing the atmos in open loop, an other sees atmos+dm after
    it shape has been computed with a geo controller, and the last one just sees
    the dm
    """
    reset_com = rtc.getCom(0)*0.
    for i in range(n):
        atm.move_atmos()
        tar.atmos_trace(0,atm,tel)
        wfs.sensors_trace(0,"atmos",tel,atm)
        wfs.sensors_compimg(0)
        ideal_bincube = wfs.get_bincubeNotNoisy(0)
        bincube = wfs.get_bincube(0)
        #noise_bincube = bincube-ideal_bincube
        
        #compute command normally
        rtc.docentroids(0)
        mesb = rtc.getCentroids(0)
        rtc.docontrol(0)
        full_com[:,i] = rtc.getCom(0)
        rtc.applycontrol(0,dms)
        rtc.setCom(0,reset_com)
        tar.dmtrace(0,dms)
        
        #compute noise command
        rtc.setCom(0,reset_com)
        wfs.set_bincube(0,ideal_bincube)
        rtc.docentroids(0)
        mes = rtc.getCentroids(0)
        noise = mesb - mes
        noise_mes[:,i] = noise
        rtc.setCentroids(0,noise)
        rtc.docontrol(0)
        noise_com[:,i] = rtc.getCom(0)
        rtc.setCom(0,reset_com)
        
        if((i+1)%1000==0):
            strehltmp = tar.get_strehl(0)
            print i+1,"\t",strehltmp[0],"\t",strehltmp[1]
        
        print "\r Recording... %d%%"%(i*100/n),    








####################################### OLD VERSIONS ###########################################################################################

###############################################################################################################################################################
# Here, we want to get centroids obtain from a SH wfs seeing an atmos screen
# and compare it to the centroids obtain from a wfs seeing the parallel component
# and a wfs seeing the orthogonal component
# For this, 1 wfs is seeing the atmos in open loop, an other sees atmos+dm after
# it shape has been computed with a geo controller, and the last one just sees
# the dm
###############################################################################################################################################################



    
def simple_loop(niter,gain,amp):
    tur = np.ones(niter)*amp # Turbulence
    com_id = np.zeros(niter) # Commande idéale, sans bruit
    mes_id = np.zeros(niter) # Mesure idéale, sans bruit
    com = np.zeros(niter) # Commande réelle
    mes = np.zeros(niter) # Mesure réelle avant ajout du bruit
    com_br = np.zeros(niter) # Commande réelle bruitée
    mes_br = np.zeros(niter) # Mesure réelle bruitée
    b = np.random.randn(niter) # Bruit
    brc = np.zeros(niter)
    err = np.zeros(niter) # Erreur due au bruit
    
    for i in range(niter):
        # Mesure classique
        mes_id[i] = tur[i] + com_id[i-1]
        mes_br[i] = tur[i] + com_br[i-1] + b[i]
        mes[i] = tur[i] + com_br[i-1]
        #Commande classique
        com_br[i] = com_br[i-1] - gain * mes_br[i]
        #Commande hybride
        com[i] = com[i-1] - gain * mes[i]
        # Commande "idéale"
        com_id[i] = com_id[i-1] - gain * mes_id[i]
        #Bruit dans la loop
        brc[i] = (1-gain)*brc[i-1] - gain*b[i]
        
        #Erreur due au bruit
        err[i] = brc[i] + com_id[i]
    
    pl.ion()
    pl.clf()
    pl.plot(com_br,color="red")
    pl.plot(com_id,color="blue")
    pl.plot(brc,color="green")
    pl.plot(com_br-com,color="blue")
    pl.plot(brc-(com_br-com_id),color="black")

#######################################################################
#  _           _       _     
# | |__   __ _| |_ ___| |__  
# | '_ \ / _` | __/ __| '_ \ 
# | |_) | (_| | || (__| | | |
# |_.__/ \__,_|\__\___|_| |_|
##########################################################################  
param,brvar,nvar,dsp_med,f_med,dsp_Kmed,f_Kmed = check_param(config,50000,save=True)




