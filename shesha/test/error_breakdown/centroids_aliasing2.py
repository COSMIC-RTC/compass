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

#get parameters from file
param_file= os.environ["SHESHA"]+"/data/par/aliasing2.py"
filename=param_file.split('/')[-1]
param_path=param_file.split(filename)[0]
sys.path.insert(0,param_path)
exec("import %s as config" % filename.split(".py")[0])
sys.path.remove(param_path)

print("param_file is",param_file)

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

#initialisation:
#   context
c=ch.naga_context()
c.set_activeDevice(0)

#    wfs
print("->wfs")
wfs=ao.wfs_init(config.p_wfss,config.p_atmos,config.p_tel,config.p_geom,config.p_target,config.p_loop, 1,0,config.p_dms)

#   atmos
print("->atmos")
atm=ao.atmos_init(c,config.p_atmos,config.p_tel,config.p_geom,config.p_loop,config.p_wfss,wfs,config.p_target,rank=0)

#   dm 
print("->dm")
dms=ao.dm_init(config.p_dms,config.p_wfss[0],config.p_geom,config.p_tel)

#   target
print("->target")
tar=ao.target_init(c,config.p_target,config.p_atmos,config.p_geom,config.p_tel,wfs,config.p_dms)

print("->rtc")
#   rtc
rtc=ao.rtc_init(wfs,config.p_wfss,dms,config.p_dms,config.p_geom,config.p_rtc,config.p_atmos,atm,config.p_tel,config.p_loop,tar,config.p_target,clean=clean,simul_name=simul_name)

print("====================")
print("init done")
print("====================")
print("objects initialzed on GPU:")
print("--------------------------------------------------------")
print(atm)
print(wfs)
print(dms)
print(tar)
print(rtc)

print("----------------------------------------------------");
print("iter# | S.E. SR | L.E. SR | Est. Rem. | framerate");
print("----------------------------------------------------");

###############################################################################################################################################################
# Here, we want to get centroids obtain from a SH wfs seeing an atmos screen
# and compare it to the centroids obtain from a wfs seeing the parallel component
# and a wfs seeing the orthogonal component
# For this, 1 wfs is seeing the atmos in open loop, an other sees atmos+dm after
# it shape has been computed with a geo controller, and the last one just sees
# the dm
###############################################################################################################################################################

def loop(n,full_slopes,ortho_slopes,parallel_slopes):
    """ Here, we want to get centroids obtain from a SH wfs seeing an atmos screen
    and compare it to the centroids obtain from a wfs seeing the parallel component
    and a wfs seeing the orthogonal component
    For this, 1 wfs is seeing the atmos in open loop, an other sees atmos+dm after
    it shape has been computed with a geo controller, and the last one just sees
    the dm
    """
    pl.ion()
    reset_com = np.zeros(rtc.getCom(0).size,dtype=np.float32)
    for i in range(n):
        for j in range(100): # to get very different realisation of input phase screen
            atm.move_atmos()
        
        # Get "full" centroids : open loop WFS centroiding
        wfs.sensors_trace(0,"atmos",atm)
        wfs.sensors_compimg(0)
        wfs.slopes_geom(0,0)
        full_slopes[:,i,1] = wfs.get_slopes(0)
        rtc.docentroids(0)
        full_slopes[:,i,0] = rtc.getCentroids(0)
        
        #Get parallel and orthogonal composantes for LS controller
        rtc.setCom(0,reset_com)
        rtc.docontrol(0)
        #rtc.applycontrol(0,dms)
        dms.set_comm("pzt",0.,rtc.getCom(0))
        dms.shape_dm("pzt",0.)
        wfs.sensors_trace(0,"all",atm,dms)
        wfs.sensors_compimg(0)
        rtc.docentroids(0)
        ortho_slopes[:,i,0] = rtc.getCentroids(0)
        
        wfs.sensors_trace(0,"dm",atm,dms,1)
        wfs.sensors_compimg(0)
        rtc.docentroids(0)
        parallel_slopes[:,i,0] = rtc.getCentroids(0)* -1.0
        
        #Get parallel and orthogonal composantes for geo controller
        tar.atmos_trace(0,atm)
        rtc.docontrol_geo(1,dms,tar,0)
        rtc.applycontrol(1,dms)
        
        wfs.sensors_trace(0,"all",atm,dms)
        wfs.sensors_compimg(0)
        wfs.slopes_geom(0,0)
        ortho_slopes[:,i,1] = wfs.get_slopes(0)
        
        wfs.sensors_trace(0,"dm",atm,dms,1)
        wfs.sensors_compimg(0)
        wfs.slopes_geom(0,0)
        parallel_slopes[:,i,1] = wfs.get_slopes(0)* -1.0
        
        
        
        print("\r Recording... %d%%"%(i*100/n), end=' ')
        
def slopes_analysis(full,ortho,parallel):
    full_ls = full_slopes[:,:,0]
    full_geo = full_slopes[:,:,1]
    ortho_ls = ortho_slopes[:,:,0]
    ortho_geo = ortho_slopes[:,:,1]
    parallel_ls = parallel_slopes[:,:,0]
    parallel_geo = parallel_slopes[:,:,1]
    
    
    err_ls = full_ls - (ortho_ls+parallel_ls)
    err_geo = full_geo - (ortho_geo+parallel_geo)
    
    print("geo error = ",np.mean(np.std(err_geo,axis=1)), " arcsec rms")
    print("ls error = ",np.mean(np.std(err_ls,axis=1))," arcsec rms")
    
def test_loop(niter,gain,amp):
    com = 0
    bk=0
    tur = np.zeros(niter)
    com = np.zeros(niter)
    brcom = np.zeros(niter)
    bk = np.zeros(niter)
    b = np.zeros(niter)
    tur[niter/2:] = amp
    err = np.zeros(niter)
    for i in range(niter):
        #tur[i] = np.random.randn(1)
        mes = tur[i]+com[i-1]  # je fais la mesure sans bruit
        b[i] = 0.01*np.random.randn(1)
        mesb = mes+b[i]  # je fais la mesure avec bruit
        
        # ici, je fais la diff des mesures pour trouver le bruit
        b[i] = mesb - mes
        # recurrence juste sur le bruit
        brcom[i] = brcom[i-1] - gain*b[i]
        # recurrence sur la mes sans bruit
        bk[i] = bk[i-1] - gain * mes
        # calcul de la commande totale
        com[i] = brcom[i] + bk[i]
        
        err[i] = com[i] - bk[i]
        
    pl.ion()
    pl.clf()
    pl.plot(err,color="blue")
    pl.plot(com,color="green")
    pl.plot(bk,color="red")
    #pl.plot(b,color="black")

    return err
    
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
    
    
#############################################################################################################
################################################# MAIN ANALYSIS #############################################
#############################################################################################################
#Allocate buffers
niter = 10000
nslopes = 2*config.p_wfss[0]._nvalid
full_slopes = np.zeros((nslopes,niter,2))
ortho_slopes = np.zeros((nslopes,niter,2))
parallel_slopes = np.zeros((nslopes,niter,2))
#loop(niter,full_slopes,ortho_slopes,parallel_slopes)
#slopes_analysis(full_slopes,ortho_slopes,parallel_slopes)





