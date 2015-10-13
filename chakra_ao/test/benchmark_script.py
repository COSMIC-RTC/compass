import os
import sys
import numpy as np
import chakra as ch
import chakra_ao as ao
import time
import datetime
from subprocess import check_output

CHAKRA_AO=os.environ.get('CHAKRA_AO')
if(CHAKRA_AO is None):
    raise EnvironmentError("Environment variable 'CHAKRA_AO' must be define")

CHAKRA_AO_SAVEPATH=CHAKRA_AO+"/data/"
PARPATH=CHAKRA_AO_SAVEPATH+"par/par4bench/"
BENCH_SAVEPATH=CHAKRA_AO_SAVEPATH+"bench-results/"

c=ch.chakra_context()


def script4bench(filename,centroider,controller, fwrite=True):
    """

    :parameters:
        filename: (str) : parameters filename

        centroider: (str) : centroider type

        controller: (str) : controller type
    """

    timer=ch.chakra_timer()

    #times measured
    synctime=0.
    move_atmos_time=0.
    t_raytrace_atmos_time=0.
    t_raytrace_dm_time=0.
    s_raytrace_atmos_time=0.
    s_raytrace_dm_time=0.
    comp_img_time=0.
    docentroids_time=0.
    docontrol_time=0.
    applycontrol_time=0.

    #reading parfile
    vessel = {}
    execfile(filename, vessel)
    p_centroiders=vessel['p_centroiders']
    p_controllers=vessel['p_controllers']
    p_loop=vessel['p_loop']
    p_wfss=vessel['p_wfss']
    p_atmos=vessel['p_atmos']
    p_tel=vessel['p_tel']
    p_geom=vessel['p_geom']
    p_target=vessel['p_target']
    p_dms=vessel['p_dms']
    p_rtc=vessel['p_rtc']

    p_centroiders[0].set_type(centroider)

    if(centroider=="tcog"): p_centroiders[0].set_thresh(9.)
    elif(centroider=="bpcog"): p_centroiders[0].set_nmax(16)
    elif(centroider=="geom"): p_centroiders[0].set_type("cog")
    elif(centroider=="wcog"):
        p_centroiders[0].set_type_fct("gauss")
        p_centroiders[0].set_width(2.0)
    elif(centroider=="corr"):
        p_centroiders[0].set_type_fct("gauss")
        p_centroiders[0].set_width(2.0)


    p_controllers[0].set_type(controller)
    if(controller=="modopti"):
        p_controllers[0].set_type("ls")
        p_controllers[0].set_modopti(1)
    

    p_loop.set_niter(2000)

    ch.threadSync()
    timer.start()
    ch.threadSync()
    synctime=timer.stop()
    timer.reset()


    #init system
    timer.start()
    wfs=ao.wfs_init(p_wfss,p_atmos,p_tel,p_geom,p_target,p_loop, 1,0,p_dms)
    ch.threadSync()
    wfs_init_time=timer.stop()-synctime
    timer.reset()

    timer.start()
    atm=p_atmos.atmos_init(c,p_tel,p_geom,p_loop)
    ch.threadSync()
    atmos_init_time=timer.stop()-synctime
    timer.reset()

    timer.start()
    dms=ao.dm_init(p_dms,p_wfss[0],p_geom,p_tel)
    ch.threadSync()
    dm_init_time=timer.stop()-synctime
    timer.reset()

    timer.start()
    target=p_target.target_init(c,p_atmos,p_geom,p_tel,p_wfss,wfs,p_dms)
    ch.threadSync()
    target_init_time=timer.stop()-synctime
    timer.reset

    timer.start()
    rtc=ao.rtc_init(wfs,p_wfss,dms,p_dms,p_geom,p_rtc,p_atmos,atm,p_tel,p_loop,p_target)
    ch.threadSync()
    rtc_init_time=timer.stop()-synctime
    timer.reset

    print "... Done with inits !"

    strehllp=[]
    strehlsp=[]

    if(controller=="modopti"):
        for zz in xrange(2048):
            atm.move_atmos()

    for cc in xrange(p_loop.niter):
        ch.threadSync()
        timer.start()
        atm.move_atmos()
        move_atmos_time+=timer.stop()-synctime
        timer.reset()

        if(p_controllers[0].type_control!="geo"):
            if((p_target is not None) and (rtc is not None)):
                for i in xrange(p_target.ntargets):
                    timer.start()
                    target.atmos_trace(i,atm)
                    ch.threadSync()
                    t_raytrace_atmos_time+=timer.stop()-synctime
                    timer.reset()

                    if(dms is not None):
                        timer.start()
                        target.dmtrace(i,dms)
                        ch.threadSync()
                        t_raytrace_dm_time+=timer.stop()-synctime
                        timer.reset()

            if(p_wfss is not None and wfs is not None):
                for i in xrange(len(p_wfss)):
                    timer.start()
                    wfs.sensors_trace(i,"atmos",atm)
                    ch.threadSync()
                    s_raytrace_atmos_time+=timer.stop()-synctime
                    timer.reset()
                    
                    if(not p_wfss[i].openloop and dms is not None):
                        timer.start()
                        wfs.sensors_trace(i,"dm",dms=dms)
                        ch.threadSync()
                        s_raytrace_dm_time+=timer.stop()-synctime
                        timer.reset()

                    timer.start()
                    wfs.sensors_compimg(i)
                    ch.threadSync()
                    comp_img_time+=timer.stop()-synctime
                    timer.reset()

            if(p_rtc is not None and
               rtc is not None and 
               p_wfss is not None and
               wfs is not None):
                if(centroider=="geom"):
                    timer.start()
                    rtc.docentroids_geom(0)
                    ch.threadSync()
                    docentroids_time+=timer.stop()-synctime
                    timer.reset()
                else:
                    timer.start()
                    rtc.docentroids(0)
                    ch.threadSync()
                    docentroids_time+=timer.stop()-synctime
                    timer.reset()

            if(dms is not None):
                timer.start()
                rtc.docontrol(0)
                ch.threadSync()
                docontrol_time+=timer.stop()-synctime
                timer.reset()

                timer.start()
                rtc.applycontrol(0,dms)
                ch.threadSync()
                applycontrol_time+=timer.stop()-synctime
                timer.reset()

        else:
            if(p_target is not None and target is not None):
                for i in xrange(p_target.ntargets):
                    timer.start()
                    target.atmos_trace(i,atm)
                    ch.threadSync()
                    t_raytrace_atmos_time+=timer.stop()-synctime
                    timer.reset()

                    if(dms is not None):
                        timer.start()
                        rtc.docontrol_geo()
                        ch.threadSync()
                        docontrol_time+=timer.stop()-synctime
                        timer.reset()

                        timer.start()
                        rtc.applycontrol(0,dms)
                        ch.threadSync()
                        applycontrol_time+=timer.stop()-synctime
                        timer.reset()

                        timer.start()
                        target.dmtrace(i,dms)
                        ch.threadSync()
                        t_raytrace_dm_time+=timer.stop()-synctime
                        timer.reset()

        strehltmp = target.get_strehl(0)
        strehlsp.append(strehltmp[0])
        if(cc>50):
            strehllp.append(strehltmp[1])


    print "\n done with simulation \n"
    print "\n Final strehl : \n", strehllp[len(strehllp)-1]

    move_atmos_time/=p_loop.niter /1000.
    t_raytrace_atmos_time/=p_loop.niter /1000.
    t_raytrace_dm_time/=p_loop.niter /1000.
    s_raytrace_atmos_time/=p_loop.niter /1000.
    s_raytrace_dm_time/=p_loop.niter /1000.
    comp_img_time/=p_loop.niter /1000.
    docentroids_time/=p_loop.niter /1000.
    docontrol_time/=p_loop.niter /1000.
    applycontrol_time/=p_loop.niter /1000.

    time_per_iter=move_atmos_time+ t_raytrace_atmos_time+\
                  t_raytrace_dm_time+ s_raytrace_atmos_time+\
                  s_raytrace_dm_time+ comp_img_time+\
                  docentroids_time+ docontrol_time+\
                  applycontrol_time


    if(p_wfss[0].gsalt>0):
        type="lgs "
    else:
        type="ngs "

    if(p_wfss[0]>3):
        type+="noisy "

    type+=p_wfss[0].type_wfs



    date=datetime.datetime.now().strftime("%a %d. %b %Y %H:%M")
    svnversion=check_output("svnversion").replace("\n","")
    hostname=check_output("hostname").replace("\n","")

    savefile=BENCH_SAVEPATH+"result_"+hostname+"_scao_revision_"+svnversion+".csv"
    SRfile=BENCH_SAVEPATH + "SR_"+hostname+"_scao_revision_"+svnversion+".csv"

    print date, svnversion, hostname, savefile, SRfile

    if(fwrite):
        print "writing files"
        if(not os.path.isfile(savefile)):
            f=open(savefile,"w")
            f.write("--------------------------------------------------------------------------\n")
            f.write("Date :"+date+"\t Revision : "+svnversion+"\n")
            f.write("--------------------------------------------------------------------------\n")
            f.write("System type\tnxsub\twfs.npix\tNphotons\tController\tCentroider\tFinal SR LE\tAvg. SR SE\trms SR SE\twfs_init\tatmos_init\tdm_init\ttarget_init\trtc_init\tmove_atmos\tt_raytrace_atmos\tt_raytrace_dm\ts_raytrace_atmos\ts_raytrace_dm\tcomp_img\tdocentroids\tdocontrol\tapplycontrol\titer_time\tAvg.gain")

        else:
            f=open(savefile,'a')


        if(controller=="modopti"):
            G=rtc.get_mgain()
        else:
            G=0.

        f.write("\n"+ type +"\t"+ str(p_wfss[0].nxsub) +"\t"+ str(p_wfss[0].npix) +"\t"+
                str(p_wfss[0]._nphotons) +"\t"+ controller +"\t"+ centroider +"\t"+
                str(strehllp[len(strehllp)-1]) +"\t"+ str(np.mean(strehlsp)) +"\t"+
                str(np.std(strehllp)) +"\t"+ str(wfs_init_time) +"\t"+
                str(atmos_init_time) +"\t"+ str(dm_init_time) +"\t"+
                str(target_init_time) +"\t"+ str(rtc_init_time) +"\t"+
                str(move_atmos_time) +"\t"+ str(t_raytrace_atmos_time) +"\t"+
                str(t_raytrace_dm_time) +"\t"+ str(s_raytrace_atmos_time) +"\t"+
                str(s_raytrace_dm_time) +"\t"+ str(comp_img_time) +"\t"+
                str(docentroids_time) +"\t"+ str(docontrol_time) +"\t"+
                str(applycontrol_time) +"\t"+str(time_per_iter) +"\t"+str(G))
        
        f.close()



        if(not os.path.isfile(SRfile)):
            f=open(SRfile,"w")
            f.write("--------------------------------------------------------------------------\n")
            f.write("Date :"+date+"\t Revision : "+svnversion+"\n")
            f.write("--------------------------------------------------------------------------\n")
            f.write("System type\tnxsub\twfs.npix\tNphotons\tController\tCentroider")

        else:
            f=open(SRfile,"a")

        f.write("\n"+ type +"\t"+ str(p_wfss[0].nxsub) +"\t"+ str(p_wfss[0].npix) +"\t"+
                str(p_wfss[0]._nphotons) +"\t"+ controller +"\t"+ centroider+"\t"+str(strehllp[len(strehllp)-1]))

        f.close()


#DONE function script4bench

if(len(sys.argv)!=4 and len(sys.argv)!=5):
    error="wrong number of argument. Got "+len(sys.argv)+" (expect 4)\ncommande line should be: 'python benchmark_script.py <filename> <centroider> <controller>"
    raise StandardError(error)

filename=PARPATH+sys.argv[1]
centroider=sys.argv[2]
controller=sys.argv[3]
if (len(sys.argv)==5):
    fwrite=int(sys.argv[4])
    print "fwrite=",fwrite
    script4bench(filename,centroider,controller,fwrite)
else:
    script4bench(filename,centroider,controller)

