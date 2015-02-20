yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

require,yoga_ao_top+"/yorick/yoga_ao.i";
//require,yoga_ao_top+"/yorick/yaokl.i";
//require,yoga_ao_top+"/ywidgets/widget_wfs.i";

#include"fits-utils.i"

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data/";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

YOGA_AO_PARPATH = YOGA_AO_SAVEPATH+"par/par4bench/";
BENCH_SAVE_PATH = YOGA_AO_SAVEPATH+"bench-results/";

mkdirp,YOGA_AO_PARPATH;
mkdirp,BENCH_SAVE_PATH;

/******************************* SCRIPT ***************************************************/

func script4bench(filename,centroider,controller){

  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;
  extern strehlsp,strehllp;

  //Profil variables
  move_atmos_time = t_raytrace_atmos_time = t_raytrace_dm_time = s_raytrace_atmos_time = s_raytrace_dm_time = comp_img_time = docentroids_time = docontrol_time = applycontrol_time = 0.;

  // reading parfile
  read_parfile,filename;
  y_centroiders(1).type = centroider;
  if(y_centroiders(1).type == "tcog") y_centroiders(1).thresh = 9.;
  if(y_centroiders(1).type == "bpcog") y_centroiders(1).nmax = 16;
  if(y_centroiders(1).type == "geom") y_centroiders(1).type = "cog";
  if(y_centroiders(1).type == "wcog") {
    y_centroiders(1).type_fct = "gauss";
    y_centroiders(1).width   = 2.0;
  }
  if(y_centroiders(1).type == "corr") {
    y_centroiders(1).type_fct = "gauss";
    y_centroiders(1).width   = 2.0;
  }
  
  y_controllers(1).type = controller;
  if(y_controllers(1).type == "modopti"){
    y_controllers(1).type = "ls";
    y_controllers(1).modopti = 1;
  }

  y_loop.niter = 2000;
 _yogaThreadSync;
  tic;
  _yogaThreadSync;
  synctime = tac();

  // init system
  tic;
  wfs_init;
_yogaThreadSync;
  wfs_init_time = tac() - synctime;

  tic;
  atmos_init;
  _yogaThreadSync;
atmos_init_time = tac()- synctime;

  tic;
  dm_init;
  _yogaThreadSync;
  dm_init_time = tac()- synctime;

  tic;
  target_init;
  _yogaThreadSync;
  target_init_time = tac()- synctime;

  tic;
  rtc_init,clean=1;
   _yogaThreadSync;
 rtc_init_time = tac()- synctime;

  write,"... Done with inits !";

  strehllp = strehlsp = [];

  for (cc=1;cc<=y_loop.niter;cc++) {
    _yogaThreadSync;
    tic,0;
    move_atmos,g_atmos;
    _yogaThreadSync;
    move_atmos_time += tac(0)- synctime;

    if(y_controllers(1).type != "geo"){
      if ((y_target != []) && (g_target != [])) {
	// loop on targets
	for (i=1;i<=y_target.ntargets;i++) {
	  tic,1;
	  target_atmostrace,g_target,i-1,g_atmos;
	    _yogaThreadSync;
	    t_raytrace_atmos_time += tac(1)- synctime;
	  if (g_dm != []) {
	    tic,2;
	    target_dmtrace,g_target,i-1,g_dm;
	      _yogaThreadSync;
  t_raytrace_dm_time += tac(2)- synctime;
	  }
	}
      }
    
      if ((y_wfs != []) && (g_wfs != [])) {
	// loop on wfs
	for (i=1;i<=numberof(y_wfs);i++) {
	  tic,3;
	  sensors_trace,g_wfs,i-1,"atmos",g_atmos;
    _yogaThreadSync;
	  s_raytrace_atmos_time += tac(3)- synctime;
	  if ((!y_wfs(i).openloop) && (g_dm != [])) {
	    tic,4;
	    sensors_trace,g_wfs,i-1,"dm",g_dm,0;
	      _yogaThreadSync;
  s_raytrace_dm_time += tac(4)- synctime;
	  }

	  if(y_wfs(i).type=="cog") {
	    tic,5;
	    sensors_compimg_tele,g_wfs,i-1;
	     _yogaThreadSync;
   comp_img_time += tac(5)- synctime;
	  } else {
	    tic,5;
	    sensors_compimg,g_wfs,i-1;
	      _yogaThreadSync;
  comp_img_time += tac(5)- synctime;
	  }
	}   
      }

      if ((y_rtc != []) && (g_rtc != [])
	  && (y_wfs != []) && (g_wfs != [])) {
	if(centroider == "geom"){
	  tic,6;
	  rtc_docentroids_geom,g_rtc,0;
	    _yogaThreadSync;
  docentroids_time += tac(6)- synctime;
	}
	else {
	  tic,6;
	  rtc_docentroids,g_rtc,0;
	     _yogaThreadSync;
 docentroids_time += tac(6)- synctime;
	}
	// compute command and apply
	if (g_dm != []) {
	  tic,7;
	  rtc_docontrol,g_rtc,0;
	    _yogaThreadSync;
  docontrol_time += tac(7)- synctime;
	  tic,8;
  rtc_applycontrol,g_rtc,0,g_dm;
	    _yogaThreadSync;
  applycontrol_time += tac(8)- synctime;
	}
      }
    }
    else{
      if ((y_target != []) && (g_target != [])) {
	for (i=1;i<=y_target.ntargets;i++) {
	  tic;
	  target_atmostrace,g_target,i-1,g_atmos;
	  _yogaThreadSync;
	  t_raytrace_atmos_time += tac();
	  if (g_dm != []) {
	    tic;
	    rtc_docontrol_geo,g_rtc,0,g_dm,g_target,0;
	    _yogaThreadSync;
	    docontrol_time += tac();
	    tic;
	    rtc_applycontrol,g_rtc,0,g_dm;
	    _yogaThreadSync;
	    applycontrol_time += tac();
	    tic;
	    target_dmtrace,g_target,i-1,g_dm;
	    _yogaThreadSync;
	    t_raytrace_dm_time += tac();
	  }
	}
      }
    }

    strehltmp = target_getstrehl(g_target,0);
    grow,strehlsp,strehltmp(1);
    if(cc>50)
      grow,strehllp,strehltmp(2);
    
  }
  
  write,"\n done with simulation \n";
  write,"\n Final strehl : \n",strehllp(0);
  
  move_atmos_time /= y_loop.niter / 1000.; // millisecond
  t_raytrace_atmos_time /= y_loop.niter / 1000.;
  t_raytrace_dm_time /= y_loop.niter / 1000.;
  s_raytrace_atmos_time  /= y_loop.niter / 1000.;
  s_raytrace_dm_time  /= y_loop.niter / 1000.;
  comp_img_time /= y_loop.niter / 1000.;
  docentroids_time  /= y_loop.niter / 1000.;
  docontrol_time /= y_loop.niter / 1000.; 
  applycontrol_time /= y_loop.niter / 1000.;
  time_per_iter = move_atmos_time + t_raytrace_atmos_time + t_raytrace_dm_time + s_raytrace_atmos_time + s_raytrace_dm_time + comp_img_time + docentroids_time + docontrol_time + applycontrol_time;

  if(y_wfs(1).gsalt > 0)
    type = "lgs";
  else type = "ngs";
  if(y_wfs(1).gsmag > 3)
    type += " noisy";
  
  tmp = timestamp();
  date = strpart(tmp,5:11)+strpart(tmp,21:);
  svnversion = rdfile(popen("svnversion",0))(1);
  savefile = BENCH_SAVE_PATH + "results_scao_revision_"+svnversion+".csv";
  SRfile = BENCH_SAVE_PATH + "SR_scao_revision_"+svnversion+".csv";

  if(!(fileExist(savefile))){
    f=open(savefile,"w");
    write,f,"--------------------------------------------------------------------------";
    write,f,format="Date : %s\t Revision : %s \t  %s \n",date,svnversion,context_get_device_name(context_getactivedevice());
    write,f,"--------------------------------------------------------------------------";
    write,f,"System type\tnxsub\twfs.npix\tNphotons\tController\tCentroider\tFinal SR LE\tAvg. SR SE\trms SR SE\twfs_init\tatmos_init\tdm_init\ttarget_init\trtc_init\tmove_atmos\tt_raytrace_atmos\tt_raytrace_dm\ts_raytrace_atmos\ts_raytrace_dm\tcomp_img\tdocentroids\tdocontrol\tapplycontrol\titer_time";
  }
  else
    f=open(savefile,"a");

  write,f,format="%s\t%d\t%d\t%f\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",type,y_wfs(1).nxsub,y_wfs(1).npix,y_wfs(1)._nphotons,controller,centroider,strehllp(0),avg(strehlsp),strehlsp(rms),wfs_init_time,atmos_init_time,dm_init_time,target_init_time,rtc_init_time,move_atmos_time,t_raytrace_atmos_time,t_raytrace_dm_time,s_raytrace_atmos_time,s_raytrace_dm_time,comp_img_time,docentroids_time,docontrol_time,applycontrol_time,time_per_iter;

  close,f;

  if(!(fileExist(SRfile))){
    f2=open(SRfile,"w");
    write,f2,"--------------------------------------------------------------------------";
    write,f2,format="Date : %s\t Revision : %s \t %s \n",date,svnversion,context_get_device_name(context_getactivedevice());
    write,f2,"--------------------------------------------------------------------------";
    write,f2,"System type\tnxsub\twfs.npix\tNphotons\tController\tCentroider";
  }
  else
    f2=open(SRfile,"a");

  write,f2,format="%s\t%d\t%d\t%f\t%s\t%s",type,y_wfs(1).nxsub,y_wfs(1).npix,y_wfs(1)._nphotons,controller,centroider;

  for(i=1 ; i<=y_loop.niter ; i++){
    write,f2,format="\t %f ",strehlsp(i);
  }
  write,f2,"";

  close,f2;

}

/******************************************* BATCH ******************************************/
if(batch()) {
  parameters = get_argv();
  filename = YOGA_AO_PARPATH + parameters(2);
  centroider = parameters(3);
  controller = parameters(4);
  script4bench,filename,centroider,controller;
} 


