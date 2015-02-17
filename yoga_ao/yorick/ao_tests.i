yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

require,yoga_ao_top+"/yorick/yoga_ao.i";
//require,yoga_ao_top+"/ywidgets/widget_wfs.i";

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

YOGA_AO_PARPATH = YOGA_AO_SAVEPATH+"/par/par4scripts/";
mkdirp,YOGA_AO_PARPATH;


/*
    _   _                       
   / \ | |_ _ __ ___   ___  ___ 
  / _ \| __| '_ ` _ \ / _ \/ __|
 / ___ \ |_| | | | | | (_) \__ \
/_/   \_\__|_| |_| |_|\___/|___/
                                

 */

func script_atmos_full(void)
{
  extern y_geom,y_tel,y_loop,y_atmos;
  extern g_atmos,g_target;
  extern ipupil;

  tabres = array(0.0f,4,4);

  filenames = [["atmos_1layer.par","atmos_1layer_vlt.par","atmos_1layer_20m.par",
                "atmos_1layer_elt.par"],
               ["atmos_1layer_5targets.par","atmos_1layer_vlt_5targets.par",
                "atmos_1layer_20m_5targets.par","atmos_1layer_elt_5targets.par"],
               ["atmos_4layers.par","atmos_4layers_vlt.par","atmos_4layers_20m.par",
                "atmos_4layers_elt.par"],
               ["atmos_12layers.par","atmos_12layers_vlt.par","atmos_12layers_20m.par",
                "atmos_12layers_elt.par"]];
                 
  for (i=1;i<=4;i++) {
    for (j=1;j<=4;j++) {
      read_parfile,YOGA_AO_PARPATH+filenames(j,i);

      if (y_loop.niter == []) y_loop.niter = 1000;
      
      //here we just simulate atmosphere so pupdiam is fixed
      // arbitrarily.
      pupdiam = long(2*y_tel.diam*20.); // we assume subaps of 50cm and 20 phase pix per subaps

      write,"Doing atmos inits on the GPU";
  
      geom_init,pupdiam;
      
      atmos_init;
      
      write,"... Done !";
      
      write,"Creating targets on the GPU";
      
      target_init;
      
      write,"... Done !";
      
      write,format="Starting loop on : %d iterations\n",y_loop.niter;
      
      mytime = tic();
      for (cc=1;cc<=y_loop.niter;cc++) {
        //tinter = tic();
        move_atmos,g_atmos;
        /*
        if (cc % 10 == 0) {
          time_move = tac(mytime)/cc;
          //write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
        }
        */
      }
      write,"... Done !";
      // first dim is the telescope size
      // second dim is the type of parfile
      tabres(j,i) = tac(mytime)/y_loop.niter;
      
      //write,format="Average atmos move gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
    }
  }
  return tabres;
}



func script_atmos(filename)
{
  extern y_geom,y_tel,y_loop,y_atmos;
  extern g_atmos,g_target;
  extern ipupil;

  if (filename == []) filename = YOGA_AO_PARPATH+"atmos_1layer.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"atmos_12layers.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 1000;
  
  //here we just simulate atmosphere so pupdiam is fixed
  // arbitrarily.
  pupdiam = long(32*y_tel.diam/4.0);

  write,"Doing atmos inits on the GPU";

  geom_init,pupdiam;

  atmos_init;
  
  write,"... Done !";
  
  write,"Creating targets on the GPU";

  target_init;

  write,"... Done !";

  write,format="Starting loop on : %d iterations\n",y_loop.niter;
  
  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    //tinter = tic();
    move_atmos,g_atmos;
    if (cc % 10 == 0) {
      time_move = tac(mytime)/cc;
      write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
    }
  }
  write,"... Done !";
  write,format="Average atmos move gpu time : %.4f s\n",tac(mytime)/y_loop.niter;

  write,"Doing atmos inits on the CPU";

  pupixsize = y_tel.diam / y_geom.pupdiam;
  
  cpu_r0 = y_atmos.r0;

  screen_cpu = create_screen(cpu_r0,pupixsize,y_geom.pupdiam,A,B,ist,100.);

  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    //tinter = tic();
    screen_cpu = extrude(screen_cpu,float(cpu_r0 /pupixsize), A, B, ist);
    if (cc % 10 == 0) {
      time_move = tac(mytime)/cc;
      write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
    }
  }
  write,"... Done !";
  write,format="Average extrude cpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  
  // now that we have both gpu and cpu inits, compare screen
  // need to compute phase structure function
  dphi_cpu = dphi_gpu1 = dphi_gpu2 = 0.0f;
  
  write,"Comparing CPU vs GPU screens";

  niter = 100000;
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(cpu_r0 /pupixsize), A, B, ist);
     dphi_cpu   += calc_dphis(screen_cpu(1:y_geom.pupdiam,1:y_geom.pupdiam));
     move_atmos,g_atmos;
     target_atmostrace,g_target,0,g_atmos;
     screen_gpu2 = target_getphase(g_target,0);
     dphi_gpu2   += calc_dphis(screen_gpu2);
     if (cc % 10 == 0) {
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
     }
  }
  
  write,"... Done !";

  fma;
  circa = circavg_quad(dphi_gpu2(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter;
  plg,6.88*(indgen(numberof(circa))/float(cpu_r0 /pupixsize))^(5./3.),marks=0,width=4;
  plg,circavg_quad(dphi_cpu(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter,color="red",marks=0,width=4;
  plg,circa,color="green",marks=0,width=4;
  pltitle,"<(phase-phase(1,1))^2^ > , GPU (g) vs CPU (r)"
}

func script_atmosgpu(filename)
{
  extern y_geom,y_tel,y_loop,y_atmos;
  extern g_atmos,g_target;
  extern ipupil;

  if (filename == []) filename = YOGA_AO_PARPATH+"atmos_1layer.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"atmos_12layers.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 1000;
  
  //here we just simulate atmosphere so pupdiam is fixed
  // arbitrarily.
  pupdiam = long(32*y_tel.diam/4.0);

  write,"Doing atmos inits on the GPU";

  geom_init,pupdiam;

  atmos_init;
  
  write,"... Done !";
  
  write,"Creating targets on the GPU";

  target_init;

  write,"... Done !";

  write,format="Starting loop on : %d iterations\n",y_loop.niter;
  
  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    //tinter = tic();
    move_atmos,g_atmos;
    if (cc % 10 == 0) {
      time_move = tac(mytime)/cc;
      write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
    }
  }
  write,"... Done !";
  write,format="Average atmos move gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  
  
  // now that we have both gpu and cpu inits, compare screen
  // need to compute phase structure function
  dphi_gpu = 0.0f;
  
  write,"Comparing CPU vs GPU screens";

  niter = 100000;
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
     move_atmos,g_atmos;
     target_atmostrace,g_target,0,g_atmos;
     screen_gpu = target_getphase(g_target,0);
     dphi_gpu   += calc_dphis(screen_gpu);
     if (cc % 10 == 0) {
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
     }
  }
  
  write,"... Done !";

  fma;
  pupixsize = y_tel.diam / y_geom.pupdiam;
  circa = circavg_quad(dphi_gpu(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter;
  plg,6.88*(indgen(numberof(circa))/float(y_atmos.r0 /pupixsize))^(5./3.),marks=0,width=4;
  plg,circa,color="green",marks=0,width=4;
  pltitle,"<(phase-phase(1,1))^2^ > , GPU (g) vs CPU (r)"
}

func script_atmoscpu(void)
{
  niter = 100000;
  pupdiam = 64;
  
  dphi_cpu0 = dphi_cpu1 = dphi_cpu2 = dphi_cpu3 = 0.0f;

  screen_cpu = create_screen(0.16,0.0313725,pupdiam,A,B,ist,100.);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
      dphi_cpu0   += calc_dphis(screen_cpu);
    time_move = tac(mytime)/cc;
    write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam+16,A,B,ist,100.);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
    dphi_cpu1   += calc_dphis(screen_cpu(1:pupdiam,1:pupdiam));
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam*2,A,B,ist,100.);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
    dphi_cpu2   += calc_dphis(screen_cpu(1:pupdiam,1:pupdiam));
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }

  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam*4,A,B,ist,100.);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
      dphi_cpu3   += calc_dphis(screen_cpu(1:pupdiam,1:pupdiam));
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  plg,6.88*(indgen(pupdiam)/5.1)^(5./3.),marks=0,width=4;
  plg,circavg_quad(dphi_cpu0(1:pupdiam,1:pupdiam))/niter,color="red",marks=0,width=4;
  plg,circavg_quad(dphi_cpu1(1:pupdiam,1:pupdiam))/niter,color="blue",marks=0,width=4;
  plg,circavg_quad(dphi_cpu2(1:pupdiam,1:pupdiam))/niter,color="green",marks=0,width=4;
  plg,circavg_quad(dphi_cpu3(1:pupdiam,1:pupdiam))/niter,color="yellow",marks=0,width=4;
}

/*
__        _______ ____  
\ \      / /  ___/ ___| 
 \ \ /\ / /| |_  \___ \ 
  \ V  V / |  _|  ___) |
   \_/\_/  |_|   |____/ 
                        

*/

func script_wfs(filename,&slpgeom1,&slpgeom2,verbose=)
{

  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc_dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 1000;

  // init system
  wfs_init;

  atmos_init;
 
  target_init;

  rtc_init;

  if (verbose) write,"... Done with inits !";
  write,"The following objects have been initialized on the GPU :";
  write,"--------------------------------------------------------";
  g_atmos;
  write,"--------------------------------------------------------";
  g_wfs;
  write,"--------------------------------------------------------";
  g_target;
  write,"--------------------------------------------------------";
  g_rtc;

  // warming up !
  // move sky
  move_atmos,g_atmos;
  
  // build wfs image
  for (i=1;i<=numberof(y_wfs);i++) {

    sensors_trace,g_wfs,i-1,"atmos",g_atmos;

    sensors_compimg,g_wfs,i-1;
  }

  // init slopes arrays
  tmp = [];
  for (i=1;i<=numberof(y_wfs);i++)
    grow,tmp,sensors_getslopes(g_wfs,i-1);
  slpgeom1 = slpgeom2 = array(0.0,numberof(tmp),y_loop.niter);
  
  centroiders = *y_rtc.centroiders;
  ncentro = numberof(centroiders);
  if (ncentro > 0) {
    tmp = [];
    for (i=1;i<=ncentro;i++)
      grow,tmp,sensors_getslopes(g_wfs,centroiders(i).nwfs-1);
    cb_slopes = array(0.0,numberof(tmp),y_loop.niter);
  }
  
  if (verbose) write,format="Starting loop on : %d iterations\n",y_loop.niter;

  mytime = tic();
  
  for (cc=1;cc<=y_loop.niter;cc++) {

  // move sky
    move_atmos,g_atmos;

    if ((y_wfs != []) && (g_wfs != [])) {
      tmp1 = tmp2 = [];
  // build wfs image & get geom slopes
      for (i=1;i<=numberof(y_wfs);i++) {

        sensors_trace,g_wfs,i-1,"atmos",g_atmos;

        sensors_compimg,g_wfs,i-1;

        slopes_geom,g_wfs,i-1,0;

        grow,tmp1,sensors_getslopes(g_wfs,i-1);

        slopes_geom,g_wfs,i-1,1;

        grow,tmp2,sensors_getslopes(g_wfs,i-1);
      }
      slpgeom1(,cc) = tmp1;
      slpgeom2(,cc) = tmp2;
      tmp = [];
  // compute slopes
      if (ncentro > 0) {
        for (i=1;i<=ncentro;i++) {
          if (centroiders(i).type == "cog") {
          sensors_compslopes,g_wfs,centroiders(i).nwfs-1,g_rtc,i-1;
          tmps = sensors_getslopes(g_wfs,centroiders(i).nwfs-1);
          tmps = (tmps - (y_wfs(centroiders(i).nwfs).npix/2.+0.5))*y_wfs(centroiders(i).nwfs).pixsize;
          }
          if (centroiders(i).type == "tcog") {
            extra = centroiders(i).thresh;
            sensors_compslopes,g_wfs,centroiders(i).nwfs-1,g_rtc,i-1,extra;
          tmps = sensors_getslopes(g_wfs,centroiders(i).nwfs-1);
          tmps = (tmps - (y_wfs(centroiders(i).nwfs).npix/2.+0.5))*y_wfs(centroiders(i).nwfs).pixsize;
          }
          if (centroiders(i).type == "bpcog") {
            extra = long(centroiders(i).nmax);
            sensors_compslopes,g_wfs,centroiders(i).nwfs-1,g_rtc,i-1,extra;
          tmps = sensors_getslopes(g_wfs,centroiders(i).nwfs-1);
          tmps = (tmps - (y_wfs(centroiders(i).nwfs).npix/2.+0.5))*y_wfs(centroiders(i).nwfs).pixsize;
          }
          if (centroiders(i).type == "wcog") {
          sensors_compslopes,g_wfs,centroiders(i).nwfs-1,g_rtc,i-1,extra;
          tmps = sensors_getslopes(g_wfs,centroiders(i).nwfs-1);
          tmps = (tmps - (y_wfs(centroiders(i).nwfs).npix/2.+0.5))*y_wfs(centroiders(i).nwfs).pixsize;
          }
          if (centroiders(i).type == "corr") {
          sensors_compslopes,g_wfs,centroiders(i).nwfs-1,g_rtc,i-1;
          tmps = sensors_getslopes(g_wfs,centroiders(i).nwfs-1);
          tmps = (tmps - (y_wfs(centroiders(i).nwfs).npix/2.+1.))*y_wfs(centroiders(i).nwfs).pixsize;
          }

          grow,tmp,tmps;
        }
        cb_slopes(,cc)=tmp;
      }
    }
    if (verbose) { 
      if (cc % 10 == 0) {
        time_move = tac(mytime)/cc;
        write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
      }
    }
  }

  if (verbose) {
    write,"... done !";
    write,format="average wfs gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  }

  return cb_slopes;
}

func reduce_script(slopes,slopesgeom,nwfs,disp=)
{
  Dssp = y_tel.diam/y_wfs(nwfs).nxsub;
  
  //in arcsec
  mesSlopes  = (slopes-(y_wfs(nwfs).npix/2.+0.5))*y_wfs(nwfs).pixsize;
//mesSlopes  = cb_slopes1;
  trueSlopes = slopesgeom;
  
  // Centroid gain : linear regression of measured data versus real ones
  tab_regressX = tab_regressY = array(0.,y_wfs(nwfs)._nvalid);
  for(j=1;j<=y_wfs(nwfs)._nvalid;j++) {
    tab_regressX(j) = regressAx(mesSlopes(j,), trueSlopes(j,)) + 1e-20;
    tab_regressY(j) = regressAx(mesSlopes(j+y_wfs(nwfs)._nvalid,), trueSlopes(j+y_wfs(nwfs)._nvalid,)) + 1e-20;
  }
  // Centroid error : this assumes centroid gain=1
  Dg = mesSlopes - trueSlopes;
  
  // computation of uncorrelated noise : the true noise (as seen by a closed loop).
  uDg = mesSlopes/_(tab_regressX,tab_regressY)(,-:1:y_loop.niter) - trueSlopes;
  // variance non centree, en arcsec
  noise_level_asec = (uDg^2)(,avg);                               
  // variance traduite en phase (rd^2)
  noiseLevel      = noise_level_asec * (2.*pi*Dssp/(y_wfs(1).lambda *1.e-6))^2 / (RASC*RASC);
  noiseLevelNmRms = sqrt(noiseLevel) * (y_wfs(nwfs).lambda*1e3/2/pi);

  "Slopes noise level (nm rms):";
  noiseLevelNmRms;

  if (disp) {
    window,0;fma;
    tmp=wfs_map(noiseLevelNmRms(1:y_wfs(nwfs)._nvalid),y_wfs(nwfs),type="subaps");
    pli,tmp;
    colorbar,min(tmp),max(tmp);
    pltitle,"x slopes noise level (nm rms)";
    
    window,1;fma;
    tmp=wfs_map(noiseLevelNmRms(y_wfs(nwfs)._nvalid+1:),y_wfs(nwfs),type="subaps");
    pli,tmp;
    colorbar,min(tmp),max(tmp);
    pltitle,"y slopes noise level (nm rms)";
    
    window,2;fma;
    tmp=wfs_map(tab_regressX,y_wfs(nwfs),type="subaps");
    pli,tmp;
    colorbar,min(tmp),max(tmp);
    pltitle,"x slopes gain";
    
    window,3;fma;
    tmp=wfs_map(tab_regressY,y_wfs(nwfs),type="subaps");
    pli,tmp;
    colorbar,min(tmp),max(tmp);
    pltitle,"y slopes gain";
  }
}

func scatter_slopes(slopes,slopesgeom,nwfs,nsub,ref=)
{
  if (ref == []) ref = 0;
  if (ref)
    mesSlopes  = (slopes-(y_wfs(nwfs).npix/2.+0.5))*y_wfs(nwfs).pixsize;
  else
    mesSlopes = slopes;
  trueSlopes = slopesgeom;
  
  // Centroid gain : linear regression of measured data versus real ones
  tab_regressX = tab_regressY = array(0.,y_wfs(nwfs)._nvalid);
  for(j=1;j<=y_wfs(nwfs)._nvalid;j++) {
    tab_regressX(j) = regressAx(mesSlopes(j,), trueSlopes(j,)) + 1e-20;
    tab_regressY(j) = regressAx(mesSlopes(j+y_wfs(nwfs)._nvalid,), trueSlopes(j+y_wfs(nwfs)._nvalid,)) + 1e-20;
  }
  
  window,0; fma;limits;
  plmk, mesSlopes(nsub,), trueSlopes(nsub,), color="red",msize=0.1;
  plmk, mesSlopes(nsub+y_wfs(nwfs)._nvalid,), trueSlopes(nsub+y_wfs(nwfs)._nvalid,), msize=0.1, color="green";
  xytitles,"True","Measured";
  a=limits();
  lim = minmax(a(1:4));
  plg,[lim(1),lim(2)],[lim(1),lim(2)],color="black",marks=0;
  //plg,[min(GG),max(GG)]*tab_regress(Jdisp),[min(GG),max(GG)],color="red",marks=0;
  plt,"y=x",0.22,0.83,tosys=0;
  gridxy,1,1;

}



func check_centroiding(filename,thresh=,nmax=)
{
  if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 1000;

  // init system
  wfs_init;

  atmos_init;
 
  target_init;

  rtc_init;

  write,"... Done with inits !";
  write,"The following objects have been initialized on the GPU :";
  write,"--------------------------------------------------------";
  g_atmos;
  write,"--------------------------------------------------------";
  g_wfs;
  write,"--------------------------------------------------------";
  g_target;
  write,"--------------------------------------------------------";
  g_rtc;

  // move sky
  move_atmos,g_atmos;
  
  // build wfs image
  for (i=1;i<=numberof(y_wfs);i++) {

    sensors_trace,g_wfs,i-1,"atmos",g_atmos;

    sensors_compimg,g_wfs,i-1;
  }

  // check geom slopes
  // gpu version
  slopes_geom,g_wfs,0,0;
  slp=sensors_getslopes(g_wfs,0);
  // yorick version
  res=reform(sensors_getdata(g_wfs,0,"phase")(*y_wfs(1)._phasemap+1),y_wfs(1)._pdiam,y_wfs(1)._pdiam,y_wfs(1)._nvalid);
  slpx = float(res(0,avg,)-res(1,avg,));
  slpy = float(res(avg,0,)-res(avg,1,));
  geom=_(slpx,slpy);
  // in arcsec :
  alpha = 0.206265 / (y_tel.diam/y_wfs(1).nxsub);
  geom *= float(alpha);
  "abs(geom - slopes): [min, max, avg]";
  tmp_array=abs(geom-slp);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  //plmk,geom,slp;
  /*
  // check geom slopes (including pupil model : real derivative)
  slopes_geom,g_wfs,0,0;
  slp=sensors_getslopes(g_wfs,0);
  slopes_geom,g_wfs,0,1;
  slp1=sensors_getslopes(g_wfs,0);
  "geom slopes 0 vs 1: slp-slp1";
  //slp/slp1;
  tmp_array=abs(slp-slp1);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  //plmk,slp1,slp,color="red";
  */
  // first put back slopes to geom to be sure d_slopes are modified in yoga_wfs
  slopes_geom,g_wfs,0,0;
  // check slopes
  sensors_compslopes,g_wfs,0,g_rtc,0;
  slp2=sensors_getslopes(g_wfs,0);
  res2=sensors_getdata(g_wfs,0,"bincube");
  tmp=indices(y_wfs(1).npix);
  slpx=(res2*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/(res2(*,)(sum,));
  slpy=(res2*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/(res2(*,)(sum,));
  //centro=_(slpx,slpy);
  centro=_((slpx-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize,(slpy-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize);
  "cog: centro-slp2";
  //centro/slp2;
  tmp_array=abs(centro-slp2);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  // to convert in arcsec : (slpx-(y_wfs(nwfs).npix/2.+0.5))*y_wfs(nwfs).pixsize

  slopes_geom,g_wfs,0,0;
  // check thresholded slopes 
  thresh = (*y_rtc.centroiders)(2).thresh;
  if (thresh == []) thresh=1000.;
  sensors_compslopes,g_wfs,0,g_rtc,1,thresh;
  slp2=sensors_getslopes(g_wfs,0);
  res2=sensors_getdata(g_wfs,0,"bincube");
  tmp=indices(y_wfs(1).npix);
  for (i=1;i<=dimsof(res2)(4);i++) {
    nn = where(res2(*,i) < thresh);
    if (numberof(nn) > 0) res2(*,i)(nn) *= 0.;
  }
  slpx=(res2*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  slpy=(res2*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  centro=_((slpx-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize,(slpy-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize);
  "tcog: centro-slp2";
  //centro/slp2;
  tmp_array=abs(centro-slp2);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  // to convert in arcsec : (slp2-(y_wfs(nwfs).npix/2.+0.5))*y_wfs(nwfs).pixsize

  // check slopes wcog
  sensors_compslopes,g_wfs,0,g_rtc,3;
  slp2=sensors_getslopes(g_wfs,0);
  res2=sensors_getdata(g_wfs,0,"bincube");
  tmp=indices(y_wfs(1).npix);
  www=*(*y_rtc.centroiders)(4).weights;
  slpx=((res2*www)*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/(res2*www)(*,)(sum,);
  slpy=((res2*www)*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/(res2*www)(*,)(sum,);
  centro=_((slpx-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize,(slpy-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize);
  "wcog";
  //centro/slp2;
  tmp_array=abs(centro-slp2);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  // to convert in arcsec : (slp2-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize
 
  slopes_geom,g_wfs,0,0;
  // check slopes nmax
  nmax = (*y_rtc.centroiders)(3).nmax;
  if (nmax == []) nmax=5;
  sensors_compslopes,g_wfs,0,g_rtc,2,long(nmax);
  slp2=sensors_getslopes(g_wfs,0);
  res2=sensors_getdata(g_wfs,0,"bincube");
  tmp=indices(y_wfs(1).npix);
  for (i=1;i<=dimsof(res2)(4);i++){
    res2(*,i)((sort(res2(*,i)))(::-1)(nmax+1:))=0.;
    res2(*,i)((sort(res2(*,i)))(::-1)(:nmax)) -= res2(*,i)((sort(res2(*,i)))(::-1)(nmax));
  }
  slpx=(res2*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  slpy=(res2*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  centro=_((slpx-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize,(slpy-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize);
  "bpcog: centro-slp2";
  //centro/slp2;
  tmp_array=abs(centro-slp2);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  // to convert in arcsec : (slp2-(y_wfs(nwfs).npix/2.+0.5))*y_wfs(nwfs).pixsize
  /*
  data_gpu = slp2(:184);
  ind_gpu = slp2(185:);
  
  res2=sensors_getdata(g_wfs,0,"bincube");
for (i=1;i<=dimsof(res2)(4);i++){
    res2(*,i)((sort(res2(*,i)))(::-1)(nmax+1:))=0.;
    res2(*,i)((sort(res2(*,i)))(::-1)(:nmax)) -= res2(*,i)((sort(res2(*,i)))(::-1)(nmax));
  }
  data_cpu = res2(*,)(sum,);
  */
  /*
  res2=sensors_getdata(g_wfs,0,"bincube");
  data_cpu = array(0.0f,184);
  ind_cpu = array(0,184);
  pos = 1;
  for(i=1 ; i<=184 ; i++){
    tmp = res2(*,i);
    data_cpu(i) = tmp(sort(tmp))(::-1)(pos)- tmp(sort(tmp))(::-1)(nmax);
    ind_cpu(i) = where(tmp == tmp(sort(tmp))(::-1)(pos));
  }
  */
 
  sensors_compslopes,g_wfs,0,g_rtc,4;
  res2=centroider_getdata(g_rtc,4,"corr");
  res3=sensors_getdata(g_wfs,0,"bincube");
  rcorr = res2 * 0.0f;
  for (i=1;i<=dimsof(res2)(4);i++) rcorr(,,i) = correlfft(res3(,,i),(*(*y_rtc.centroiders)(5).weights)(,,i));
  mcorr = array(0.0f,dimsof(res2)(4));
  for (i=1;i<=dimsof(res2)(4);i++) mcorr(i) = wheremax(rcorr(,,i));
  //for (i=1;i<=dimsof(res2)(4);i++) mcorr(i) = wheremax(res2(8:8+14,8:8+14,i));
  "corr max location: abs(res-mcorr)";
  res=centroider_getdata(g_rtc,4,"corrmax");
  //res-mcorr;
  tmp_array=abs(res-mcorr);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  
  slp = array(0.0f,2*dimsof(res3)(4));
  imat = create_interp_mat(3,3);
  for (i=1;i<=dimsof(res3)(4);i++) {
    tmp = findMaxParaboloid(rcorr(,,i),3,3,imat);
    slp(i) = tmp(1);
    slp(i+dimsof(res3)(4)) = tmp(2);
  }
  
  slp2=sensors_getslopes(g_wfs,0);
  // to convert in arcsec : (slp2-(y_wfs(1).npix/2.+1))*y_wfs(1).pixsize
  "corr: ((slp-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize)-slp2;";
  //((slp-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize)/slp2;
  tmp_array=abs(((slp-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize)-slp2);
  [min(tmp_array),max(tmp_array),avg(tmp_array)];
  //error;
}

/*
tesla 2070
  canaray_elt : 
Average wfs gpu time : 0.1088 s => 9.19118 Hz
  canaray_elt + centroid : 
Average wfs gpu time : 0.1340 s => 7.46269 Hz
  1wfs80x80_1layer :
Average wfs gpu time : 0.0253 s => 39.5257 Hz

gtx 480
  1wfs80x80_1layer :
Average wfs gpu time : 0.0154 s => 64.9351 Hz
  1wfs80x80_1layer + centroid :
Average wfs gpu time : 0.0200 s => 50 Hz
  canary_ngs :
Average wfs gpu time : 0.0016 s => 625 Hz
  canary_ngs + centroid :
Average wfs gpu time : 0.0018 s => 555.556 Hz
  1wfs16x16_1layer :
Average wfs gpu time : 0.0007 s => 1428.57 Hz
  1wfs16x16_1layer + centroid :
Average wfs gpu time : 0.0008 s => 1250 Hz
*/


/*
// results
naos-like
diam=[4.,8.,20.,40.]
["1wfs8x8_1layer","1wfs1x16_1layer","1wfs40x40_1layer","1wfs80x80_1layer",]
4m          - 8m       - 20m      - 40m
naos1 = [0.000454419, 0.00117409, 0.0052497,  0.0253246]
naos2 = [0.000467373, 0.0011995,  0.00536267, 0.0257187]
naos3 = [0.000485916, 0.00128136, 0.00581861, 0.027374]
naos4 = [0.000470978, 0.00120419, 0.00538718, 0.0264424]
naos5 = [0.000496864, 0.00126634, 0.00575755, 0.0314704]
naos6 = [0.000473632, 0.001206,   0.00539444, 0.0263731]

naos-lgs-like
["1wfs8x8_1layer_lgs","1wfs1x16_1layer_lgs","1wfs40x40_1layer_lgs","1wfs80x80_1layer_lgs",]
4m          - 8m       - 20m     - 40m
naosLgs1 = [0.000763339, 0.00229221, 0.0121887, 0.218376]
naosLgs2 = [0.00077659 , 0.00231758, 0.0123006, 0.21854]
naosLgs3 = [0.00079519 , 0.0024001 , 0.0127567, 0.220232]
naosLgs4 = [0.000780166, 0.00232536, 0.01233  , 0.219253]
naosLgs5 = [0.000806483, 0.00238409, 0.0127017, 0.224476]
naosLgs6 = [0.000781823, 0.00232588, 0.0123366, 0.219242]

sphere-like
["sphere4m_1layer","sphere_1layer","sphere20m_1layer","sphere40m_1layer",]
4m          - 8m        - 20m      - 40m
sphere1 = [0.000343465, 0.000730022, 0.00327722, 0.0116978]
sphere2 = [0.000368606, 0.000797327, 0.0036923 , 0.0133116]
sphere3 = [0.000372487, 0.000863224, 0.00408113, 0.0148787]
sphere4 = [0.000374008, 0.000853045, 0.00405633, 0.014764]
sphere5 = [0.000432512, 0.0012298  , 0.00637453, 0.0240691]
sphere6 = [0.000377101, 0.000860718, 0.00405923, 0.0147827]


canary-like
["canary_1layer","canary8m_1layer","canary20m_1layer","canary40m_1layer",]
4m         - 8m       - 20m     - 40m
canary1 = [0.0021683 , 0.00573891, 0.0286964, 0.109454]
canary2 = [0.00222212, 0.00583664, 0.0291034, 0.110974]
canary3 = [0.00230654, 0.00608467, 0.0308139, 0.117758]
canary4 = [0.0022471 , 0.0059265 , 0.0298238, 0.113827]
canary5 = [0.00246509, 0.00669579, 0.03499  , 0.134702]
canary6 = [0.0022514 , 0.0059298 , 0.0297773, 0.113609]


*/

func script_testwfs(filename,&slpgeom1,&slpgeom2,verbose=)
{

  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1subap_1layer_rtc_dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;
  if (y_loop.niter == []) y_loop.niter = 1000;

  // init system
  wfs_init;

  atmos_init;
 
  dm_init;
 
  target_init;

  rtc_init;
  
  if (verbose) write,"... Done with inits !";
  write,"The following objects have been initialized on the GPU :";
  write,"--------------------------------------------------------";
  g_atmos;
  write,"--------------------------------------------------------";
  g_wfs;
  write,"--------------------------------------------------------";
  g_dm;
  write,"--------------------------------------------------------";
  g_target;
  write,"--------------------------------------------------------";
  g_rtc;

  // warming up !
  // move sky
  move_atmos,g_atmos;
  
  // build wfs image
  for (i=1;i<=numberof(y_wfs);i++) {

    sensors_trace,g_wfs,i-1,"atmos",g_atmos;

    sensors_compimg,g_wfs,i-1;
  }

  // init slopes arrays
  tmp = [];
  for (i=1;i<=numberof(y_wfs);i++)
    grow,tmp,sensors_getslopes(g_wfs,i-1);
  slpgeom1 = slpgeom2 = array(0.0,numberof(tmp),y_loop.niter);
  
  rtc_docentroids,g_rtc,g_wfs,0;
  tmp=rtc_getcentroids(g_rtc,0);
  cb_slopes = array(0.0,numberof(tmp),y_loop.niter);

  if (verbose) write,format="Starting loop on : %d iterations\n",y_loop.niter;

  nphase = sensors_getdata(g_wfs,0,"phase")*0.0f;
  sensors_setphase,g_wfs,0,nphase;
  sensors_compimg,g_wfs,0;
  ref_fctn = sensors_getdata(g_wfs,0,"fttotim")(1,,,);
  
  mytime = tic();
  
  for (cc=1;cc<=y_loop.niter;cc++) {

  // move sky
    move_atmos,g_atmos;

    if ((y_wfs != []) && (g_wfs != [])) {
      tmp1 = tmp2 = [];
  // build wfs image & get geom slopes
      for (i=1;i<=numberof(y_wfs);i++) {
        
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;

        if ((!y_wfs(i).openloop) && (g_dm != [])){
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
        }

        sensors_compimg,g_wfs,i-1;

        slopes_geom,g_wfs,i-1,0;

        grow,tmp1,sensors_getslopes(g_wfs,i-1);

        slopes_geom,g_wfs,i-1,1;

        grow,tmp2,sensors_getslopes(g_wfs,i-1);
      }
      
      slpgeom1(,cc) = tmp1;
      slpgeom2(,cc) = tmp2;

      //window,2;pli,roll(sensors_getdata(g_wfs,0,"fttotim")(1,,,1));
      rtc_docentroids,g_rtc,g_wfs,0;
      cb_slopes(,cc) = rtc_getcentroids(g_rtc,0);
      /*
      res3=sensors_getdata(g_wfs,0,"bincube");
      tmp=indices(y_wfs(1).npix);
      slpx=(res3*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res3(*,)(sum,);
      slpy=(res3*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res3(*,)(sum,);
      //centro=_(slpx,slpy);
      centro=_((slpx-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize,(slpy-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize);
      cb_slopes(,cc) = centro;
      */
      /*
      res3=sensors_getdata(g_wfs,0,"fttotim")(1,,,);
      rcorr = array(0.,2*dimsof(res3)(2)-1,2*dimsof(res3)(3)-1,dimsof(res3)(4));
      
      for (i=1;i<=dimsof(res3)(4);i++) {
        rcorr(,,i) = correlfft(res3(,,i),ref_fctn(,,i));
        rcorr(,,i) /= max(rcorr(,,i));
      }
      
      mcorr = array(0.0f,dimsof(res3)(4));
      for (i=1;i<=dimsof(res3)(4);i++) mcorr(i) = wheremax(rcorr(,,i));
  
      slp = array(0.0f,2*dimsof(res3)(4));
      //imat = create_interp_mat(3,3);//*(*y_rtc.centroiders)(1).interpmat;//;
      imat = [[18,-36,18,18,-36,18,18,-36,18],[18,18,18,-36,-36,-36,18,18,18],[27,0,-27,0,0,0,-27,0,27],[-18,0,18,-18,0,18,-18,0,18],[-18,-18,-18,0,0,0,18,18,18],[-12,24,-12,24,60,24,-12,24,-12]]/108.;
      for (i=1;i<=dimsof(res3)(4);i++) {
        tmp = findMaxParaboloid(rcorr(,,i),3,3,imat);
        slp(i) = tmp(1);
        slp(i+dimsof(res3)(4)) = tmp(2);
      }
      cb_slopes(,cc) = ((slp-(dimsof(res3)(2)/2.+1))*y_wfs(1)._qpixsize);
      */
      /*
      tmp=indices(y_wfs(1).npix);
      slpx=(res3*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res3(*,)(sum,);
      slpy=(res3*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res3(*,)(sum,);
      //centro=_(slpx,slpy);
      centro=_((slpx-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize,(slpy-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize);

      slp2=sensors_getslopes(g_wfs,0);
      // to convert in arcsec : (slp2-(y_wfs(1).npix/2.+1))*y_wfs(1).pixsize
      write,format="iter %d\n",cc;
      cb_slopes(,cc);
      ((slp-(y_wfs(1).npix/2.+0.5))*y_wfs(1).pixsize);
      centro;
      
      res=reform(sensors_getdata(g_wfs,0,"phase")(*y_wfs(1)._phasemap+1),y_wfs(1)._pdiam,y_wfs(1)._pdiam,y_wfs(1)._nvalid);
      slpx = float(res(0,avg,)-res(1,avg,));
      slpy = float(res(avg,0,)-res(avg,1,));
      geom=_(slpx,slpy);
      lam_over_D = RASC * y_wfs(1).lambda *1.e-6/ (y_tel.diam/y_wfs(1).nxsub);
      geom *= float(lam_over_D/2/pi);

      geom;
      
      slopes_geom,g_wfs,0,0;
      sensors_getslopes(g_wfs,0);
      */
      /*
      write,format="iter %d\n",cc;
      slopes_geom,g_wfs,0,0;
      sensors_getslopes(g_wfs,0);
      cb_slopes(,cc);
      */
      //error;
     //rtc_docontrol,g_rtc,0,g_dm;

    }
  
    if (verbose) { 
      if (cc % 10 == 0) {
        time_move = tac(mytime)/cc;
        write,format="\v",;
        write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
      }
    }
  }

  if (verbose) {
    write,"... done !";
    write,format="average wfs gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  }

  return cb_slopes;
}

func getslopes_fromsens(slopes,nsensor)
{
  if (nsensor == 1) istart = 1;
  else istart = y_wfs(1:nsensor-1)._nvalid(sum)*2+1;
  if (dimsof(slopes)(1) == 1)
    return slopes(istart:istart+y_wfs(nsensor)._nvalid*2-1);
  if (dimsof(slopes)(1) == 2)
    return slopes(istart:istart+y_wfs(nsensor)._nvalid*2-1,);
}

/*
WFS # | Nsubaps | Nvalid | Npix | Nphase | Nfft | Nrebin | Ntot | Npup  Bias | Coor bias | wcog bias |
    1 |   8x8   |     44 |   15 |     18 |   64 |      5 |   76 |  148   No  |       No  |        No | 
    2 |   8x8   |     44 |   14 |     18 |   64 |      5 |   72 |  148   Yes |      Yes  |       Yes |
    3 |   8x8   |     44 |   15 |     18 |   64 |     10 |  152 |  148   Yes |       No  |        No |
    4 |   8x8   |     44 |   14 |     18 |   64 |     10 |  142 |  148   Yes |      Yes  |       Yes |
    5 |   8x8   |     44 |   15 |     18 |   64 |      9 |  136 |  148   No  |       No  |        No |
    6 |   8x8   |     44 |   14 |     18 |   64 |      9 |  126 |  148   Yes |      Yes  |       Yes | 
    7 |  16x16  |    188 |   15 |      9 |   32 |      2 |   32 |  148   Yes |       No  |        No |
    8 |  15x15  |    162 |   14 |     10 |   32 |      2 |   32 |  148   Yes |      Yes  |       Yes |


 */

/*

  //tests des lgs

nwfs=1
nphase = sensors_getdata(g_wfs,nwfs-1,"phase")*0.0f;
sensors_setphase,g_wfs,nwfs-1,nphase;
sensors_compimg,g_wfs,nwfs-1;
tmp=sensors_getdata(g_wfs,nwfs-1,"fttotim"); // img hr of lgs spot
tmp2 = sensors_getdata(g_wfs,nwfs-1,"bincube"); // binned version
ref_fctn = sensors_getdata(g_wfs,nwfs-1,"lgskern"); // spot not rotated
ampli = sensors_getdata(g_wfs,nwfs-1,"amplifoc"); // subaps psf
spot = ref_fctn(,,16);
psf = roll(ampli(1,,,16));
psf2=spot*0.;
npix = dimsof(psf)(2);
Npix = dimsof(spot)(2);
offs = (Npix-npix)/2;
psf2(offs+1:offs+npix,offs+1:offs+npix) = psf
tmp3=fft(fft(roll(psf2))*fft(spot),-1).re;
pli,tmp(1,,,16);
pli,tmp3;
offset = (y_wfs(nwfs)._Ntot-y_wfs(nwfs)._nrebin*y_wfs(nwfs).npix)/2;
rr = offset+1 : offset + y_wfs(nwfs)._nrebin*y_wfs(nwfs).npix;
tmp4 = float(tmp(1,rr,rr,)(cum,cum,)(::y_wfs(nwfs)._nrebin,::y_wfs(nwfs)._nrebin,)(dif,dif,));
centro=*y_rtc.centroiders;
tmp5 = *centro(1).weights;

 */

func script_modopti(Gmin,Gmax,N)
{
  //activeDevice,1;

  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;
  
  if (filename == []) filename = YOGA_AO_SAVEPATH+"/par/1wfs8x8_1layer_rtc_modopti_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"1pyr32x32_1layer_rtc_dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;
  read_parfile,filename;
  
  G = span(Gmin,Gmax,N);
  strehl = array(0.0f,N,y_loop.niter - 50);
  strehl_opti = array(0.0f,5,y_loop.niter - 50);
  phi_res = array(0.0f,N);
  SRgi = phi_res;
  modG = array(0.0f,5*(y_loop.niter/y_controllers(1).nrec + 1),y_controllers(1).nmodes);
  restore,openb("S2M");
  for(Gi=1;Gi<=N+5;Gi++){
  
  // reading parfile
  read_parfile,filename;
  
  if (Gi<=N){
    y_controllers(1).modopti = 0;
    y_controllers(1).gain = G(Gi);
  }
  else{
    y_controllers(1).modopti = 1;
    cpt_mG=0;
  }
  if (y_loop.niter == []) y_loop.niter = 100000;

  // init system
  wfs_init;

  atmos_init;

  dm_init;

  target_init;
  rtc_init,clean=clean;

  if (Gi==1) {
  write,"The following objects have been initialized on the GPU :";
  write,"--------------------------------------------------------";
  g_atmos;
  write,"--------------------------------------------------------";
  g_wfs;
  write,"--------------------------------------------------------";
  g_target;
  write,"--------------------------------------------------------";
  g_dm;
  write,"--------------------------------------------------------";
  g_rtc;
  }
  
  slopes = array(0.0f,y_wfs(1)._nvalid*2,y_loop.niter);
  /*
                 _         _                   
 _ __ ___   __ _(_)_ __   | | ___   ___  _ __  
| '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \ 
| | | | | | (_| | | | | | | | (_) | (_) | |_) |
|_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/ 
                                        |_|    

   */
  //Align atmos with modal optimisation mode
  if(Gi<=N){
  for (k=1;k<=y_controllers(1).nrec;k++){
    move_atmos,g_atmos;
  }
  }
  for (cc=1;cc<=y_loop.niter;cc++) {
    
    move_atmos,g_atmos;  
 
    if ((y_target != []) && (g_target != [])) {
      // loop on targets
      for (i=1;i<=y_target.ntargets;i++) {
        target_atmostrace,g_target,i-1,g_atmos;
        if (g_dm != []) {
          target_dmtrace,g_target,i-1,g_dm;
        }
      }
    }
    
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        if ((!y_wfs(i).openloop) && (g_dm != [])) {
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
        }

	if(y_wfs(i).type=="cog") {
	  sensors_compimg_tele,g_wfs,i-1;
	} else {
	  sensors_compimg,g_wfs,i-1;
	}
      }     
    }

    if ((y_rtc != []) && (g_rtc != [])
        && (y_wfs != []) && (g_wfs != [])) {
      rtc_docentroids,g_rtc,g_wfs,0;
      slopes(,cc) = rtc_getcentroids(g_rtc,0);
      //rtc_docentroids_geom,g_rtc,g_wfs,0; 
      // compute command and apply
      if (g_dm != []) {
	rtc_docontrol,g_rtc,0;
	rtc_applycontrol,g_rtc,0,g_dm;
      }

    }
    
    if(cc>50 && Gi<=N)
      strehl(Gi,cc-50) = target_getstrehl(g_target,0)(2);
    if(Gi>N){
      strehl_opti(Gi-N,cc-50) = target_getstrehl(g_target,0)(2);
      if((cc-2)/y_controllers(1).nrec >= cpt_mG){
	cpt_mG++;
	modG((Gi-N-1)*dimsof(modG)(2)/5+cpt_mG,) = rtc_getmgain(g_rtc,0); 
      }
    }
  }
  //S2M = rtc_getS2M(g_rtc,0);
  res = S2M(,+) * slopes(+,);
  phi_res(Gi) = sum(res(,rms)^2);
  SRgi(Gi) = avg(strehl(Gi,));
  }
  write,"\n";
  write,"------------------------------------";
  write,"Gain   | Final SR | Avg SR | Min SR | Max SR";
  write,"------------------------------------";
  for(k=1;k<=N;k++){
    write,format=" %5.4f   %2.4f    %2.4f    %2.4f    %2.4f\n",G(k),strehl(k,y_loop.niter-50),avg(strehl(k,)),min(strehl(k,)),max(strehl(k,));
  }

  nref = y_loop.niter/y_controllers(1).nrec;
  write,"\n";
  write,"------------------------------------";
  write,"Avg. Init   | Avg. refresh  |  Avg. SR";
  write,"------------------------------------";
  for(k=1;k<=5;k++){
    write,format="%2.4f        ",avg(modG((k-1)*nref + 1,));
    for(m=1;m<=nref;m++){
      write,format="%2.4f  ",avg(modG((k-1)*nref + m,));
    }
    write,format="%2.4f",avg(strehl_opti(k,));
    write,"\n";
  }

  write,format="Best gain LS : %5.4f for avg SR = %5.4f \n ",G(where(strehl(,avg)==max(strehl(,avg)))),avg(strehl(where(strehl(,avg)==max(strehl(,avg))),)); 
  write,"\n";

  window,0;
  plg,res,G,color="blue";
  plg,SRgi*100,G,color="red";
  error;
}
