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

YOGA_AO_PARPATH = YOGA_AO_SAVEPATH+"par/";
mkdirp,YOGA_AO_PARPATH;

//activeDevice,1;
func script_system(filename,verbose=,strehl=,r0=,clean=)
{
  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  if (strehl == []) strehl = 0;
  if (r0 == []) r0 = 0;
  if (clean == []) clean = 1;

  if (strehl) {
    extern strehlsp,strehllp,mimg;
  }
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"1pyr32x32_1layer_rtc_dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 100000;

  if (r0 > 0) y_atmos.r0 = r0;

  // init system
  wfs_init;

  atmos_init;

  dm_init;

  target_init;
  rtc_init,clean=clean;

  if (verbose) write,"... Done with inits !";
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
  
  /*
                 _         _                   
 _ __ ___   __ _(_)_ __   | | ___   ___  _ __  
| '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \ 
| | | | | | (_| | | | | | | | (_) | (_) | |_) |
|_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/ 
                                        |_|    

   */
  //yoga_start_profiler;
  
  time_move = 0;
  mytime = tic();
  /*
  mspec=0;
  mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(1));
  nxscreen = dimsof(mscreen)(2);
  tst=(mscreen*0.)(,,-:1:y_loop.niter);
  */

  if (strehl) {
    mimg = 0.; // initializing average image
    strehllp = strehlsp = [];
    write,"\n";
    write,"----------------------------------------------------";
    write,"iter# | S.E. SR | L.E. SR | Est. Rem. | framerate";
    write,"----------------------------------------------------";
  }
  for (cc=1;cc<=y_loop.niter;cc++) {
    
    if (g_target != []) move_sky,g_atmos,g_target;
    else move_atmos,g_atmos;
    /*
    mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(1));
    tst(,,cc)=mscreen;
    mspec += circavg(abs(fft(mscreen)/nxscreen/nxscreen)^2);
    */
   
 
    if ((y_target != []) && (g_target != [])) {
      // loop on targets
      for (i=1;i<=y_target.ntargets;i++) {
        target_atmostrace,g_target,i-1,g_atmos;
        if (g_dm != []) {
          target_dmtrace,g_target,i-1,g_dm;
        }
      }
      //saving average image from target #1
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
      
      // do centroiding
    }
    
    if ((y_rtc != []) && (g_rtc != [])
        && (y_wfs != []) && (g_wfs != [])) {
      rtc_docentroids,g_rtc,g_wfs,0;
      // compute command and apply
      if (g_dm != []) rtc_docontrol,g_rtc,0,g_dm;
    }
    
    if (verbose) {
      subsample=100.;
      if (cc % subsample == 0) {
        timetmp = time_move*(cc-subsample);
        time_move = tac(mytime)/cc;
        timetmp -= time_move*cc;
        if (strehl) {
          //error;
          strehltmp = target_getstrehl(g_target,0);
          grow,strehlsp,strehltmp(1);
          grow,strehllp,strehltmp(2);
          write,format=" %5i    %5.2f     %5.2f     %5.2f s   %5.2f it./s\n",
            cc,strehlsp(0),strehllp(0),(y_loop.niter - cc)*time_move, -1/timetmp*subsample; 
        } else {
          write,format="\v",;
          write,format="\r Estimated remaining time : %.2f s (%.2f it./s)",(y_loop.niter - cc)*time_move, -1/timetmp*subsample;
        }
      }
    } 
  }

  //yoga_stop_profiler;
  
  write,"\n done with simulation \n";
  write,format="simulation time : %f sec. per iteration\n",tac(mytime)/y_loop.niter;
    if (strehl) 
      return strehllp(0);
  //mimg /= y_loop.niter;
  //window,1;fma;pli,mimg; 
  //error;
}

if(batch()) {
  testname=get_argv();
  nb_tests=numberof(testname);
  for(i=2; i<=nb_tests; i++){
    /* valid_rtc stuf
     *
     * pos = strfind("/", testname(i), back=1)(2)+1;
     * output_dir=testname(i);
     * if(pos){ // "/" find
     *   output_dir=strpart(testname(i), pos:);
     * }
     * write, "test de "+testname(i)+", output="+output_dir;
     * script_valid_rtc,testname(i), output=output_dir;
     */
    script_system,testname(i),strehl=1;
  }
} else {
  tmp = get_argv();
  if (numberof(tmp) > 1) {
    if (numberof(tmp) < 3) {
      filename = tmp(2);
      script_system,filename,strehl=1;
    }
    if (numberof(tmp) > 3) {
      filename = tmp(4);
      script_system,filename,strehl=1;
    }
  }
}

func compare_yao(filename)
{
  nr0 = 11;
  
  seeingy=indgen(nr0)/10.+0.2;
  r0y=0.5e-6/seeingy*206265;

  mr0 = array(0.0f,nr0);
  
  mr0(1)=script_system(filename,strehl=1,r0=r0y(1));

  for (i=2;i<=nr0;i++)
    mr0(i)=script_system(filename,strehl=1,r0=r0y(i),clean=0);

  return mr0;
}

func script_profile(filename,verbose=)
{
  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"3wfs8x8_4layers_rtc_2dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 100000;

  // init system
  activeDevice,1;
  
  wfs_init;

  atmos_init;
  
  dm_init;
 
  target_init;
  
  activeDevice,0;
  
  rtc_init;

  activeDevice,1;

  if (verbose) write,"... Done with inits !";
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

  /*
                 _         _                   
 _ __ ___   __ _(_)_ __   | | ___   ___  _ __  
| '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \ 
| | | | | | (_| | | | | | | | (_) | (_) | |_) |
|_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/ 
                                        |_|    

   */
  prof = 0.; // initializing average image
  
  for (cc=1;cc<=y_loop.niter;cc++) {
    tmp = [];
    mytime = tic();
    
    if (g_target != []) move_sky,g_atmos,g_target;
    else move_atmos,g_atmos;
    grow,tmp,tac(mytime);
    
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        grow,tmp,tac(mytime);
        
        if ((!y_wfs(i).openloop) && (g_dm != [])) {
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
          grow,tmp,tac(mytime);
        }
        sensors_compimg_tele,g_wfs,i-1;
        grow,tmp,tac(mytime);
      }
      
      // do centroiding
      if ((y_rtc != []) && (g_rtc != [])) {
        rtc_docentroids,g_rtc,g_wfs,0;
        grow,tmp,tac(mytime);
        // compute command and apply
        if (g_dm != []) rtc_docontrol,g_rtc,0,g_dm;
        grow,tmp,tac(mytime);
      }
    }

    if ((y_target != []) && (g_target != [])) {
      // loop on targets
      for (i=1;i<=y_target.ntargets;i++) {
        target_atmostrace,g_target,i-1,g_atmos;
        grow,tmp,tac(mytime);
        if (g_dm != []) {
          target_dmtrace,g_target,i-1,g_dm;
          grow,tmp,tac(mytime);
        }
      }
      prof += tmp;
    }
    
  }
  write,"\n done with simulation \n";
  write,"simulation profile :\n";
  prof /= y_loop.niter;
  prof;
  return prof;
  //mimg /= y_loop.niter;
  //window,1;fma;pli,mimg; 
}

func script_imat(filename,verbose=,strehl=,r0=,clean=)
{
  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_wfs,y_dm,y_rtc,y_target;
  extern g_atmos,g_target,g_wfs,g_dm,g_rtc,g_target;
  extern ipupil;

  if (verbose == []) verbose = 1;
  if (strehl == []) strehl = 0;
  if (r0 == []) r0 = 0;
  if (clean == []) clean = 1;

  if (strehl) {
    extern strehlsp,strehllp,mimg;
  }
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_geo_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"1pyr32x32_1layer_rtc_dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  // init system
  wfs_init;
  
  dm_init;

  target_init;

  imat = imat_geom();

  correct_dm,imat;

  imat = imat_geom();

 
  if (verbose) write,"... Done with inits !";
  write,"The following objects have been initialized on the GPU :";
  write,"--------------------------------------------------------";
  g_wfs;
  write,"--------------------------------------------------------";
  g_target;
  write,"--------------------------------------------------------";
  g_dm;

  com = random_n(y_dm(1)._ntotact)*y_dm(1).push4imat;
  yoga_setcomm,g_dm,y_dm(1).type,y_dm(1).alt,com;
  yoga_shapedm,g_dm,y_dm(1).type,y_dm(1).alt;
              
  dm_shape = yoga_getdm(g_dm,y_dm(1).type,y_dm(1).alt);

  target_dmtrace,g_target,0,g_dm,1;

  //pli,eclat(target_getimage(g_target,0,"se"));
  pli, yoga_getdm
  error;
}

func script_pyr(filename,verbose=)
{
  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"3wfs8x8_4layers_rtc_2dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 100000;

  // init system
  wfs_init;
             
  atmos_init;
  
  dm_init;

  target_init;
 
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
  
  
  sensors_compimg_tele,g_wfs,0;
 pup=*y_geom._mpupil
phase = pup *0.
phase_rad = pup *0.
  npup = dimsof(phase_rad)(2);
wfs_pupil = roll( pup * exp(1i*phase_rad) );
xy = indices(npup)-(npup+1)/2.;
phase_shift = roll( exp(1i*2*pi*(0.5*xy(,,sum))/npup) );
complex_amplitude = wfs_pupil*phase_shift;
res=sensors_getdata(g_wfs,0,"amplipup");
 complex_amplitude = fft(wfs_pupil*phase_shift,-1);
  complex_amplitude *= *y_wfs(1)._submask;
res=sensors_getdata(g_wfs,0,"amplifoc");
 pyr_npix = y_wfs(1)._Nfft;
 cx = *y_wfs(1)._pyr_cx;
 cy = *y_wfs(1)._pyr_cy;
  xoffset = [1,0,1,0]*pyr_npix;
  yoffset = [1,1,0,0]*pyr_npix;
  reimaged_pupil = array(0.,[3,pyr_npix,pyr_npix,4]);
  pshift = *y_wfs(1)._pyr_offsets;
  for (k=1;k<=numberof(cx);k++) {
    for (i=1;i<=4;i++) {
      ca = roll(complex_amplitude,[xoffset(i)+cx(k),yoffset(i)+cy(k)]);
 small_comp_amp = roll(ca(1:pyr_npix,1:pyr_npix));
 reimaged_pupil(,,i) += abs(fft(small_comp_amp*roll(pshift),1))^2;
    }
  }
res=sensors_getdata(g_wfs,0,"hrimg");

/*
 yoga_oneactu,g_dm,"pzt",0,225,100;
sensors_trace,g_wfs,0,"dm",g_dm;
sensors_compimg,g_wfs,0;
*/
 
 tmp = array(0.,[3,pyr_npix/y_wfs(1)._nrebin,pyr_npix/y_wfs(1)._nrebin,4]);
 for (i=1;i<=4;i++) tmp(,,i) = bin2d(reimaged_pupil(,,i),y_wfs(1)._nrebin);
 reimaged_pupil = tmp;

     pup = *y_geom._mpupil;
    pupreb = bin2d(pup*1.,y_wfs(1).npix)/y_wfs(1).npix^2.;
    wsubok = where(pupreb>=y_wfs(1).fracsub);
    
  pixels = reimaged_pupil(*,)(wsubok,);

  sigx = (pixels(,[2,4])(,sum)-pixels(,[1,3])(,sum))/(pixels(,sum)+1e-6);
  sigy = (pixels(,[3,4])(,sum)-pixels(,[1,2])(,sum))/(pixels(,sum)+1e-6);

res=sensors_getdata(g_wfs,0,"bincube");

  pixels = res(*,)(wsubok,);

  sigx2 = (pixels(,[2,4])(,sum)-pixels(,[1,3])(,sum))/(pixels(,sum)+1e-6);
  sigy2 = (pixels(,[3,4])(,sum)-pixels(,[1,2])(,sum))/(pixels(,sum)+1e-6);


 error;

  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    
    move_atmos,g_atmos;
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        sensors_compimg_tele,g_wfs,i-1;
      }
    }
    error;
    window,1;fma;pli,sensors_getdata(g_wfs,0,"phase");
    window,2;fma;pli,sensors_getdata(g_wfs,0,"bincube")(,,1);

    hitReturn;
  }
}

/*
activeDevice,1;res1=script_profile("1wfs8x8_1layer_rtc_dm.par");res2=script_profile("1wfs16x16_1layer_rtc_dm.par");res3=script_profile("1wfs40x40_1layer_rtc_dm.par");res4=script_profile("1wfs60x60_1layer_rtc_dm.par");res5=script_profile("1wfs80x80_1layer_rtc_dm.par");
prof=[res1,res2,res3,res4,res5];
diam=[4.,8.,20.,30.,40.];

1.70972e-05,1.49579e-05,0.000176983,6.05841e-05,0.000225203,1.83723e-05, 1.33798e-05
1.73852e-05,1.50435e-05,0.000176213,2.17693e-05,0.00012555,1.78308e-05, 1.32682e-05

2.77002e-05,2.27838e-05,0.000238173,4.21488e-05,0.0393674,2.89803e-05, 2.11103e-05

30m :
fast
[0.00022509,0.000241423,0.000254502,0.000359016,0.000382632,0.0456178,0.0456362,0.0456483]
1.63331e-05,1.30796e-05,0.000104513,2.36161e-05,0.0452352,1.83451e-05,1.20678e-05
uva
[0.000243724,0.000260524,0.000273978,0.000385617,0.000410179,0.0457285,0.0457473,0.0457598]
[1.63908e-05,1.30837e-05,0.000100648,2.34144e-05,0.0453188,1.84021e-05,1.19252e-05]

*/
func script_relax(filename,verbose=)
{
  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  
  if (filename == []) filename = YOGA_AO_PARPATH+"relax.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc_dm_lgs.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"3wfs8x8_4layers_rtc_2dm.par";
  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 100000;

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
  g_target;
  write,"--------------------------------------------------------";
  g_dm;
  write,"--------------------------------------------------------";
  g_rtc;

  /*
                 _         _                   
 _ __ ___   __ _(_)_ __   | | ___   ___  _ __  
| '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \ 
| | | | | | (_| | | | | | | | (_) | (_) | |_) |
|_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/ 
                                        |_|    

   */
  cb_mir        = yoga_getdm(g_dm,"pzt",0.)(,,-:1:1); 
  cb_tt         = yoga_getdm(g_dm,"tt",0.)(,,-:1:1); 
  cb_phase      = sensors_getdata(g_wfs,0,"phasetele")(,,-:1:1);
  cb_phase_nott = sensors_getdata(g_wfs,0,"phasetele")(,,-:1:1);
  
  tt_influ      = (*y_dm(2)._influ)(8:-7,8:-7,);
  pup           = *y_geom._mpupil;
  val_pts       = where(pup);
  tmp           = tt_influ(*,)(val_pts,);
  mat_pass      = (tmp(+,)*tmp(+,))(,+)*tmp(,+);
  
  cb_centro     = sensors_getslopes(g_wfs,0)(,,-:1:1);
  slopes_geom,g_wfs,0,1;
  cb_slps       = sensors_getslopes(g_wfs,0)(,,-:1:1);

  
  mytime = tic();
  pli,cb_mir(,,1);
  
  for (cc=1;cc<=y_loop.niter;cc++) {
    
    if (g_target != []) move_sky,g_atmos,g_target;
    else move_atmos,g_atmos;
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        if ((!y_wfs(i).openloop) && (g_dm != [])) {
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
        }
        sensors_compimg_tele,g_wfs,i-1;
      }
      
      // do centroiding
      if ((y_rtc != []) && (g_rtc != [])) {
        rtc_docentroids,g_rtc,g_wfs,0;
        // compute command and apply
        if (g_dm != []) rtc_docontrol,g_rtc,0,g_dm;
      }
    }
    
    tmp  = yoga_getdm(g_dm,"pzt",0.);
    grow,cb_mir,tmp;
    tmp  +=  yoga_getdm(g_dm,"tt",0.);
    //grow,cb_tt,tmp;
    
    tmp2 = sensors_getdata(g_wfs,0,"phasetele");
    grow,cb_phase,tmp2;

    nx = dimsof(tmp2)(2);
    Nx = dimsof(tmp)(2);
    
    fma;pli,tmp((Nx-nx)/2+1:Nx-(Nx-nx)/2,(Nx-nx)/2+1:Nx-(Nx-nx)/2)+tmp2;
    //grow,cb_phase_nott,tmp-(tmp(*)(val_pts)(+)*mat_pass(,+))(+)*tt_influ(,,+);

    

    /*
    grow,cb_centro,sensors_getslopes(g_wfs,0);
    slopes_geom,g_wfs,0,1;
    grow,cb_slps,sensors_getslopes(g_wfs,0);
    */
    yoga_resetdm,g_dm,"pzt",0.;
    yoga_resetdm,g_dm,"tt",0.;
    
    
     if (verbose) { 
      if (cc % 10 == 0) {
        time_move = tac(mytime)/cc;
        write,format="\v",;
        write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
      }
     } 
  }

  write,"\n done with simulation \n";
  write,format="simulation time : %f sec. per iteration\n",tac(mytime)/y_loop.niter;
  cb_mir        = cb_mir(8:-7,8:-7,2:);
  cb_tt         = cb_tt(8:-7,8:-7,2:);
  cb_phase      = cb_phase(,,2:);
  cb_phase_nott = cb_phase_nott(,,2:);
  cb_centro     = cb_centro(,2:);
  cb_slps       = cb_slps(,2:);

  error;
  fits_write,strip_file_extension(filename,".par")+"_cb_mir.fits",cb_mir;
  fits_write,strip_file_extension(filename,".par")+"_cb_tt.fits",cb_tt;
  fits_write,strip_file_extension(filename,".par")+"_cb_phase.fits",cb_phase;
  fits_write,strip_file_extension(filename,".par")+"_cb_phase_nott.fits",cb_phase_nott;
  fits_write,strip_file_extension(filename,".par")+"_cb_centro.fits",cb_centro;
  fits_write,strip_file_extension(filename,".par")+"_cb_slps.fits",cb_slps;

  g_atmos = g_wfs = g_target = g_dm = g_rtc = [];
  y_geom = y_tel = y_atmos = y_target = y_loop = y_wfs =
    y_dm = y_controllers = y_centroiders = y_rtc = [];

}

func script_valid_rtc(filename,verbose=,strehl=,r0=,clean=, output=)
{
  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  if (strehl == []) strehl = 0;
  if (r0 == []) r0 = 0;
  if (clean == []) clean = 1;

  if (strehl) {
    extern strehlsp,strehllp,mimg;
  }
  
  if (filename == []) filename = YOGA_AO_PARPATH+"/1wfs8x8_1layer_rtc_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"/1pyr32x32_1layer_rtc_dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 100000;

  if (r0 > 0) y_atmos.r0 = r0;

  // init system
  wfs_init;

  atmos_init;
  
  dm_init;

  target_init;

  rtc_init,clean=clean;

  if (verbose) write,"... Done with inits !";
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
  
  y_loop.niter = 1000;
  
  /*
                 _         _                   
 _ __ ___   __ _(_)_ __   | | ___   ___  _ __  
| '_ ` _ \ / _` | | '_ \  | |/ _ \ / _ \| '_ \ 
| | | | | | (_| | | | | | | | (_) | (_) | |_) |
|_| |_| |_|\__,_|_|_| |_| |_|\___/ \___/| .__/ 
                                        |_|    

   */
  
  time_move = 0;
  mytime = tic();
  /*
  mspec=0;
  mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(1));
  nxscreen = dimsof(mscreen)(2);
  tst=(mscreen*0.)(,,-:1:y_loop.niter);
  */

  if (strehl) {
    mimg = 0.; // initializing average image
    strehllp = strehlsp = [];
    airy = roll(abs(fft(*y_geom._ipupil*exp(*y_geom._ipupil*1i*0.)))^2)/numberof(*y_geom._ipupil);
    sairy = max(airy);
    write,"\n";
    write,"----------------------------------------------------";
    write,"iter# | S.E. SR | L.E. SR | Est. Rem. | framerate";
    write,"----------------------------------------------------";
  }

  nb_pix = y_wfs(1).nxsub* y_wfs(1).npix;
  img_cube=array(0.f,nb_pix,nb_pix,y_loop.niter);
  binimg_cube=array(0.f,y_wfs(1).npix,y_wfs(1).npix,y_wfs._nvalid(1),y_loop.niter);

  centro_cube=array(0.f, numberof(rtc_getcentroids(g_rtc, 0)), y_loop.niter);
  slopes_cube=centro_cube;

  com_cube=array(0.f, numberof(controller_getdata(g_rtc, 0, "com")), y_loop.niter);

  yoga_start_profile;
  for (cc=1;cc<=y_loop.niter;cc++) {
    
    if (g_target != []) move_sky,g_atmos,g_target;
    else move_atmos,g_atmos;
    /*
    mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(1));
    tst(,,cc)=mscreen;
    mspec += circavg(abs(fft(mscreen)/nxscreen/nxscreen)^2);
    */
    
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        if ((!y_wfs(i).openloop) && (g_dm != [])) {
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
        }
        sensors_compimg_tele,g_wfs,i-1;
    img_cube(,,cc) = sensors_getdata(g_wfs,i-1,"imgtele");
    binimg_cube(,,,cc) = sensors_getdata(g_wfs,i-1,"bincube");
      }
      // do centroiding
    }
    
    if ((y_rtc != []) && (g_rtc != [])
        && (y_wfs != []) && (g_wfs != [])) {
      rtc_docentroids,g_rtc,g_wfs,0;
      // compute command and apply
      centro_cube(,cc) = rtc_getcentroids(g_rtc, 0);
      slopes_cube(,cc) = sensors_getdata(g_wfs, 0, "slopes");
      if (g_dm != []) rtc_docontrol,g_rtc,0,g_dm;
      com_cube(,cc) = controller_getdata(g_rtc, 0, "com");
    }
    
    if ((y_target != []) && (g_target != [])) {
      // loop on targets
      for (i=1;i<=y_target.ntargets;i++) {
        target_atmostrace,g_target,i-1,g_atmos;
        if (g_dm != []) {
          target_dmtrace,g_target,i-1,g_dm;
        }
      }
      //saving average image from target #1
    }

    if (verbose) {
      subsample=100.;
      if (cc % subsample == 0) {
        timetmp = time_move*(cc-subsample);
        time_move = tac(mytime)/cc;
        timetmp -= time_move*cc;
        if (strehl) {
          strehltmp = target_getstrehl(g_target,0);
          grow,strehlsp,strehltmp(1);
          grow,strehllp,strehltmp(2);
          write,format=" %5i    %5.2f     %5.2f     %5.2f s   %5.2f it./s\n",
            cc,strehlsp(0),strehllp(0),(y_loop.niter - cc)*time_move, -1/timetmp*subsample; 
        } else {
          write,format="\v",;
          write,format="\r Estimated remaining time : %.2f s (%.2f it./s)",(y_loop.niter - cc)*time_move, -1/timetmp*subsample;
        }
      }
    } 
  }

  yoga_stop_profile;
  
  write,"\n done with simulation \n";
  write,format="simulation time : %f sec. per iteration\n",tac(mytime)/y_loop.niter;
  mc = rtc_getcmat(g_rtc, 0);
  //error;
  if(output!=[]){
    mkdirp,output;
    writefits, output+"/img_cube.fits", img_cube;
    writefits, output+"/binimg_cube.fits", binimg_cube;
    writefits, output+"/slopes_cube.fits", slopes_cube;
    //fits_write, output+"/centro_cube.fits", centro_cube;
    writefits, output+"/com_cube.fits", com_cube;
    writefits, output+"/mc.fits", mc;
    writefits, output+"/config.fits", transpose(*y_wfs._validsubs(1)-1), ["NPIX", "NVALID", "NXSUB"], [y_wfs.npix(1), y_wfs._nvalid(1),y_wfs.nxsub(1)];
  }

  //return strehllp(0);
  //mimg /= y_loop.niter;
  //window,1;fma;pli,mimg; 
  //error;
}



func autocorr(tmp)
{
  tmp1 = 0.;
  if (dimsof(tmp)(1) == 1) {
    nx = numberof(tmp);
    aa = array(0.,2*nx);
    aa(1:nx) = tmp-tmp(*)(avg);
    bb = aa;
    tmp1 = roll(fft(fft(aa)*conj(fft(bb)),-1).re)(2:);
  }
  
  return tmp1;
}

func script_pyr_diff(filename,verbose=)
{
  //activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs, g_dm, g_rtc;
  extern ipupil;

  if (verbose == []) verbose = 1;
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1pyr16x16_1layer_rtc_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"3wfs8x8_4layers_rtc_2dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 100000;

  // init system
  wfs_init;
             
  atmos_init;
  
  dm_init;

  target_init;
 
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
  
  error;
move_sky,g_atmos,g_target;
sensors_trace,g_wfs,0,"atmos",g_atmos;
sensors_compimg,g_wfs,0;
 res=sensors_getdata(g_wfs,0,"amplifoc");
 pli,roll(res(1,,)^2+res(2,,)^2)^0.2;

    xy = indices(Nfft)-(Nfft+1)/2.;
 npup = y_wfs(1)._Ntot/2;

 xy = indices(npup)/(npup)-0.5;
 pshift = roll(exp(1i*2*pi*xy(,,sum)));

 phasemap_diff = array(complex,2*npup,2*npup);
 phasemap_diff(1:npup,1:npup)   = pshift;
 phasemap_diff(npup+1:,1:npup)  = transpose(pshift)(::-1,);
 phasemap_diff(npup+1:,npup+1:)  = pshift(::-1,::-1);
 phasemap_diff(1:npup,npup+1:) = transpose(pshift)(,::-1);


  sensors_compimg,g_wfs,0;
 pup=*y_geom._mpupil
phase = pup *0.
phase_rad = pup *0.
  npup = dimsof(phase_rad)(2);
wfs_pupil = roll( pup * exp(1i*phase_rad) );
xy = indices(npup)-(npup+1)/2.;
phase_shift = roll( exp(1i*2*pi*(0.5*xy(,,sum))/npup) );
complex_amplitude = wfs_pupil*phase_shift;
res=sensors_getdata(g_wfs,0,"amplipup");
 complex_amplitude = fft(wfs_pupil*phase_shift,-1);
  complex_amplitude *= *y_wfs(1)._submask;
res=sensors_getdata(g_wfs,0,"amplifoc");
 pyr_npix = y_wfs(1)._Nfft;
 cx = *y_wfs(1)._pyr_cx;
 cy = *y_wfs(1)._pyr_cy;
  xoffset = [1,0,1,0]*pyr_npix;
  yoffset = [1,1,0,0]*pyr_npix;
  reimaged_pupil = array(0.,[3,pyr_npix,pyr_npix,4]);
  pshift = *y_wfs(1)._pyr_offsets;
  for (k=1;k<=numberof(cx);k++) {
    for (i=1;i<=4;i++) {
      ca = roll(complex_amplitude,[xoffset(i)+cx(k),yoffset(i)+cy(k)]);
 small_comp_amp = roll(ca(1:pyr_npix,1:pyr_npix));
 reimaged_pupil(,,i) += abs(fft(small_comp_amp*roll(pshift),1))^2;
    }
  }
res=sensors_getdata(g_wfs,0,"hrimg");

/*
 yoga_oneactu,g_dm,"pzt",0,225,100;
sensors_trace,g_wfs,0,"dm",g_dm;
sensors_compimg,g_wfs,0;
*/
 
 tmp = array(0.,[3,pyr_npix/y_wfs(1)._nrebin,pyr_npix/y_wfs(1)._nrebin,4]);
 for (i=1;i<=4;i++) tmp(,,i) = bin2d(reimaged_pupil(,,i),y_wfs(1)._nrebin);
 reimaged_pupil = tmp;

     pup = *y_geom._mpupil;
    pupreb = bin2d(pup*1.,y_wfs(1).npix)/y_wfs(1).npix^2.;
    wsubok = where(pupreb>=y_wfs(1).fracsub);
    
  pixels = reimaged_pupil(*,)(wsubok,);

  sigx = (pixels(,[2,4])(,sum)-pixels(,[1,3])(,sum))/(pixels(,sum)+1e-6);
  sigy = (pixels(,[3,4])(,sum)-pixels(,[1,2])(,sum))/(pixels(,sum)+1e-6);

res=sensors_getdata(g_wfs,0,"bincube");

  pixels = res(*,)(wsubok,);

  sigx2 = (pixels(,[2,4])(,sum)-pixels(,[1,3])(,sum))/(pixels(,sum)+1e-6);
  sigy2 = (pixels(,[3,4])(,sum)-pixels(,[1,2])(,sum))/(pixels(,sum)+1e-6);


 error;

  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    
    move_atmos,g_atmos;
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        sensors_compimg_tele,g_wfs,i-1;
      }
    }
    error;
    window,1;fma;pli,sensors_getdata(g_wfs,0,"phase");
    window,2;fma;pli,sensors_getdata(g_wfs,0,"bincube")(,,1);

    hitReturn;
  }
}

//script_system, YOGA_AO_PARPATH+"1wfs40x40_1layer_rtc_dm.par",strehl=1;
