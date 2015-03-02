#include"svipc.i"
shm_key=1234
sem_key=2345
write, "je suis dans Yorick!"

yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

require,yoga_ao_top+"/yorick/yoga_ao.i";
//require,"prana.i";

#include "fits-utils.i"

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data/";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

YOGA_AO_PARPATH = YOGA_AO_SAVEPATH+"par/";
mkdirp,YOGA_AO_PARPATH;

activeDevice,0;
func script_system(filename,verbose=,strehl=,r0=,clean=,force_niter=)
{
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  yoga_stop_profile;
  
  activeDeviceForce,0;
  if (verbose == []) verbose = 1;
  if (strehl == []) strehl = 1;
  if (r0 == []) r0 = 0;
  if (clean == []) clean = 1;

  if (strehl) {
    extern strehlsp,strehllp,mimg;
  }
  
  if (filename == []) filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"1wfs40x40_1layer_rtc_dm.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"1pyr32x32_1layer_rtc_dm.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  // reading parfile
  write, "reading filename:"+filename
  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 100000;

  if (r0 > 0) y_atmos.r0 = r0;

  // init system
  wfs_init;

  tab_type= strchar(y_wfs.type);
  nb_wfs=[numberof(y_wfs)];
  tab_nxsub=y_wfs.nxsub;
  tab_nvalid=y_wfs._nvalid;
  tab_npix=y_wfs.npix;
  tab_nphase=y_wfs._pdiam;
  tab_nrebin=y_wfs._nrebin;
  tab_nfft=y_wfs._Nfft;
  tab_ntot=y_wfs._Ntot;
  tab_npup=[y_geom._n];
  tab_pdiam=float(y_wfs._subapd);
  tab_nphot=float(y_wfs._nphotons);
  tab_lgs=y_wfs.gsalt>0;

  //shm_write, shm_key, "type", &tab_type;
  shm_write, shm_key, "nwfs", &nb_wfs;
  shm_write, shm_key, "nbsub", &tab_nxsub;
  shm_write, shm_key, "nvalid", &tab_nvalid;
  shm_write, shm_key, "npix", &tab_npix;
  shm_write, shm_key, "nphase", &tab_nphase;
  shm_write, shm_key, "nrebin", &tab_nrebin;
  shm_write, shm_key, "nfft", &tab_nfft;
  shm_write, shm_key, "ntot", &tab_ntot;
  shm_write, shm_key, "npup", &tab_npup;
  shm_write, shm_key, "pdiam", &tab_pdiam;
  shm_write, shm_key, "nphot", &tab_nphot;
  shm_write, shm_key, "lgs", &tab_lgs;

  sem_give,sem_key,0;

  size=(y_geom._n)(-::numberof(y_wfs)-1);
  shm_write, shm_key, "xpos", &(y_wfs.xpos);
  shm_write, shm_key, "ypos", &(y_wfs.ypos);
  shm_write, shm_key, "lambda", &(y_wfs.lambda);
  shm_write, shm_key, "gsmag", &(y_wfs.gsmag);
  shm_write, shm_key, "size", &size;
  shm_write, shm_key, "noise", &(y_wfs.noise);
  
  sem_give,sem_key,0;

  hrmap=int([*y_wfs(1)._hrmap]);
  validsubsx= int((*y_wfs(1)._validsubs)(1,)-1);
  validsubsy= int((*y_wfs(1)._validsubs)(2,)-1);
  istart=int((*y_wfs(1)._istart)+1);
  jstart=int((*y_wfs(1)._jstart)+1);
  fluxPerSub=(*y_wfs(1)._fluxPerSub)(where(*y_wfs(1)._isvalid));

  shm_write, shm_key, "phasemap",y_wfs(1)._phasemap;
  shm_write, shm_key, "hrmap", &hrmap;
  shm_write, shm_key, "binmap", y_wfs(1)._binmap;
  shm_write, shm_key, "offsets", y_wfs(1)._halfxy;
  shm_write, shm_key, "pupil", y_geom._mpupil;
  shm_write, shm_key, "fluxPerSub", &fluxPerSub
  shm_write, shm_key, "isvalid", y_wfs(1)._isvalid;
  shm_write, shm_key, "validsubsx",&validsubsx;
  shm_write, shm_key, "validsubsy", &validsubsy;
  shm_write, shm_key, "istart", &istart;
  shm_write, shm_key, "jstart", &jstart;
  shm_write, shm_key, "ftkernel", y_wfs(1)._ftkernel;
  sem_give,sem_key,0;
  
}

func ajeter( void ){
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

  use_gpu=int([1]);
  p_rtc = prana_rtc(use_gpu, y_wfs(1).nxsub*y_wfs(1).npix,y_wfs(1).nxsub*y_wfs(1).npix, 
    y_wfs(1)._nvalid*2, y_dm._ntotact(sum),y_wfs(1).nxsub,y_wfs(1).nxsub, y_controllers(1).gain, 
    transpose(short(*y_wfs(1)._validsubs-1)),rtc_getcmat(g_rtc, 0));
  use_p_rtc=1;
  //p_rtc=[];
    
  if(p_rtc!=[] && use_p_rtc)
    prana_start, p_rtc;
  
  //prana_set_image, p_rtc, g_wfs,0; //yoga_obj(float(random(64, 64)));
  //p_rtc=[];
  
/*               _         _                   
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
  
  if(force_niter!=[])
    y_loop.niter = force_niter;
  loop_timer = array(float, [2, y_loop.niter, 7]);
  
  //yoga_start_profile;
  
  for (cc=1;cc<=y_loop.niter;cc++) {
    // error;    
    //write, format="beg loop%d\n", cc;

    if (g_target != []) move_sky,g_atmos,g_target;
    else move_atmos,g_atmos;

    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        if ((!y_wfs(i).openloop) && (g_dm != [])) {
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
        }
        if(y_wfs(i).type=="sh") {
          sensors_compimg,g_wfs,i-1;
        } else {
          sensors_compimg,g_wfs,i-1;
        }
      } 
    }

    if(p_rtc!=[] && use_p_rtc) {
      //pli,sensors_getimg(g_wfs,0); pause, 100;
      //image = sensors_getdata(g_wfs,0,"imgtele");
 
      // write, "row1";
      // Image(,1);
      // write, "col1";
      // image(1,);
      //pli, image; pause, 100;
    
      //write, "debug: prana_set_image"
      loop_timer(cc, 1) = prana_set_image(p_rtc, g_wfs,0);

      //d_com = yoga_obj(array(0.f, y_dm._ntotact(sum)));

      //write, "debug: prana_apply_command"
      loop_timer(cc, 2:) = prana_apply_commands(p_rtc, g_dm);
      //write, "**** new loop ****";
      //image(*)(:64);
      
      //controller_getdata(g_rtc, 0, "centroids");
    
      //compare commands computed with yorick and rtc-main
      //com;
      //rtc_getcom(g_rtc,0);
    } else {
      if ((y_rtc != []) && (g_rtc != [])
          && (y_wfs != []) && (g_wfs != [])) {
        rtc_docentroids,g_rtc,g_wfs,0;
        // compute command and apply
        if (g_dm != []) rtc_docontrol,g_rtc,0,g_dm;
      }
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
  } //fin main loop

  yoga_stop_profile;
  
  write,"\n done with simulation \n";
  write,format="simulation time : %f sec. per iteration\n",tac(mytime)/y_loop.niter;
  //error;
  //mimg /= y_loop.niter;
  //window,1;fma;pli,mimg; 
  //error;
  if(p_rtc!=[] && use_p_rtc)
    prana_stop, p_rtc;

  return loop_timer;
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
/*
write,"script_system";
write,"loop_timer=script_system(YOGA_AO_PARPATH+\"1wfs40x40_1layer_rtc_dm.par\", force_niter=100000)"
write,"loop_timer=script_system(YOGA_AO_PARPATH+\"1wfs80x80_1layer_rtc_dm.par\", force_niter=100000)"
*/