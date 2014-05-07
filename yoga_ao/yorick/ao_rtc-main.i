
yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

require,yoga_ao_top+"/yorick/yoga_ao.i";
//require,yoga_ao_top+"/yorick/yaokl.i";
//require,yoga_ao_top+"/ywidgets/widget_wfs.i";

#include "fits-utils.i"
#include "svipc.i"

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data/";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

YOGA_AO_PARPATH = YOGA_AO_SAVEPATH+"par/";
mkdirp,YOGA_AO_PARPATH;

activeDevice,0;
func script_system(filename,verbose=,strehl=,r0=,clean=)
{
  //activeDevice,1;

  shm_init,1234,slots=10;
  sem_init,2345,nums=2;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  if (strehl == []) strehl = 1;
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

  //y_loop.niter=5;
  
  config=[y_wfs.nxsub*y_wfs.npix,y_wfs.nxsub*y_wfs.npix,y_wfs._nvalid*2,y_dm._ntotact(sum),y_wfs.nxsub,y_wfs.nxsub,y_loop.niter,0,0,0](*);
  shm_write,1234,"config",&config;
  
  mc = rtc_getcmat(g_rtc, 0);
  shm_write,1234,"mc",&mc;
  validsubs = transpose(short(*y_wfs(1)._validsubs-1));
  shm_write,1234,"validsubs",&validsubs;
  
  write,"wait for rtc-main-4yorick";
  sem_give,2345,0;
  sem_take,2345,1;

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

    //write, format="loop%d\n", cc;
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
        if(y_wfs(i).type=="sh") {
	  sensors_compimg,g_wfs,i-1;
	} else {
	  sensors_compimg,g_wfs,i-1;
	}
      }
      
    }
    // do centroiding

    image = sensors_getimg(g_wfs,0);
    //image = sensors_getdata(g_wfs,0,"imgtele");
 
    // write, "row1";
    // image(,1);
    // write, "col1";
    // image(1,);
    //pli, image; pause, 100;
    shm_write,1234,"image",&image;
    sem_give,2345,0;
    sem_take,2345,1;

    com = shm_read(1234, "com");
    shm_free, 1234, "com";
    
    if ((y_rtc != []) && (g_rtc != [])
        && (y_wfs != []) && (g_wfs != [])) {
      
      //rtc_docentroids,g_rtc,g_wfs,0;
      // compute command and apply
      //if (g_dm != []) rtc_docontrol,g_rtc,0,g_dm;
      
      if (g_dm != []) dms_comp_shape,g_dm, com;
    }
    
    //write, "**** new loop ****";
    //image(*)(:64);
      
    //controller_getdata(g_rtc, 0, "centroids");
    
    //compare commands computed with yorick and rtc-main
    //com;
    //rtc_getcom(g_rtc,0);
    
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
  }

  //yoga_stop_profiler;
  
  write,"\n done with simulation \n";
  write,format="simulation time : %f sec. per iteration\n",tac(mytime)/y_loop.niter;
  //error;
  //mimg /= y_loop.niter;
  //window,1;fma;pli,mimg; 
  //error;

  shm_cleanup,1234;
  sem_cleanup,2345;
  error;
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
write,"script_system";

