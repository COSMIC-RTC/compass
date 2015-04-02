/*
*/
// Environment check
yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

/*
plug_path = yoga_ao_top;
grow,plug_path,"/home/brujo/yorick/yoga/trunk";
plug_dir,plug_path;
*/

// load necessary files : yorick-python wrapper and styc_utils
require,yoga_ao_top+"/ywidgets/pyk.i";
require,yoga_ao_top+"/ywidgets/widgets_utils.i";
require,yoga_ao_top+"/ywidgets/atmos_utils.i";
require,yoga_ao_top+"/ywidgets/target_utils.i";
require,yoga_ao_top+"/ywidgets/wfs_utils.i";
require,yoga_ao_top+"/ywidgets/rtc_utils.i";
require,yoga_ao_top+"/ywidgets/dm_utils.i";
require,yoga_ao_top+"/yorick/yoga_ao.i";


func load_default_system(tconf,teldiam)
{
  if (tconf == 1) {
    load_default_atmos,1,,teldiam;
    load_default_wfs,1,teldiam;
    
    load_default_dm,1,max(y_wfs.nxsub);
    load_default_target,1;
    load_default_centro,"cog",;
    load_default_control,"ls";
  }
  if (tconf == 2) {
    load_default_atmos,1,,teldiam;
    load_default_wfs,2,teldiam;
    load_default_dm,1,max(y_wfs.nxsub);
    load_default_target,1;
    load_default_centro,"cog",;
    load_default_control,"ls";
  }
  if (tconf == 3) {
    load_default_atmos,1,,teldiam;
    load_default_wfs,3,teldiam;
    load_default_dm,1,max(y_wfs.nxsub);
    load_default_target,1;
    load_default_centro,"pyr",;
    load_default_control,"ls";
  }
  if (tconf == 4) {
    load_default_atmos,1,,teldiam;
    load_default_wfs,4,teldiam;
    load_default_dm,1,max(y_wfs.nxsub);
    load_default_target,1;
    load_default_centro,"pyr",;
    load_default_control,"ls";
  }
  if (tconf == 5) {
    load_default_atmos,1,,teldiam;
    load_default_wfs,5,teldiam;
    load_default_dm,1,max(y_wfs.nxsub);
    load_default_target,1;
    load_default_centro,"cog",;
    load_default_control,"ls";
  }
  if (tconf == 6) {
    load_default_atmos,1,,teldiam;
    load_default_wfs,6,teldiam;
    load_default_dm,1,max(y_wfs.nxsub);
    load_default_target,1;
    load_default_centro,"cog",;
    load_default_control,"ls";
  }
  if (tconf == 7) {
    load_default_atmos,4,,teldiam;
    load_default_wfs,7,teldiam;
    load_default_dm,2,max(y_wfs.nxsub);
    load_default_target,3;
    load_default_centro,"cog",;
    load_default_control,"ls";
  }
  
}



func load_parfile(parfile,filename)
{
  extern yoga_parfilename;
  extern y_wfs,y_atmos,y_rtc,y_loop,y_tel,y_geom,y_target;
  
  yoga_parfilename = filename;

  y_wfs = y_atmos = y_rtc = y_loop = y_tel = y_geom = y_target = [];
  g_wfs = g_atmos = g_rtc = g_target = [];
  
  read_parfile,parfile;

  //y_atmos.nscreens;
  if (y_tel != []) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('teldiam').set_value(%f)",y_tel.diam);
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('cobs').set_value(%f)",y_tel.cobs);
  }
  if (y_atmos != []) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('nlayers').set_value(%d)",y_atmos.nscreens);
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('r0').set_value(%f)",y_atmos.r0);
    pyk,swrite(format=wfs_disp._cmd+"y_layer_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('layer_select').clear()");
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('layer_select').insert_text(%d,'%s')",0,
                 swrite(format="Layer # %d",numberof(y_wfs)-cc+1));
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('layer_select').set_active(%d)",0);
  }
  if (y_wfs != []) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('nwfs').set_value(%d)",numberof(y_wfs));
    pyk,swrite(format=wfs_disp._cmd+"y_wfs_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('wfs_select').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",0);
    //if (y_rtc != []) 
    // pyk,swrite(format=ao_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",0);
    //else 
    //  pyk,swrite(format=ao_disp._cmd+"glade.get_widget('default_centro').set_active(%d)",0);
  }
  if (y_target != []) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('ntargets').set_value(%d)",y_target.ntargets);
    pyk,swrite(format=wfs_disp._cmd+"y_target_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('target_select').clear()");
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('target_select').set_active(%d)",0);
  }
  if (y_dm != []) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('ndm').set_value(%d)",numberof(y_dm));
    pyk,swrite(format=wfs_disp._cmd+"y_dm_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('target_select').clear()");
    for (cc=1;cc<=numberof(y_dm);cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('dm_select').insert_text(%d,'%s')",0,
                 swrite(format="DM # %d",numberof(y_dm)-cc+1));
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('dm_select').set_active(%d)",0);
  }
  if (y_rtc != []) {
    pyk,swrite(format=wfs_disp._cmd+"y_centro_clear(%d)",1);    
    //pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').insert_text(%d,'%s')",0,swrite(format="centro # %d",numberof(y_wfs)-cc+1));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",0);
    type = (*y_rtc.centroiders)(1).type;
    if (type == "wcog") {
      if (type_fct == "gauss") {
        pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('type_func').set_active(%d)",0);
      }
    }
    update_control_prop; 
  }
  
}

func ao_win_init(pid1)
{
  extern ao_disp;

  if (numberof(*ao_disp._xid)==0) ao_disp._xid = &[0,0,0,0,0,0];
  (*ao_disp._xid)(1) = pid1;
  
  if (catch(0x08)) {
    (*ao_disp._winits)(1) = 1;
  }
  
  if (!(*ao_disp._winits)(1)) {
    stylename = yoga_ao_top+"/gs/yoga_noaxis.gs";
    window,(*ao_disp._wins)(1),dpi=ao_disp._defaultdpi,width=0,height=0,
      xpos=-4,ypos=-4,style=stylename,parent=pid1;
    limits,square=1;
  }

  (*ao_disp._winits)(1) = 1;
}

func turbu_win_init(pid2)
{
  extern ao_disp;

  if (numberof(*ao_disp._xid)==0) ao_disp._xid = &[0,0,0,0,0,0];
  (*ao_disp._xid)(2) = pid2;
  
  if (catch(0x08)) {
    (*ao_disp._winits)(2) = 1;
  }
  
  if (!(*ao_disp._winits)(2)) {
    stylename = yoga_ao_top+"/gs/yoga_ao.gs";
    window,(*ao_disp._wins)(2),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid2;
    limits,square=1;
    palette,"gray.gp";
  }

  (*ao_disp._winits)(2) = 1;
}
    
func lgs_win_init(pid3)
{
  extern ao_disp;

  if (numberof(*ao_disp._xid)==0) ao_disp._xid = &[0,0,0,0,0,0];
  (*ao_disp._xid)(3) = pid3;
  
  if (catch(0x08)) {
    (*ao_disp._winits)(3) = 1;
  }
  
  if (!(*ao_disp._winits)(3)) {
    stylename = yoga_ao_top+"/gs/yoga_ao.gs"; 
    window,(*ao_disp._wins)(3),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid3;
    limits,square=1;
    palette,"gray.gp"; 
  }

  (*ao_disp._winits)(3) = 1;
}
    
func centro_win_init(pid4)
{
  extern ao_disp;
  
  if (numberof(*ao_disp._xid)==0) ao_disp._xid = &[0,0,0,0,0,0];
  (*ao_disp._xid)(4) = pid4;
  
  if (catch(0x08)) {
    (*ao_disp._winits)(4) = 1;
  }
  
  if (!(*ao_disp._winits)(4)) {
    stylename = yoga_ao_top+"/gs/yoga_noaxis.gs"; 
    window,(*ao_disp._wins)(4),dpi=27,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid4;
    limits,square=1;
    palette,"gray.gp"; 
  }
  
  (*ao_disp._winits)(4) = 1;
}


func control_win_init(pid5)
{
  extern ao_disp;
  
  if (numberof(*ao_disp._xid)==0) ao_disp._xid = &[0,0,0,0,0,0];
  (*ao_disp._xid)(5) = pid5;
  
  if (catch(0x08)) {
    (*ao_disp._winits)(5) = 1;
  }

  if (!(*ao_disp._winits)(5)) {
    stylename = yoga_ao_top+"/gs/yoga_ao.gs"; 
    window,(*ao_disp._wins)(5),dpi=27,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid5;
    limits,square=1;
    palette,"gray.gp"; 
  }

  (*ao_disp._winits)(5) = 1;

}
func target_win_init(pid6)
{
  extern ao_disp;
  
  if (numberof(*ao_disp._xid)==0) ao_disp._xid = &[0,0,0,0,0,0];
  (*ao_disp._xid)(6) = pid6;
  
  if (catch(0x08)) {
    (*ao_disp._winits)(6) = 1;
  }
  
  if (!(*ao_disp._winits)(6)) {
    stylename = yoga_ao_top+"/gs/yoga_ao.gs"; 
    window,(*ao_disp._wins)(6),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid6;
    limits,square=1;
    palette,"gray.gp"; 
  }

  (*ao_disp._winits)(6) = 1;
}

func reset_dm(ndm)
{
  if (g_dm == []) return;
  yoga_resetdm,g_dm,y_dm(ndm).type,y_dm(ndm).alt;
}

func start_ao_loop
{
  extern aoiter,aoloop, aotimer;

  aoiter = 0;
  aotimer = array(0., 100);
  aoloop = 1;
  animate,1;
  ao_loop;
}

func ao_loop(one)
{
  extern aoloop,aodisp_type,aodisp_numn, aotimer;
  extern g_atmos;
  extern y_atmos;
  extern time_move;

  if (aoloop) {
    if (g_atmos == []) return;

    mytime = tic();

    move_atmos,g_atmos;
    
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
        //sensors_trace,g_wfs,i-1,"all",g_atmos,g_dm;
	
        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        if ((!y_wfs(i).openloop) && (g_dm != [])) {
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
        }
	
        
        sensors_compimg,g_wfs,i-1;      
      }      
    }
    
    if ((y_wfs != []) && (g_wfs != [])) {
      // do centroiding
      if ((y_rtc != []) && (g_rtc != [])) {
        rtc_docentroids,g_rtc,0;
        // compute command and apply
        rtc_docontrol,g_rtc,0;
        if (g_dm != []) rtc_applycontrol,g_rtc,0,g_dm;
      }
    }
      
    if (enable_display == 1)
      update_main,aodisp_type,aodisp_num;    
/*  
    if (wfsiter < 1) time_move = tac(mytime);
    else {
      time_move += 0.01*tac(mytime);
      time_move /= 1.01;
    }
*/
    time_move = tac(mytime);
    aoiter ++;
    if(sum(aotimer) == 0) aotimer=array(time_move,1000);
    else                  aotimer(aoiter)=time_move; 
    aoiter = aoiter%1000;

    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('progressbar_wfs').set_fraction(%f)",float((aoiter%1000)/1000.));
    //mtext = swrite(format="Framerate : %.2f",1./time_move);
    strehltmp = target_getstrehl(g_target,0);
    mtext = swrite(format="Framerate : %.2f L.E. Strehl : %.2f, S.E. Strehl %.2f",1./aotimer(avg),strehltmp(2),strehltmp(1));
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('progressbar_wfs').set_text('%s')",mtext);
    
    if (!one) after,0.001,ao_loop;
  } else return;
}

func init_all(zenith,teldiam,cobs,r0,freq)
{
  extern y_atmos,y_wfs;
  extern g_wfs;
  extern imlp,strehllp,strehlsp,airy,sairy,niterok;

  if ((y_wfs == []) || (y_atmos == [])) return;
      
  //pyk,swrite(format=wfs_disp._cmd+"y_win_clear(%d)",1);
  
  pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('initall').set_sensitive(%d)",0);
  
  y_atmos.r0 = r0;
  if (y_tel == []) y_tel   = tel_struct();
  if (y_loop == []) y_loop = loop_struct();
  
  y_tel.diam         = teldiam;
  y_tel.cobs         = cobs;
  y_loop.ittime      = 1.0/freq;

  wfs_init;

  atmos_init;
  //init_gatmos,y_geom.pupdiam,zenith,teldiam,cobs,r0,freq;

  dm_init;
 
  target_init;

  rtc_init;
  
  update_main_display1,["atmos","wfs","image","slopes","centroids","dm","targetp","target"];
  
  write,"... Done with inits !";
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
  
  pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('initall').set_sensitive(%d)",1);
  
  imlp = 0.; // initializing average image
  strehllp = strehlsp = [];
  airy = roll(abs(fft(*y_geom._ipupil*exp(*y_geom._ipupil*1i*0.)))^2)/numberof(*y_geom._ipupil);
  sairy = max(airy);
  niterok=0;
  //y_target = [];
}
/*
 __  __       _                                _ 
|  \/  | __ _(_)_ __    _ __   __ _ _ __   ___| |
| |\/| |/ _` | | '_ \  | '_ \ / _` | '_ \ / _ \ |
| |  | | (_| | | | | | | |_) | (_| | | | |  __/ |
|_|  |_|\__,_|_|_| |_| | .__/ \__,_|_| |_|\___|_|
                       |_|                       
*/
func update_main_display1(type)
{
  pyk,swrite(format=wfs_disp._cmd+"y_win_clear(%d)",1);    
  //    
  //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').clear()");
  
  if (anyof(type == "atmos")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",0,
               "Phase - Atmos");
  }
  
  if (anyof(type == "wfs")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",1,
               "Phase - WFS");
  }
  
  if (anyof(type == "image")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",2,
               "Image - WFS");
  }
  
  if (anyof(type == "slopes")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",3,
               "Slopes - WFS");
  }
  
  if (anyof(type == "centroids")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",4,
               "Centroids - WFS");
  }
  if (anyof(type == "dm")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",5,
               "Phase - DM");
  }
  if (anyof(type == "targetp")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",5,
               "Phase - Target");
  }
  if (anyof(type == "target")) {
    pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",6,
               "Image - Target");
  }
  
  pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_type').set_active(%d)",0);
}

func update_main_display2(type)
{
  extern y_atmos;

  if (type == "Phase - Atmos") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Layer # %d",y_atmos.nscreens-cc+1));  
  }
  
  if (type == "Phase - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  if (type == "Image - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  
  if (type == "Slopes - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
     //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  if (type == "Centroids - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
    //pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  if (type == "Phase - DM") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
    //pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_dm);cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="DM # %d",numberof(y_dm)-cc+1));  
  }
  if (type == "Phase - Target") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
    //pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));  
  }
  if (type == "Image - Target") {
    pyk,swrite(format=wfs_disp._cmd+"y_winnum_clear(%d)",1);    
    //pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));  
  }  
  pyk,swrite(format=ao_disp._cmd+"glade.get_widget('winselect_number').set_active(%d)",0);
}

func update_main(type,nlayer)
{
  extern y_atmos,g_atmos,ao_disp;
  if (nlayer < 0) return;
  
  if (type == "Phase - Atmos") {
    if (nlayer > numberof(*y_atmos.alt)) return;
    mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(nlayer+1));
    window,(*ao_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase on layer # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*ao_disp._defaultdpi/200.);
  }
  
  if (type == "Phase - WFS") {
    if (nlayer > numberof(y_wfs)) return;
    mscreen = sensors_getdata(g_wfs,nlayer,"phase");

    //Pupil on top of phase screen
    if(y_wfs(nlayer+1).atmos_seen == 0){
      pup=float(make_pupil(y_geom.pupdiam,y_geom.pupdiam,type_ap=y_tel.type_ap,angle=y_tel.pupangle,spiders_type=y_tel.spiders_type,t_spiders=y_tel.t_spiders,nbr_miss_seg=y_tel.nbrmissing,std_ref_err=y_tel.referr,xc=cent,yc=cent,real=,cobs=));
      pupcustom = pad_array(pup,y_geom._n);
    }
    else
      pupcustom = *y_geom._mpupil; // modif: lecture pupille
    mscreen = mscreen*pupcustom; // modif: multiplication mscreen par pupille

    //mscreen = sensors_getdata(g_wfs,nlayer,"phasetele");
    window,(*ao_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase seen from WFS # %d",nlayer+1);

    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*ao_disp._defaultdpi/200.);
    rho=indgen(1000)/1000.*2*pi;
    plg,y_geom.pupdiam/2*cos(rho)+y_geom._n/2.,y_geom.pupdiam/2*sin(rho)+y_geom._n/2.,color="red",marks=0,width=3;
  }
  
  if (type == "Image - WFS") {
    if (nlayer > numberof(y_wfs)) return;
    
    if (y_wfs(nlayer+1).type == "sh") mimg = sensors_getdata(g_wfs,nlayer,"imgtele");
    else mimg = sensors_getdata(g_wfs,nlayer,"bincube");
    
    window,(*ao_disp._wins)(1);fma;
    if (y_wfs(nlayer+1).type == "sh") pli,mimg;
    else {
      pli,get_pyrimg(mimg);
    }
    myxtitle = swrite(format="image on WFS # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*ao_disp._defaultdpi/200.);
  }
  if (type == "Slopes - WFS") {
    if (nlayer > numberof(y_wfs)) return;
    //sensors_trace,g_wfs,nlayer,"atmos",g_atmos;
    slopes_geom,g_wfs,nlayer,0;
    slps = sensors_getslopes(g_wfs,nlayer);
    window,(*ao_disp._wins)(1);fma;
    display_slopes,slps,nlayer+1,"Phase Difference";
  }
  if (type == "Centroids - WFS") {
    //slps = sensors_getslopes(g_wfs,nlayer);
    slps = rtc_getcentroids(g_rtc,0);
    window,(*ao_disp._wins)(1);fma;
    display_slopes,slps,nlayer+1,"Centroids";
  }
  if (type == "Phase - DM") {
    mscreen = yoga_getdm(g_dm,y_dm(nlayer+1).type,y_dm(nlayer+1).alt);
    window,(*ao_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase on DM # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*ao_disp._defaultdpi/200.);
    nx = dimsof(mscreen)(2);
    rho=indgen(1000)/1000.*2*pi;
    plg,y_geom.pupdiam/2*cos(rho)+nx/2.,y_geom.pupdiam/2*sin(rho)+nx/2.,color="red",marks=0,width=3;
  }
  if (type == "Phase - Target") {
    //pupcustom = *y_geom._apodizer; // modif: lecture pupille
    pupcustom = *y_geom._spupil;
      //target_atmostrace,g_target,nlayer,g_atmos;
    if (g_dm != []) {
      //target_dmtrace,g_target,nlayer,g_dm;
      //mscreen = target_getphasetele(g_target,nlayer);
      mscreen = target_getphase(g_target,nlayer);
      mscreen = mscreen*pupcustom; // modif: multiplication mscreen par pupille
    } else mscreen = target_getphase(g_target,nlayer);
    window,(*ao_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase seen from target # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*ao_disp._defaultdpi/200.);
    nx = dimsof(mscreen)(2);
    rho=indgen(1000)/1000.*2*pi;
    plg,y_geom.pupdiam/2*cos(rho)+nx/2.,y_geom.pupdiam/2*sin(rho)+nx/2.,color="red",marks=0,width=3;
  }
  if (type == "Image - Target") {
    mimg = target_getimage(g_target,nlayer,"se");
    window,(*ao_disp._wins)(1);fma;
    Di=dimsof(mimg);
    mimg=roll(mimg);
    //OWA=2*y_dm.nact(1);
    OWA = Di(2) > 128 ? 64 : Di(2);
    mimg=mimg(Di(2)/2-1-OWA:Di(2)/2+OWA,Di(2)/2-1-OWA:Di(2)/2+OWA);
    mimg=mimg*(mimg > max(mimg)/100000)+max(mimg)/100000*(mimg < max(mimg)/100000);
    //pli,log10(mimg);    
    pli,mimg;
    myxtitle = swrite(format="image for target # %d", nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*ao_disp._defaultdpi/200.);
  }

}

func ao_unzoom
{
  animate,0;
  for (i=1;i<=numberof(*ao_disp._wins);i++) {
    if ((*ao_disp._winits)(i)) {
      window,(*ao_disp._wins)(i);
      unzoom;
    }
  }
}

//////////////////////////////////////////////////////////
//              **********************                 //
////////         END OF ROUTINES DEFINITION      ////////
//              **********************                 //
//////////////////////////////////////////////////////////

// start standalone version if called from shell  
//pyk_debug=1;

arg_wfs = get_argv();
ao_disp = display_struct();

ao_disp._cmd = "wao.";

if (anyof(strmatch(arg_wfs,"widget_ao.i")) || get_env("EMACS")=="t" ) {
  ao_disp._cmd = "";
  python_exec = yoga_ao_top+"/widgets/widget_ao.py";
  pyk_cmd=[python_exec];
  if (!_pyk_proc) _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  write,"widget_wfs  ready";
  write,"standalone version";
}

ao_disp._wins            = &[13,14,15,16,17,18];
ao_disp._winits          = &[0,0,0,0,0,0];
ao_disp._defaultdpi      = 130;    
ao_disp._ncolors         = 200;
ao_disp._lut             = 0;     // default LUT index [0-41]
ao_disp._xytitles_adjust = &[0.012,0.019]; // X and Y notch axis titles in main area
ao_disp._invertlut       = 0;
ao_disp._gui_realized    = 0;
ao_disp._zoom            = 1;

pldefault,opaque=1,marks = 1;

aodisp_type = "";
aodisp_num  = 0;
aoloop      = 0;
aoiter      = 0;
aotimer     = array(0., 100);


enable_display = 1;

atmos_disp = ao_disp;
sky_disp = ao_disp;
dm_disp = ao_disp;
wfs_disp = ao_disp;

yoga_parfilename = "";

nmax_device = yoga_getnDevice();
pyk,swrite(format=ao_disp._cmd+"glade.get_widget('device').set_range(%d,%d)",0,nmax_device-1);  
//pyk, ao_disp._cmd+"on_load_defaults_clicked()";
