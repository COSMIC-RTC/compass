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


func load_parfile(parfile,filename)
{
  extern yoga_parfilename;
  extern y_wfs,y_atmos,y_rtc,y_loop,y_tel,y_geom,y_target;
  
  yoga_parfilename = filename;

  y_wfs = y_atmos = y_rtc = y_loop = y_tel = y_geom = y_target = [];
  read_parfile,parfile;
  //y_atmos.nscreens;
  if (y_tel != []) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('teldiam').set_value(%f)",y_tel.diam);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('cobs').set_value(%f)",y_tel.cobs);
  }
  if (y_atmos != []) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('nlayers').set_value(%d)",y_atmos.nscreens);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('r0').set_value(%f)",y_atmos.r0);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('layer_select').clear()");
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('layer_select').insert_text(%d,'%s')",0,
                 swrite(format="Layer # %d",numberof(y_wfs)-cc+1));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('layer_select').set_active(%d)",0);
  }
  if (y_wfs != []) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('nwfs').set_value(%d)",numberof(y_wfs));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",0);
    if (y_rtc != []) pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",0);
    else pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('default_centro').set_active(%d)",0);
  }
  if (y_target != []) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('ntargets').set_value(%d)",y_target.ntargets);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('target_select').clear()");
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('target_select').set_active(%d)",0);
  }
}

func wfs_win_init(pid1,pid2,pid3,pid4,pid5,pid6)
{
  extern wfs_disp;

  wfs_disp._xid=&[pid1,pid2,pid3,pid4,pid5,pid6];

  if (catch(0x08)) {
    wfs_disp._gui_realized = 1;
  }
  
  if (!wfs_disp._gui_realized) {
    stylename = yoga_ao_top+"/gs/yoga_noaxis.gs";
    window,(*wfs_disp._wins)(1),dpi=wfs_disp._defaultdpi,width=0,height=0,
      xpos=-4,ypos=-4,style=stylename,parent=pid1;
    limits,square=1;

    stylename = yoga_ao_top+"/gs/yoga_ao.gs";
    window,(*wfs_disp._wins)(2),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid2;
    limits,square=1;
    palette,"gray.gp"; 
   
    window,(*wfs_disp._wins)(3),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid3;
    limits,square=1;
    palette,"gray.gp"; 
   
    window,(*wfs_disp._wins)(4),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid4;
    limits,square=1;
    palette,"gray.gp"; 
   
    window,(*wfs_disp._wins)(5),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid5;
    limits,square=1;
    palette,"gray.gp"; 
   
    window,(*wfs_disp._wins)(6),dpi=25,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid6;
    limits,square=1;
    palette,"gray.gp"; 
   
  }

  wfs_disp._gui_realized = 1;
}



func start_wfs_loop
{
  extern wfsiter,wfsloop, wfstimer;

  wfsiter = 0;
  wfstimer = array(0., 100);
  wfsloop = 1;
  animate,1;
  wfs_loop;
}

func wfs_loop(one)
{
  extern wfsloop,wfsdisp_type,wfsdisp_numn, wfstimer;
  extern g_atmos;
  extern y_atmos;
  extern time_move;

  if (wfsloop) {
    if (g_atmos == []) return;

    mytime = tic();

    if (g_target != []) move_sky,g_atmos,g_target;
    else move_atmos,g_atmos;
    
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (i=1;i<=numberof(y_wfs);i++) {

        sensors_trace,g_wfs,i-1,"atmos",g_atmos;
        if ((!y_wfs(i).openloop) && (g_dm != [])) {
          sensors_trace,g_wfs,i-1,"dm",g_dm,0;
        }
        
        if ((enable_display == 1) && (wfsdisp_type == "Image - WFS") && (i == wfsdisp_num + 1)) {
            sensors_compimg_tele,g_wfs,i-1;
        } else sensors_compimg,g_wfs,i-1;
        
        //sensors_compimg,g_wfs,i-1;
      }
      
      // do centroiding
      if ((y_rtc != []) && (g_rtc != [])) {
        rtc_docentroids,g_rtc,g_wfs,0;
        // compute command and apply
        if (g_dm != []) rtc_docontrol,g_rtc,0,g_dm;
      }
    }
    
    if (enable_display == 1)
      update_main,wfsdisp_type,wfsdisp_num;    
/*  
    if (wfsiter < 1) time_move = tac(mytime);
    else {
      time_move += 0.01*tac(mytime);
      time_move /= 1.01;
    }
*/
    time_move = tac(mytime);
    wfsiter ++;
	if(sum(wfstimer) == 0) wfstimer=array(time_move,100);
	else                  wfstimer(wfsiter)=time_move; 
	wfsiter = wfsiter%100;

    //pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('progressbar_wfs').set_fraction(%f)",float((wfsiter%1000)/1000.));
    //mtext = swrite(format="Framerate : %.2f",1./time_move);
    mtext = swrite(format="Framerate : %.2f",1./wfstimer(avg));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('progressbar_wfs').set_text('%s')",mtext);
    
    if (!one) after,0.001,wfs_loop;
  } else return;
}

func init_all(zenith,teldiam,cobs,r0,freq)
{
  extern y_atmos,y_wfs;
  extern g_wfs;

  if ((y_wfs == []) || (y_atmos == [])) return;
      
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
  pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').clear()");
  
  if (anyof(type == "atmos")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",0,
               "Phase - Atmos");
  }
  
  if (anyof(type == "wfs")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",1,
               "Phase - WFS");
  }
  
  if (anyof(type == "image")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",2,
               "Image - WFS");
  }
  
  if (anyof(type == "slopes")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",3,
               "Slopes - WFS");
  }
  
  if (anyof(type == "centroids")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",4,
               "Centroids - WFS");
  }
  if (anyof(type == "dm")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",5,
               "Phase - DM");
  }
  if (anyof(type == "targetp")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",5,
               "Phase - Target");
  }
  if (anyof(type == "target")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",6,
               "Image - Target");
  }
  
  pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').set_active(%d)",0);
}

func update_main_display2(type)
{
  extern y_atmos;

  if (type == "Phase - Atmos") {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Layer # %d",y_atmos.nscreens-cc+1));  
  }
  
  if (type == "Phase - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  if (type == "Image - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  
  if (type == "Slopes - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  if (type == "Centroids - WFS") {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  if (type == "Phase - DM") {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=numberof(y_dm);cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="DM # %d",y_target.ntargets-cc+1));  
  }
  if (type == "Phase - Target") {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));  
  }
  if (type == "Image - Target") {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));  
  }  
  pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').set_active(%d)",0);
}

func update_main(type,nlayer)
{
  extern y_atmos,g_atmos,wfs_disp;
  if (nlayer < 0) return;
  
  if (type == "Phase - Atmos") {
    if (nlayer > numberof(*y_atmos.alt)) return;
    mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(nlayer+1));
    window,(*wfs_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase on layer # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*wfs_disp._defaultdpi/200.);
  }
  
  if (type == "Phase - WFS") {
    if (nlayer > numberof(y_wfs)) return;
    //mscreen = sensors_getdata(g_wfs,nlayer,"phase");
    mscreen = sensors_getdata(g_wfs,nlayer,"phasetele");
    window,(*wfs_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase seen from WFS # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*wfs_disp._defaultdpi/200.);
  }
  
  if (type == "Image - WFS") {
    if (nlayer > numberof(y_wfs)) return;
    //mimg = sensors_getimg(g_wfs,nlayer);
    mimg = sensors_getdata(g_wfs,nlayer,"imgtele");
    window,(*wfs_disp._wins)(1);fma;
    pli,mimg;
    myxtitle = swrite(format="image on WFS # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*wfs_disp._defaultdpi/200.);
  }
  if (type == "Slopes - WFS") {
    if (nlayer > numberof(y_wfs)) return;
    //sensors_trace,g_wfs,nlayer,"atmos",g_atmos;
    slopes_geom,g_wfs,nlayer,0;
    slps = sensors_getslopes(g_wfs,nlayer);
    window,(*wfs_disp._wins)(1);fma;
    display_slopes,slps,nlayer+1,"Phase Difference";
  }
  if (type == "Centroids - WFS") {
    slps = sensors_getslopes(g_wfs,nlayer);
    window,(*wfs_disp._wins)(1);fma;
    display_slopes,slps,nlayer+1,"Centroids";
  }
  if (type == "Phase - DM") {
    mscreen = yoga_getdm(g_dm,y_dm(nlayer+1).type,y_dm(nlayer+1).alt);
    window,(*wfs_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase on DM # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*wfs_disp._defaultdpi/200.);
  }
  if (type == "Phase - Target") {
    target_atmostrace,g_target,nlayer,g_atmos;
    if (g_dm != []) {
      target_dmtrace,g_target,nlayer,g_dm;
      mscreen = target_getphasetele(g_target,nlayer);
    } else mscreen = target_getphase(g_target,nlayer);
    window,(*wfs_disp._wins)(1);fma;
    pli,mscreen;
    myxtitle = swrite(format="phase seen from target # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*wfs_disp._defaultdpi/200.);
  }
  if (type == "Image - Target") {
    if (g_dm == []) mimg = target_getimage(g_target,nlayer,"se");
    else mimg = target_getimage(g_target,nlayer,"se");
    
    window,(*wfs_disp._wins)(1);fma;
    pli,roll(mimg);
    myxtitle = swrite(format="image for target # %d",nlayer+1);
    port= viewport();
    plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
      font="helvetica", justify="CB", height=long(pltitle_height*wfs_disp._defaultdpi/200.);
  }

}

func wfs_unzoom
{
  animate,0;
  for (i=1;i<=numberof(*wfs_disp._wins);i++) {
    window,(*wfs_disp._wins)(i);
    unzoom;
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
wfs_disp = display_struct();

wfs_disp._cmd = "wfs.";

if (anyof(strmatch(arg_wfs,"widget_wfs.i")) || get_env("EMACS")=="t" ) {
  wfs_disp._cmd = "";
  python_exec = yoga_ao_top+"/widgets/widget_wfs.py";
  pyk_cmd=[python_exec];
  if (!_pyk_proc) _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  write,"widget_wfs  ready";
  write,"standalone version";
}

wfs_disp._wins            = &[13,14,15,16,17,18];
wfs_disp._defaultdpi      = 130;    
wfs_disp._ncolors         = 200;
wfs_disp._lut             = 0;     // default LUT index [0-41]
wfs_disp._xytitles_adjust = &[0.012,0.019]; // X and Y notch axis titles in main area
wfs_disp._invertlut       = 0;
pldefault,opaque=1,marks  = 1;
wfs_disp._gui_realized    = 0;
wfs_disp._zoom            = 1;

wfsdisp_type = "";
wfsdisp_num  = 0;
wfsloop      = 0;
wfsiter      = 0;
wfstimer     = array(0., 100);


enable_display = 1;

atmos_disp = wfs_disp;
sky_disp = wfs_disp;
dm_disp = wfs_disp;

yoga_parfilename = "";

nmax_device = yoga_getnDevice();
pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('device').set_range(%d,%d)",0,nmax_device-1);  
pyk, wfs_disp._cmd+"on_load_defaults_clicked()";
