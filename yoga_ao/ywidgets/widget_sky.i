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
require,yoga_ao_top+"/yoga_ao.i";

//////////////////////////////////////////////////////////
// A basic set of display functions
//////////////////////////////////////////////////////////

func sky_win_init(pid1,pid2,pid3)
{
  extern sky_disp;

  sky_disp._xid=&[pid1,pid2,pid3];

  if (catch(0x08)) {
    sky_disp._gui_realized = 1;
  }
  
  if (!sky_disp._gui_realized) {
    stylename = yoga_ao_top+"/gs/yoga_ao.gs";
    window,(*sky_disp._wins)(1),dpi=sky_disp._defaultdpi,width=0,height=0,
      xpos=-2,ypos=-2,style=stylename,parent=pid1;
    limits,square=1;
    //palette,"gray.gp"; // need this if loadct is used!?
    /*
    widget_set_lut,sky_disp._lut,(*sky_disp._wins)(1),sky_disp._lut,
                   sky_disp._ncolors,sky_disp._itt,sky_disp._log_itt_dex,
                   sky_disp._rlut,sky_disp._glut,sky_disp._blut;
    */

    window,(*sky_disp._wins)(2),dpi=35,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid2;
    limits,square=1;
    palette,"gray.gp"; 
   
    window,(*sky_disp._wins)(3),dpi=37,width=0,height=0,
      xpos=-2,ypos=-25,style=stylename,parent=pid3;
    limits,square=1;
    palette,"gray.gp"; 
   
  }

  sky_disp._gui_realized = 1;
}


func init_all(pupdiam,zenith,teldiam,cobs,r0,freq)
{
  extern y_atmos,y_target;
  extern g_target;
  
  init_gatmos,pupdiam,zenith,teldiam,cobs,r0/100.,freq;

  target_init;
  
  update_main_display1,["atmos","target","image"];
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
  if (anyof(type == "atmos")) {
    for (cc=1;cc<=4;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').remove_text(%d)",0);
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",0,
               "Phase - Atmos");
  }
  if (anyof(type == "target")) {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",1,
               "Phase - Target");
  }
  if (anyof(type == "image")) {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",2,
               "Image - Target");
  }
  
  pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').set_active(%d)",0);
}

func update_main_display2(type)
{
  extern y_atmos;

  if (type == "Phase - Atmos") {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').clear()");
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Layer # %d",y_atmos.nscreens-cc+1));  
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
  pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').set_active(%d)",0);
}

func update_main(type,nlayer)
{
  extern y_atmos,g_atmos,sky_disp;
  if (nlayer < 0) return;
  
  if (type == "Phase - Atmos") {
    mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(nlayer+1));
    window,(*sky_disp._wins)(1);fma;limits;
    pli,mscreen;
  }
  if (type == "Phase - Target") {
    target_atmostrace,g_target,nlayer,g_atmos;
    mscreen = target_getphase(g_target,nlayer);
    window,(*sky_disp._wins)(1);fma;limits;
    pli,mscreen;
  }
  if (type == "Image - Target") {
    mscreen = target_getimage(g_target,nlayer,"se");
    window,(*sky_disp._wins)(1);fma;limits;
    pli,eclat(mscreen);
  }
}

func start_sky_loop
{
  extern skyiter,skyloop;

  skyiter = 0;
  skyloop = 1;
  animate,1;
  sky_loop;
}

func sky_loop(one)
{
  extern skyloop,skydisp_type,skydisp_num;
  extern g_atmos,g_target;
  extern y_atmos;
  extern time_move;

  if (skyloop) {
    if (g_atmos == []) return;

    mytime = tic();

    //if (g_target == []) extrude_tscreen,g_atmos,(*y_atmos.alt)(1);
    if (g_target == []) move_atmos,g_atmos;
    else move_sky,g_atmos,g_target;

    
    if (skydisp_type != [])
      update_main,skydisp_type,skydisp_num;    
    
    if (skyiter < 1) time_move = tac(mytime);
    else {
      time_move += 0.01*tac(mytime);
      time_move /= 1.01;
    }

    skyiter ++;

    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('progressbar_sky').set_fraction(%f)",float((skyiter%1000)/1000.));
    mtext = swrite(format="Framerate : %.2f",1./time_move);
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('progressbar_sky').set_text('%s')",mtext);
    
    if (!one) after,0.001,sky_loop;
  } else return;
}


//////////////////////////////////////////////////////////
//              **********************                 //
////////         END OF ROUTINES DEFINITION      ////////
//              **********************                 //
//////////////////////////////////////////////////////////

// start standalone version if called from shell  
//pyk_debug=1;

arg_sky = get_argv();
sky_disp = display_struct();

sky_disp._cmd = "sky.";

if (anyof(strmatch(arg_sky,"widget_sky.i")) || get_env("EMACS")=="t" ) {
  sky_disp._cmd = "";
  python_exec = yoga_ao_top+"/widgets/widget_sky.py";
  pyk_cmd=[python_exec];
  if (!_pyk_proc) _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  write,"widget_sky  ready";
  write,"standalone version";
}

sky_disp._wins            = &[10,11,12];
sky_disp._defaultdpi      = 85;    // change size of spydr graphic area
sky_disp._ncolors         = 200;
sky_disp._lut             = 0;     // default LUT index [0-41]
sky_disp._xytitles_adjust = &[0.012,0.019]; // X and Y notch axis titles in main area
sky_disp._invertlut       = 0;
pldefault,opaque=1,marks  =1;
sky_disp._gui_realized    = 0;
sky_disp._zoom            = 1;

skydisp_type = "";
skydisp_num  = 0;
skyloop      = 0;
skyiter      = 0;

atmos_disp = sky_disp;
