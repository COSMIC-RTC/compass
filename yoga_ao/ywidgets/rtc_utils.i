func update_centro_prop(ncentro)
{
  extern wfs_disp;
  
  if (y_rtc == []) {
    pyk,swrite(format=wfs_disp._cmd+"y_init_centro(%d)",1);
  } else {
    if (y_rtc.nwfs == 0) pyk,swrite(format=wfs_disp._cmd+"y_init_centro(%d)",1);
    else {
      if (ncentro > numberof(y_wfs)) error,"not a valid centroider";
      else {
        pyk,swrite(format=wfs_disp._cmd+"y_update_centro_gui('%s', %d, %f, '%s', %f)",
                   (*y_rtc.centroiders)(ncentro).type,(*y_rtc.centroiders)(ncentro).nmax,
                   (*y_rtc.centroiders)(ncentro).thresh,(*y_rtc.centroiders)(ncentro).type_fct,
                   (*y_rtc.centroiders)(ncentro).width);
      }
    }
  }
}

func init_centro_prop(ncentro,type,nmax,thresh,type_fct,width)
{
  extern y_rtc;

  if (y_rtc == []) return;
  else {
    (*y_rtc.centroiders)(ncentro).type     = type;
    (*y_rtc.centroiders)(ncentro).nmax     = nmax;
    (*y_rtc.centroiders)(ncentro).thresh   = thresh;
    (*y_rtc.centroiders)(ncentro).width    = width;
    (*y_rtc.centroiders)(ncentro).type_fct = type_fct;
  }
}


func load_default_centro(type,type_fct)
{
  extern y_rtc;
  
  if (y_wfs == []) {
    write,"warning : no wfs defined, centroider was not initizalized";
    return;
  }
  
  ncentroiders   = numberof(y_wfs);
  y_centroiders  = array(centroider_struct(),ncentroiders);
  for (i=1;i<=ncentroiders;i++) {
    if (y_wfs(i).type == "pyr") type = "pyr";
    y_centroiders(i).nwfs      = i;
    y_centroiders(i).type   = type;
    if (type == "tcog") y_centroiders(i).thresh = 100.;
    else y_centroiders(i).thresh = 0;
    if (type == "bpcog") y_centroiders(i).nmax   = 10;
    else y_centroiders(i).nmax = 0;
    if (type == "wcog") {
      y_centroiders(i).type_fct = type_fct;
      if (type_fct == "gauss") {
        y_centroiders(i).width   = 2.;
      }
    } else y_centroiders(i).nmax = 0;
  }
  
  if (y_rtc == []) y_rtc    = rtc_struct();

  y_rtc.nwfs = ncentroiders;
  y_rtc.centroiders = &(y_centroiders);

  pyk,swrite(format=wfs_disp._cmd+"y_centro_clear(%d)",1);    
  //pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').clear()");
  for (cc=1;cc<=numberof(y_wfs);cc++)
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').insert_text(%d,'%s')",0,swrite(format="centro # %d",numberof(y_wfs)-cc+1));
  pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",0);
  if (type == "wcog") {
    if (type_fct == "gauss") {
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('type_func').set_active(%d)",0);
    }
  }
}



func update_control_prop(void)
{
  extern wfs_disp;
  
  if (y_rtc == []) {
    pyk,swrite(format=wfs_disp._cmd+"y_init_control(%d)",1);
  } else {
    pyk,swrite(format=wfs_disp._cmd+"y_update_control_gui('%s', %f, %f, %f)",
               (*y_rtc.controllers)(1).type,(*y_rtc.controllers)(1).maxcond,
               (*y_rtc.controllers)(1).delay,(*y_rtc.controllers)(1).gain);
  }
}



func init_control_prop(type,maxcond,delay,gain)
{
  extern y_rtc;

  if (y_rtc == []) return;
  else {
    if (typeof(type) == "string") {
      if (type == "LS") type = "ls";
      if (type == "MV") type = "mv";
      if (type == "L&A") type = "la";
      if (type == "LQG") type = "lqg";
    } else {
      if (type == 1) type = "ls";
      if (type == 2) type = "mv";
      if (type == 3) type = "la";
      if (type == 4) type = "lqg";
    }
    (*y_rtc.controllers)(1).type    = type;
    (*y_rtc.controllers)(1).maxcond = maxcond;
    (*y_rtc.controllers)(1).delay   = delay;
    (*y_rtc.controllers)(1).gain    = gain;
    if (g_rtc != []) rtc_setgain,g_rtc,0,gain;
  }
}


func load_default_control(type)
{
  extern y_rtc,y_controllers;

  if (y_wfs == []) {
    write,"warning : no wfs defined, controller was not initizalized";
    return;
  }
  
  if (y_dm == []) {
    write,"warning : no dm defined, controller was not initizalized";
    return;
  }
  
  if (type == "LS") type = "ls";
  if (type == "MV") type = "mv";
  if (type == "L&A") type = "la";
  if (type == "LQG") type = "lqg";
  
  y_controllers = array(controller_struct(),1);
  y_controllers(1).type    = type;
  y_controllers(1).nwfs    = &(indgen(numberof(y_wfs)));
  y_controllers(1).ndm     = &(indgen(numberof(y_dm)));
  y_controllers(1).maxcond = 20.;
  y_controllers(1).delay   = 1;
  y_controllers(1).gain    = 0.4;

  if (y_rtc == []) y_rtc    = rtc_struct();
  y_rtc.controllers = &(y_controllers);

  update_control_prop;
}


func update_control(type,maxcond,delay,gain)
{
  extern y_rtc;

  if (y_rtc == []) return;
  else {
    if (typeof(type) == "string") {
      if (type == "LS") type = "ls";
      if (type == "MV") type = "mv";
      if (type == "L&A") type = "la";
      if (type == "LQG") type = "lqg";
    } else {
      if (type == 1) type = "ls";
      if (type == 2) type = "mv";
      if (type == 3) type = "la";
      if (type == 4) type = "lqg";
    }
    (*y_rtc.controllers)(1).type    = type;
    (*y_rtc.controllers)(1).maxcond = maxcond;
    (*y_rtc.controllers)(1).delay   = delay;
    (*y_rtc.controllers)(1).gain    = gain;
    if (g_rtc != []) {
      cmat_update,0,maxcond;
      rtc_setdelay,g_rtc,0,delay;
      rtc_setgain,g_rtc,0,gain;
    }
  }
}



