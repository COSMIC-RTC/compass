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

func update_control_plot(type,ncontrol)
{
  extern y_rtc;
  window,(*ao_disp._wins)(5);fma;logxy,0,0;
  if (y_rtc == []) return;
  else {
    if(type == 0) {// Imat eigenvalues
      if((*y_rtc.controllers)(ncontrol).type == "mv") eigenv = rtc_getDeigenvals(g_rtc,0);
      else eigenv = controller_getdata(g_rtc,ncontrol,"eigenvals");
      maxcond = (*y_rtc.controllers)(ncontrol).maxcond;
      if (eigenv(1) < eigenv(0)) mfilt = where((eigenv/eigenv(-2)) < 1./maxcond);
      else mfilt = where(1./(eigenv/eigenv(3)) > maxcond);
      nfilt = numberof(mfilt);
      if ( (ao_disp!=[]) && (numberof(*ao_disp._winits) > 0)) {
        if ((*ao_disp._winits)(5)) {
          logxy,0,1;
          if (eigenv(1) < eigenv(0)) {
            plg, eigenv(::-1), marks=0;
            plmk, eigenv(::-1), msize = 0.3, marker=4;
          } else {
            plg, eigenv, marks=0;
            plmk, eigenv, msize = 0.3, marker=4;
          }
          x0 = numberof(eigenv) - nfilt + 0.5;
          pldj, x0 ,min(eigenv), x0, max(eigenv), color="red";
        }
      }
    } else if(type == 1) { //Cmm
      if((*y_rtc.controllers)(ncontrol).type == "mv"){
        maxcond = (*y_rtc.controllers)(ncontrol).maxcond;
        Cphim = mat_cphim_gpu(0);
        Cmm = rtc_getCmm(g_rtc,0);
        pli,Cmm;
        rtc_buildcmatmv,g_rtc,0,maxcond;
        rtc_filtercmatmv,g_rtc,0,maxcond;
      }

    } else if(type == 2) { //Cmm inverse
      if((*y_rtc.controllers)(ncontrol).type == "mv"){
        Cmm = rtc_getCmm(g_rtc,0);
        pli,Cmm
      }
    } else if(type == 3) { //Cmm eigenvalues
      if((*y_rtc.controllers)(ncontrol).type == "mv"){
        eigenv = rtc_getCmmeigenvals(g_rtc,0);
        maxcond = (*y_rtc.controllers)(ncontrol).maxcond;
        if (eigenv(1) < eigenv(0)) mfilt = where((eigenv/eigenv(-2)) < 1./maxcond);
        else mfilt = where(1./(eigenv/eigenv(3)) > maxcond);
        nfilt = numberof(mfilt);
        if ( (ao_disp!=[]) && (numberof(*ao_disp._winits) > 0)) {
          if ((*ao_disp._winits)(5)) {
            logxy,0,1;
            if (eigenv(1) < eigenv(0)) {
              plg, eigenv(::-1), marks=0;
              plmk, eigenv(::-1), msize = 0.3, marker=4;
            } else {
              plg, eigenv, marks=0;
              plmk, eigenv, msize = 0.3, marker=4;
            }
            x0 = numberof(eigenv) - nfilt + 0.5;
            pldj, x0 ,min(eigenv), x0, max(eigenv), color="red";
          }
        }
      }
    } else if(type == 4) { //Cphim
      Cphim = rtc_getCphim(g_rtc,0);
      pli,Cphim;
    } else if(type == 5) { //cmat
      cmat = rtc_getcmat(g_rtc,0);
      pli,cmat;
    }
  }
}

func cmat_update(ncontrol,maxcond)
{
  extern y_rtc;
  (*y_rtc.controllers)(ncontrol).maxcond = maxcond;

  //error;
  if((*y_rtc.controllers)(ncontrol).type == "ls"){
    eigenv = controller_getdata(g_rtc,ncontrol,"eigenvals");
  
    if (eigenv(1) < eigenv(0)) mfilt = where((eigenv/eigenv(-2)) < 1./maxcond);
    else mfilt = where(1./(eigenv/eigenv(3)) > maxcond);
    //nfilt = numberof(mfilt)+2;
    nfilt = numberof(mfilt);

    write,format="nb modes filtered : %d",nfilt;
    type = 0;
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('control_plot').set_active(%d)",type);
    update_control_plot(type,ncontrol);
  
    write,"building cmat";
    tic;
    rtc_buildcmat,g_rtc,ncontrol,nfilt;
    write,format="cmat time %f\n",tac();

    cmat = rtc_getcmat(g_rtc,ncontrol);
  
    controllers = *y_rtc.controllers;
    controllers(ncontrol).cmat = &float(cmat);
    y_rtc.controllers = &controllers;
  }
  else if((*y_rtc.controllers)(ncontrol).type == "mv"){
    Cphim = mat_cphim_gpu(ncontrol-1);
    write,"Building cmat...";
    rtc_buildcmatmv,g_rtc,ncontrol-1,maxcond;
    rtc_filtercmatmv,g_rtc,ncontrol-1,maxcond;
    write,"cmat done"
    type = 3;
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('control_plot').set_active(%d)",type);
    update_control_plot(type,ncontrol);
  }
}

