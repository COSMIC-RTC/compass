require,yoga_ao_top+"/ywidgets/widgets_utils.i";

/*
    _   _                       
   / \ | |_ _ __ ___   ___  ___ 
  / _ \| __| '_ ` _ \ / _ \/ __|
 / ___ \ |_| | | | | | (_) \__ \
/_/   \_\__|_| |_| |_|\___/|___/
                                
*/

func update_layer_prop(nlayer)
{
  extern atmos_disp,y_atmos;
  
  if (y_atmos == []) {
    pyk,swrite(format=atmos_disp._cmd+"y_init_atmos(%d)",1);
  } else {
    if (nlayer + 1 > y_atmos.nscreens) error,"not a valid layer";
    else {
      if (numberof((*y_atmos.dim_screens)) == 0) dims = 0;
      else dims = (*y_atmos.dim_screens)(nlayer+1);
      pyk,swrite(format=atmos_disp._cmd+"y_update_layer_gui(%d, %f, %f, %f, %f, '%s')",
                 long((*y_atmos.alt)(nlayer+1)),float((*y_atmos.frac)(nlayer+1)),
                 float((*y_atmos.L0)(nlayer+1)),float((*y_atmos.windspeed)(nlayer+1)),float((*y_atmos.winddir)(nlayer+1)),
                 swrite(format="%d",long(dims)));
      pyk,swrite(format=atmos_disp._cmd+"y_layers_plot_update(%d)",1);
    }
  }
}

func create_atmos(teldiam,pupdiam,zenith,nlayers,r0)
{
  extern y_atmos;
  
  if (y_geom != []) pupdiam = y_geom.pupdiam;
  
  y_atmos  = atmos_struct();
  
  y_atmos.r0           = r0;
  y_atmos.nscreens     = nlayers;
  y_atmos.frac         = &(array(1./nlayers,nlayers));
  y_atmos.L0           = &(array(100.0f,nlayers));
  y_atmos.alt          = &((indgen(nlayers)-1)*1000.);
  y_atmos.windspeed    = &(indgen(nlayers)*10.);
  y_atmos.winddir      = &(indgen(nlayers)*0.);
  y_atmos.pupixsize    = teldiam / pupdiam;

  if (y_target == []) y_atmos.dim_screens =&(array(pupdiam+4,nlayers));
  else {
    max_size = max(abs(*y_target.xpos,*y_target.ypos));
    patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
    patch_diam += patch_diam % 2;
    y_atmos.dim_screens = &patch_diam;
  }
  pyk,swrite(format=atmos_disp._cmd+"y_layer_clear(%d)",1);
  //pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').clear()");
  for (cc=1;cc<=nlayers;cc++)
    pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').insert_text(%d,'%s')",0,swrite(format="Layer # %d",nlayers-cc+1));
  
  update_layer_prop,0;
  pyk,swrite(format=atmos_disp._cmd+"y_layers_plot_update(%d)",1);
}

func update_nlayers(pupdiam,nlayers)
{
  extern y_atmos;
  
  if (y_atmos != []) {
    if (y_geom != []) pupdiam = y_geom.pupdiam;
    
    if (nlayers > y_atmos.nscreens) {
      y_atmos.nscreens     += 1;
      y_atmos.frac         = &(_((*y_atmos.frac),0.1));
      *y_atmos.frac       /= sum(*y_atmos.frac);
      y_atmos.L0           = &(_((*y_atmos.L0),100));
      y_atmos.alt          = &(_((*y_atmos.alt),10000));
      y_atmos.windspeed    = &(_((*y_atmos.windspeed),10));
      y_atmos.winddir      = &(_((*y_atmos.winddir),0));
      if (y_target == []) y_atmos.dim_screens =&(_((*y_atmos.dim_screens),pupdiam+4));
      else {
        max_size = max(abs(*y_target.xpos,*y_target.ypos));
        patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
        patch_diam += patch_diam % 2;
        y_atmos.dim_screens = &patch_diam;
      }
      pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').insert_text(%d,'%s')",
                 nlayers-1,swrite(format="Layer # %d",y_atmos.nscreens));
      pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').set_active(%d)",y_atmos.nscreens-1);
    
    } else if (nlayers <  y_atmos.nscreens) {
      y_atmos.nscreens     -= 1;
      y_atmos.frac         = &((*y_atmos.frac)(1:-1));
      *y_atmos.frac       /= sum(*y_atmos.frac);
      y_atmos.L0           = &((*y_atmos.L0)(1:-1));
      y_atmos.alt          = &((*y_atmos.alt)(1:-1));
      y_atmos.windspeed    = &((*y_atmos.windspeed)(1:-1));
      y_atmos.winddir      = &((*y_atmos.winddir)(1:-1));
      if (y_target == []) y_atmos.dim_screens =&((*y_atmos.dim_screens)(1:-1));
      else {
        max_size = max(abs(*y_target.xpos,*y_target.ypos));
        patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
        patch_diam += patch_diam % 2;
        y_atmos.dim_screens = &patch_diam;
      }
      pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').remove_text(%d)",y_atmos.nscreens);
      pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').set_active(%d)",y_atmos.nscreens-1);
    }
    pyk,swrite(format=atmos_disp._cmd+"y_layers_plot_update(%d)",1);
    if (nlayers != y_atmos.nscreens) update_nlayers,pupdiam,nlayers;
  }
}

func init_layer_prop(nscreen,alt,frac,l0,windspeed,winddir,pupdiam)
{
  extern y_atmos;

  if (y_atmos == []) return;
  else {
    if (y_geom != []) pupdiam = y_geom.pupdiam;
    
    (*y_atmos.alt)(nscreen) = alt;
    (*y_atmos.frac)(nscreen) = frac;
    *y_atmos.frac       /= sum(*y_atmos.frac);
    (*y_atmos.L0)(nscreen) = l0;
    (*y_atmos.windspeed)(nscreen) = windspeed;
    (*y_atmos.winddir)(nscreen) = winddir;
    if (y_target == []) (*y_atmos.dim_screens)(nscreen) = pupdiam+4;
    else {
      max_size = max(abs(*y_target.xpos,*y_target.ypos));
      patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)(nscreen)/y_atmos.pupixsize);
      patch_diam += patch_diam % 2;
      (*y_atmos.dim_screens)(nscreen) = patch_diam;
    }
    pyk,swrite(format=atmos_disp._cmd+"y_layers_plot_update(%d)",1);
  }
}

func remove_layer(nscreen)
{
  extern y_atmos;

  if (y_atmos == []) return;
  else {
    if (y_atmos.nscreens > 1) {
      //      y_atmos.nscreens -=1;      
      if (nscreen != y_atmos.nscreens) {
        if (nscreen == 1) {
          *y_atmos.alt = roll(*y_atmos.alt,-1);
          *y_atmos.frac = roll(*y_atmos.frac,-1);
          *y_atmos.L0 = roll(*y_atmos.L0,-1);
          *y_atmos.windspeed = roll(*y_atmos.windspeed,-1);
          *y_atmos.winddir = roll(*y_atmos.winddir,-1);
        } else {
          alt = *y_atmos.alt;
          y_atmos.alt = &(_(alt(1:nscreen-1),alt(nscreen+1:),alt(nscreen)));
          frac = *y_atmos.frac;
          y_atmos.frac = &(_(frac(1:nscreen-1),frac(nscreen+1:),frac(nscreen)));
          l0 = *y_atmos.L0;
          y_atmos.L0 = &(_(l0(1:nscreen-1),l0(nscreen+1:),l0(nscreen)));
          windspeed = *y_atmos.windspeed;
          y_atmos.windspeed = &(_(windspeed(1:nscreen-1),windspeed(nscreen+1:),windspeed(nscreen)));
          winddir = *y_atmos.winddir;
          y_atmos.winddir = &(_(winddir(1:nscreen-1),winddir(nscreen+1:),winddir(nscreen)));
          dim_screens = *y_atmos.dim_screens;
          y_atmos.dim_screens = &(_(dim_screens(1:nscreen-1),dim_screens(nscreen+1:),dim_screens(nscreen)));
        }
      }
      //update_nlayers(pupdiam,y_atmos.nscreens-1); 
      pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('nlayers').set_value(%d)",y_atmos.nscreens-1);
      pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').set_active(%d)",nscreen-1);
      //pyk,swrite(format=atmos_disp._cmd+"y_layers_plot_update(%d)",1);
    }
  }
}

func load_default_atmos(tconf,pupdiam,teldiam)
{
  extern y_atmos;
  
  if (y_geom != []) pupdiam = y_geom.pupdiam;

  if ((pupdiam == 0) || (pupdiam == [])) pupdiam =132;

  y_atmos  = atmos_struct(); // clean start

  if (tconf == 1) { // one layer ...  
    y_atmos.nscreens     = 1;
    y_atmos.frac         = &([1.0]);
    y_atmos.alt          = &([0.0]);
    y_atmos.L0           = &([100.0]);
    y_atmos.windspeed    = &([10.]);
    y_atmos.winddir      = &([0.]);
    y_atmos.dim_screens  = &([pupdiam+4]);
    y_atmos.pupixsize    = teldiam / pupdiam;
  }
  if (tconf == 2) { // two layers ...  
    y_atmos.nscreens     = 2;
    y_atmos.frac         = &([.8,0.2]);
    y_atmos.alt          = &([0.,10000.]);
    y_atmos.L0           = &([100.,100.]);
    y_atmos.windspeed    = &([10.,35.]);
    y_atmos.winddir      = &([0.,15.]);
    y_atmos.pupixsize    = teldiam / pupdiam;
    if (y_target == []) y_atmos.dim_screens =&(array(pupdiam+4,2));
    else {
      max_size = max(abs(*y_target.xpos,*y_target.ypos));
      patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
      patch_diam += patch_diam % 2;
      y_atmos.dim_screens = &patch_diam;
    }
  }
  if (tconf == 3) { // three layers ...  
    y_atmos.nscreens     = 3;
    y_atmos.frac         = &([.6,0.25,0.15]);
    y_atmos.alt          = &([0.,5000.,10000.]);
    y_atmos.L0           = &([100.,100.,100.]);
    y_atmos.windspeed    = &([10.,20.,35.]);
    y_atmos.winddir      = &([0.,25.,15.]);
    y_atmos.pupixsize    = teldiam / pupdiam;
    if (y_target == []) y_atmos.dim_screens =&(array(pupdiam+4,3));
    else {
      max_size = max(abs(*y_target.xpos,*y_target.ypos));
      patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
      patch_diam += patch_diam % 2;
      y_atmos.dim_screens = &patch_diam;
    }
  }
  if (tconf == 4) { // four layers ...  
    y_atmos.nscreens     = 4;
    y_atmos.frac         = &([.4,0.2,0.25,0.15]);
    y_atmos.alt          = &([0.,1000.,5000.,10000.]);
    y_atmos.L0           = &([100.,100.,100.,100.]);
    y_atmos.windspeed    = &([10.,15.,20.,35.]);
    y_atmos.winddir      = &([0.,10.,25.,15.]);
    y_atmos.pupixsize    = teldiam / pupdiam;
    if (y_target == []) y_atmos.dim_screens =&(array(pupdiam+4,4));
    else {
      max_size = max(abs(*y_target.xpos,*y_target.ypos));
      patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
      patch_diam += patch_diam % 2;
      y_atmos.dim_screens = &patch_diam;
    }
  }
  if (tconf == 5) { // five layers ...  
    y_atmos.nscreens     = 5;
    y_atmos.frac         = &([.4,0.2,0.10,0.15,0.15]);
    y_atmos.alt          = &([0.,1000.,5000.,10000.,15000.]);
    y_atmos.L0           = &([100.,100.,100.,100.,100.]);
    y_atmos.windspeed    = &([10.,15.,20.,35.,50.]);
    y_atmos.winddir      = &([0.,10.,25.,15.,35.]);
    y_atmos.pupixsize    = teldiam / pupdiam;
    if (y_target == []) y_atmos.dim_screens =&(array(pupdiam+4,5));
    else {
      if (y_target != []) 
        max_size = max(abs(*y_target.xpos,*y_target.ypos));
      else max_size = 0.;
      patch_diam = long(pupdiam+4+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
      patch_diam += patch_diam % 2;
      y_atmos.dim_screens = &patch_diam;
    }
  }
  if (tconf < 6) {
    pyk,swrite(format=atmos_disp._cmd+"y_layer_clear(%d)",1);
    //pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').clear()");
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').insert_text(%d,'%s')",0,swrite(format="Layer # %d",y_atmos.nscreens-cc+1));
    pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layer_select').set_active(%d)",0);
    pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('nlayers').set_value(%d)",y_atmos.nscreens);
    
    pyk,swrite(format=atmos_disp._cmd+"y_layers_plot_update(%d)",1);
  }
}

func select_plot_layers(type_plot,nlayer)
{
  extern atmos_disp;
  if (type_plot == -1) pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('layers_plot').set_active(%d)",0);
  if (type_plot == 0) plot_cn2,(*atmos_disp._wins)(2),nlayer=nlayer;
  if (type_plot == 1) plot_speed,(*atmos_disp._wins)(2),nlayer=nlayer;
  if (type_plot == 2) plot_dir,(*atmos_disp._wins)(2),nlayer=nlayer;
}

func plot_cn2(win,nlayer=)
{
  extern y_atmos;

  window,win;
  fma;limits;
  plg,*y_atmos.alt/1000.,*y_atmos.frac,marks=0,width=4;
  plmk,*y_atmos.alt/1000.,*y_atmos.frac,marker=4;
  range,-1,20;limits,0.,0.8;
  xytitles,"R_0_ fraction","Altitude (km)",[0.02,0.02];
  yoga_pltitle,"C_n_^2^",[0.,-0.015];

  if (nlayer != []) {
    plmk,(*y_atmos.alt)(nlayer)/1000.,(*y_atmos.frac)(nlayer),marker=4,color="red";
  }
}


func plot_speed(win,nlayer=)
{
  extern y_atmos;

  window,win;
  fma;limits;
  plvf,(*y_atmos.alt/1000.)*0.0f,*y_atmos.windspeed,*y_atmos.alt/1000.,(*y_atmos.windspeed)*0.0f,hsize=0.1,prop=1,width=4;
  range,-1,20;limits,0.,max(*y_atmos.windspeed)*1.3;
  xytitles,"Wind speed (m/s)","Altitude (km)",[0.015,0.01];
  yoga_pltitle,"Wind Speed",[0.,-0.01];
  
  if (nlayer != [])
    plvf,((*y_atmos.alt)(nlayer)/1000.)*0.0f,(*y_atmos.windspeed)(nlayer),(*y_atmos.alt)(nlayer)/1000.,(*y_atmos.windspeed)(nlayer)*0.0f,hsize=0.1,prop=1,width=4,color="red";
}

func plot_dir(win,nlayer=)
{
  extern y_atmos;

  window,win;
  fma;limits;
  plvf,(*y_atmos.windspeed)*sin(*y_atmos.winddir*dtor),(*y_atmos.windspeed)*cos(*y_atmos.winddir*dtor),(*y_atmos.windspeed)*0.0f,(*y_atmos.windspeed)*0.0f,hsize=0.1,prop=1,width=4;
  range,-1.,max(*y_atmos.windspeed*sin(*y_atmos.winddir*dtor))*1.3;limits,0.,max(*y_atmos.windspeed*cos(*y_atmos.winddir*dtor))*1.3;
  yoga_pltitle,"Wind Directions",[0.,-0.02];

  if (nlayer != [])
    plvf,(*y_atmos.windspeed)(nlayer)*sin((*y_atmos.winddir)(nlayer)*dtor),(*y_atmos.windspeed)(nlayer)*cos((*y_atmos.winddir)(nlayer)*dtor),(*y_atmos.windspeed)(nlayer)*0.0f,(*y_atmos.windspeed)(nlayer)*0.0f,hsize=0.1,prop=1,width=4,color="red";
    
}

func init_gatmos(pupdiam,zenith,teldiam,cobs,r0,freq,only)
{
  extern y_geom,y_tel,y_loop,y_atmos;
  extern g_atmos,g_target;
  extern ipupil;

  if (y_atmos == []) return;
  if (only == []) only = 0;
  
  //pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('init_atmos').set_sensitive(%d)",0);
  
  write,"Doing atmos inits on the GPU";
  
  if (y_geom == []) y_geom   = geom_struct();
  if (y_tel == []) y_tel    = tel_struct();
  if (y_loop == []) y_loop   = loop_struct();

  y_geom.zenithangle = zenith;

  y_tel.diam         = teldiam;
  y_tel.cobs         = cobs;

  y_atmos.r0         = r0;
 
  y_loop.ittime      = 1.0/freq;

  geom_init,pupdiam;

  atmos_init;
  
  write,"... Done with atmos inits on the GPU";

  //pyk,swrite(format=atmos_disp._cmd+"glade.get_widget('init_atmos').set_sensitive(%d)",1);

  update_main_display1,"atmos";

  if (only) g_target = [];
}

atmos_disp = display_struct();
