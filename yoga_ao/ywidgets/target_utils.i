/*
 _____                    _   
|_   _|_ _ _ __ __ _  ___| |_ 
  | |/ _` | '__/ _` |/ _ \ __|
  | | (_| | | | (_| |  __/ |_ 
  |_|\__,_|_|  \__, |\___|\__|
               |___/          
*/
func target_plot_update(ntarget,winnum)
{
  if (y_target == []) return;
  if (*sky_disp._winits == []) return;
  
  if (winnum == []) winnum = (*sky_disp._wins)(6);
  if (!(*sky_disp._winits)(6)) return;
  
  window,winnum;fma;
  plmk,*y_target.ypos,*y_target.xpos,marker=5;
  plmk,(*y_target.ypos)(ntarget),(*y_target.xpos)(ntarget),marker=5,color="red";
  limits,-40,40;
  range,-40,40;
  xytitles,"xpos (arcsec)","ypos (arcsec)",[0.02,0.02];
  yoga_pltitle,"Asterism",[0.,-0.015];
  
  
}
func create_target(ntargets)
{
  extern y_target;
  
  y_target  = target_struct();
  
  y_target.ntargets = ntargets;
  if (ntargets > 1) {
    y_target.xpos     = &(random_n(ntargets));
    y_target.ypos     = &(random_n(ntargets));
  } else {
    y_target.xpos     = &([0.]);
    y_target.ypos     = &([0.]);
  }
  y_target.mag      = &(array(5.,ntargets));
  y_target.lambda   = &(array(1.6,ntargets));

  pyk,swrite(format=sky_disp._cmd+"y_target_clear(%d)",1);    
    //pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').clear()");
  for (cc=1;cc<=ntargets;cc++)
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",0,
               swrite(format="Source # %d",ntargets-cc+1));
  
  update_target_prop,0;
  pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
}

func update_target_prop(ntarget)
{
  extern y_target;
  
  if (y_target == []) {
    pyk,swrite(format=sky_disp._cmd+"y_init_target(%d)",1);
  } else {
    if (ntarget > y_target.ntargets) error,"not a valid target";
    else {
      pyk,swrite(format=sky_disp._cmd+"y_update_target_gui(%f, %f, %f, %f)",
                 (*y_target.xpos)(ntarget),(*y_target.ypos)(ntarget),
                 (*y_target.lambda)(ntarget),(*y_target.mag)(ntarget));
      pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
    }
  }
}

func update_ntargets(ntargets)
{
  extern y_target;
  
  if (y_target != []) {
    if (ntargets > y_target.ntargets) {
      y_target.ntargets += 1;
      y_target.xpos      = &(_((*y_target.xpos),0.));
      y_target.ypos      = &(_((*y_target.ypos),0.));
      y_target.lambda    = &(_((*y_target.lambda),1.6));
      y_target.mag       = &(_((*y_target.mag),5));
      
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",
                 ntargets-1,swrite(format="Source # %d",y_target.ntargets));
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').set_active(%d)",y_target.ntargets-1);
    
    } else if (ntargets <  y_target.ntargets) {
      y_target.ntargets -= 1;
      y_target.xpos      = &((*y_target.xpos)(1:-1));
      y_target.ypos      = &((*y_target.ypos)(1:-1));
      y_target.lambda    = &((*y_target.lambda)(1:-1));
      y_target.mag       = &((*y_target.mag)(1:-1));

      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').remove_text(%d)",y_target.ntargets);
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').set_active(%d)",y_target.ntargets-1);
    }
    pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
    if (ntargets != y_target.ntargets) update_ntargets,ntargets;
  }
}

func init_target_prop(ntarget,mag,xpos,ypos,lambda)
{
  extern y_target;

  if (y_target == []) return;
  else {
    (*y_target.mag)(ntarget) = mag;
    (*y_target.xpos)(ntarget) = xpos;
    (*y_target.ypos)(ntarget) = ypos;
    (*y_target.lambda)(ntarget) = lambda;
    pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
  }
}

func remove_target(ntarget)
{
  extern y_target;

  if (y_target == []) return;
  else {
    if (y_target.ntargets > 1) {
      //      y_atmos.nscreens -=1;      
      if (ntarget != y_target.ntargets) {
        if (ntarget == 1) {
          *y_target.xpos = roll(*y_target.xpos,-1);
          *y_target.ypos = roll(*y_target.ypos,-1);
          *y_target.lambda = roll(*y_target.lambda,-1);
          *y_target.mag = roll(*y_target.mag,-1);
        } else {
          xpos = *y_target.xpos;
          y_target.xpos = &(_(xpos(1:ntarget-1),xpos(ntarget+1:),xpos(ntarget)));
          ypos = *y_target.ypos;
          y_target.ypos = &(_(ypos(1:ntarget-1),ypos(ntarget+1:),ypos(ntarget)));
          lambda = *y_target.lambda;
          y_target.lambda = &(_(lambda(1:ntarget-1),lambda(ntarget+1:),lambda(ntarget)));
          mag = *y_target.mag;
          y_target.mag = &(_(mag(1:ntarget-1),mag(ntarget+1:),mag(ntarget)));
        }
      }
      //update_nlayers(pupdiam,y_atmos.nscreens-1); 
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('ntargets').set_value(%d)",y_target.ntargets-1);
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').set_active(%d)",ntarget-1);
      //pyk,swrite(format=sky_disp._cmd+"y_layers_plot_update(%d)",1);
    }
  }
}

func load_default_target(tconf)
{
  extern y_target;
  
  y_target  = target_struct(); // clean start

  if (tconf == 1) { // one source ...  
    y_target.ntargets = 1;
    y_target.xpos     = &([0.0]);
    y_target.ypos     = &([0.0]);
    y_target.lambda   = &([1.6]);
    y_target.mag      = &([5.]);
  }
  if (tconf == 2) { // two sources ...  
    y_target.ntargets = 2;
    y_target.xpos     = &([20.0,-20.0]);
    y_target.ypos     = &([0.0,0.0]);
    y_target.lambda   = &([1.6,1.6]);
    y_target.mag      = &([5.,5.]);
  }
  if (tconf == 3) { // three sources ...  
    y_target.ntargets = 3;
    y_target.xpos     = &([20.0,-10.0,-10.0]);
    y_target.ypos     = &([0.0,17.3,-17.3]);
    y_target.lambda   = &([1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.]);
  }
  if (tconf == 4) { // four sources ...  
    y_target.ntargets = 4;
    y_target.xpos     = &([20.0,0.0,-20.0,0.0]);
    y_target.ypos     = &([0.0,20.0,0.0,-20.0]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.]);
  }
  if (tconf == 5) { // five sources ...  
    y_target.ntargets = 5;
    y_target.xpos     = &([20.0,6.2,-16.2,-16.2,6.2]);
    y_target.ypos     = &([0.0,19.0,11.7,-11.7,-19.0]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.,5.]);
  }
  if (tconf == 6) { // six sources ...  
    y_target.ntargets = 6;
    y_target.xpos     = &([20.0,10.0,-10.0,-20.0,-10.,10.]);
    y_target.ypos     = &([0.0,17.3,17.3,0.0,-17.3,-17.3]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.,5.,5.]);
  }
  if (tconf == 7) { // seven sources ...  
    y_target.ntargets = 7;
    y_target.xpos     = &([20.0,10.,-10.,-20.0,-10.,10.,0.0]);
    y_target.ypos     = &([0.0,17.3,17.3,0.0,-17.3,-17.3,0.0]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.,5.,5.,5.]);
  }
  if (tconf < 8) {
    pyk,swrite(format=sky_disp._cmd+"y_target_clear(%d)",1);    
    //pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').clear()");
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",0,swrite(format="Source # %d",y_target.ntargets-cc+1));
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').set_active(%d)",0);
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('ntargets').set_value(%d)",y_target.ntargets);
    
    pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
  }
}
