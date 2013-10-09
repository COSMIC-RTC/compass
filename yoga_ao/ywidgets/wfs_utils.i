func lgs_plot_update(ntarget,plot_type,winnum)
{
  if (y_wfs == []) return;
  if (!(y_wfs(ntarget).gsalt > 0)) return;
  
  if (!(*sky_disp._winits)(3)) return;
  if (winnum == []) winnum = (*sky_disp._wins)(3);
  
  window,winnum;fma;limits;
  if (plot_type == 1) {
    plg,*y_wfs(ntarget)._profna,*y_wfs(ntarget)._altna,marks=0;
    limits,min(*y_wfs(ntarget)._altna),max(*y_wfs(ntarget)._altna);
    range,min(*y_wfs(ntarget)._profna),max(*y_wfs(ntarget)._profna);
    xytitles,"Atitude (m)","Na density",[0.02,0.02];
    yoga_pltitle,"Asterism",[0.,-0.015];
  }
  if (plot_type == 2)
    pli,(*y_wfs(ntarget)._beam)(,-) * (*y_wfs(ntarget)._beam)(-,);
  if (plot_type == 3) {
    mimg = sensors_getdata(g_wfs,ntarget-1,"bincube");
    pli,mimg(,,1);
  }
}

func update_wfs_prop(numwfs)
{
  extern wfs_disp,y_wfs;
  
  if (y_wfs == []) {
    pyk,swrite(format=wfs_disp._cmd+"y_init_wfs(%d)",1);
  } else {
    if (numwfs > numberof(y_wfs)) error,"not a valid wfs";
    else {
      if (y_wfs(numwfs).gsalt > 0.)  pyk,swrite(format=wfs_disp._cmd+"y_update_lgs(%d)",1);
      else pyk,swrite(format=wfs_disp._cmd+"y_update_lgs(%d)",0);
      if (y_wfs(numwfs).type == "sh") typewfs = 0;
      if (y_wfs(numwfs).type == "pyr") typewfs = 1;
      
      pyk,swrite(format=wfs_disp._cmd+"y_update_wfs_gui(%d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f)",
                 typewfs,(y_wfs.nxsub)(numwfs),(y_wfs.npix)(numwfs),
                 (y_wfs.pixsize)(numwfs),(y_wfs.fracsub)(numwfs),
                 (y_wfs.xpos)(numwfs),(y_wfs.ypos)(numwfs),(y_wfs.lambda)(numwfs),
                 (y_wfs.gsmag)(numwfs),log10((y_wfs.zerop)(numwfs)),
                 (y_wfs.optthroughput)(numwfs),(y_wfs.noise)(numwfs));
      
    }
  }
}

func init_wfs_prop(typewfs,numwfs,nsub,npix,pixsize,mag,xpos,ypos,lambda,frac,zp,throughput,noise)
{
  extern y_wfs;

  if (y_wfs == []) return;
  else {
    if (typewfs == 0) y_wfs(numwfs).type = "sh";
    else y_wfs(numwfs).type = "pyr";
    
    y_wfs(numwfs).nxsub         = nsub;
    y_wfs(numwfs).npix          = npix;
    y_wfs(numwfs).pixsize       = pixsize;
    y_wfs(numwfs).fracsub       = frac;
    y_wfs(numwfs).xpos          = xpos;
    y_wfs(numwfs).ypos          = ypos;
    y_wfs(numwfs).lambda        = lambda;
    y_wfs(numwfs).gsmag         = mag;
    y_wfs(numwfs).optthroughput = throughput;
    y_wfs(numwfs).zerop       = 10^zp;
    y_wfs(numwfs).noise       = noise;
  }
}

func init_wfs_prop_lgs(gsalt,lltx,llty,power,wreturn,proftype,beam)
{
  extern y_wfs;

  if (y_wfs == []) return;
  else {
    y_wfs(numwfs).gsalt            = gsalt*1.e3;
    y_wfs(numwfs).lltx             = lltx;
    y_wfs(numwfs).llty             = llty;
    y_wfs(numwfs).laserpower       = power;
    y_wfs(numwfs).lgsreturnperwatt = wreturn;
    y_wfs(numwfs).proftype         = proftype;
    y_wfs(numwfs).beamsize         = beam;
  }
}

func create_wfs(numwfs,typewfs,teldiam)
{
  extern y_wfs;
  
  y_wfs  = array(wfs_struct(),numwfs);
  
  if (numwfs > 1) {
    xpos     = random_n(numwfs);
    ypos     = random_n(numwfs);
  } else {
    xpos     = 0.;
    ypos     = 0.;
  }
  for (i=1;i<=numwfs;i++) {
    if (typewfs == 0) y_wfs(i).type = "sh";
    else y_wfs(i).type = "pyr";

    y_wfs(i).nxsub         = long(teldiam*1./0.5);
    y_wfs(i).npix          = 6;
    y_wfs(i).pixsize       = 0.3;
    y_wfs(i).fracsub       = 0.8;
    y_wfs(i).xpos          = xpos(i);
    y_wfs(i).ypos          = ypos(i);
    y_wfs(i).lambda        = 0.5;
    y_wfs(i).gsmag         = 5.;
    y_wfs(i).optthroughput = 0.5;
    y_wfs(i).zerop         = 1.e11;
    y_wfs(i).noise         = -1;
  }
  
  pyk,swrite(format=wfs_disp._cmd+"y_wfs_clear(%d)",1);    
  //pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').clear()");
  for (cc=1;cc<=numwfs;cc++)
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",
               0,swrite(format="WFS # %d",numwfs-cc+1));
  
  update_wfs_prop,0;
  if (y_wfs(1).type == "pyr") load_default_centro,"pyr";
  else load_default_centro,"cog";
}

func update_nwfs(numwfs,teldiam)
{
  extern y_wfs;
  
  if (y_wfs != []) {
    if (numwfs > numberof(y_wfs)) {
      tmp = y_wfs;
      tmp2 = array(wfs_struct,numberof(y_wfs)+1);
      tmp2(1:numberof(y_wfs)) = y_wfs;
      y_wfs = tmp2;
      
      if (typewfs == 0) y_wfs(0).type = "sh";
      else y_wfs(0).type = "pyr";

      y_wfs(0).nxsub         = long(teldiam * 1/0.5);
      y_wfs(0).npix          = 6;
      y_wfs(0).pixsize       = 0.3;
      y_wfs(0).fracsub       = 0.8;
      y_wfs(0).xpos          = 0.0;
      y_wfs(0).ypos          = 0.0;
      y_wfs(0).lambda        = 0.5;
      y_wfs(0).gsmag           = 5.;
      y_wfs(0).optthroughput = 0.5;
      y_wfs(0).zerop         = 1.e11;
      y_wfs(0).noise         = -1;
      
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",
                 numberof(y_wfs)-1,swrite(format="WFS # %d",numberof(y_wfs)));
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",numberof(y_wfs)-1);
    
      if (y_rtc != []) {
        tmp = (*y_rtc.centroiders);
        tmp2 = array(centroider_struct,numberof(tmp)+1);
        tmp2(1:-1) = (*y_rtc.centroiders);
        if (typewfs == 0) tmp2(0).type   = "cog";
        else tmp2(0).type   = "pyr";
        
        tmp2(0).nwfs   = numberof(y_wfs);
        y_rtc.centroiders = &tmp2;
        y_rtc.nwfs       += 1;
        pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').insert_text(%d,'%s')",
                   numberof(y_wfs)-1,swrite(format="centro # %d",numberof(y_wfs)));
        pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",numberof(y_wfs)-1);
      }
    } else if (numwfs <  numberof(y_wfs)) {
      y_wfs = y_wfs(1:-1);

      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').remove_text(%d)",numberof(y_wfs));
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",numberof(y_wfs)-1);
      if (y_rtc != []) {
        tmp = (*y_rtc.centroiders)(1:-1);
        y_rtc.centroiders = &tmp;
        y_rtc.nwfs       -= 1;
        pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",numberof(y_wfs)-1);
      }
    } else {
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",0);
      if (y_rtc != []) {
        pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",0);
      }
    }
    if (numwfs != numberof(y_wfs)) update_nwfs,numwfs,teldiam;
  }
}

func remove_wfs(numwfs)
{
  extern y_wfs,y_rtc;

  if (y_wfs == []) return;
  else {
    if (numberof(y_wfs)>1) {
      if (numwfs != numberof(y_wfs)) {
        if (numwfs == 1) {
          tmp = y_wfs(2:);
          y_wfs = tmp;
          if (y_rtc != []) {
            tmp2 = (*y_rtc.centroiders)(2:);
            y_rtc.centroiders = &tmp2;
          }
        } else {
          tmp = y_wfs(1:-1);
          tmp(1:numwfs-1) = y_wfs(1:numwfs-1);
          tmp(numwfs:-1) = y_wfs(numwfs+1:);
          y_wfs = tmp;
          if (y_rtc != []) {
            tmp2 = (*y_rtc.centroiders)(1:-1);
            tmp2(1:numwfs-1)  = (*y_rtc.centroiders)(1:numwfs-1);
            tmp2(numwfs:-1)   = (*y_rtc.centroiders)(numwfs+1:);
            y_rtc.centroiders = &tmp2;
            y_rtc.nwfs       -= 1;
          }
        }
      }
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('nwfs').set_value(%d)",numberof(y_wfs)-1);
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",numberof(y_wfs)-1);
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('centro_select').set_active(%d)",numberof(y_wfs)-1);
    }
  }
}


func load_default_wfs(tconf,teldiam)
{
  extern y_wfs;
  
  if (tconf == 1) { // one wfs 50cm / subap  
    y_wfs                  = array(wfs_struct(),1); // clean start
    y_wfs(1).type          = "sh";
    y_wfs(1).nxsub         = long(teldiam * 1/0.5);
    y_wfs(1).npix          = 6;
    y_wfs(1).pixsize       = 0.3;
    y_wfs(1).fracsub       = 0.8;
    y_wfs(1).xpos          = 0.0;
    y_wfs(1).ypos          = 0.0;
    y_wfs(1).lambda        = 0.5;
    y_wfs(1).gsmag         = 5.;
    y_wfs(1).optthroughput = 0.5;
    y_wfs(1).zerop         = 1.e11;
    y_wfs(1).noise         = -1;
    y_wfs(1).gsalt         = 0.;
  }
  if (tconf == 2) { // one wfs 60cm / subap  
    y_wfs                  = array(wfs_struct(),1); // clean start
    y_wfs(1).type          = "sh";
    y_wfs(1).nxsub         = long(teldiam * 1/0.6);
    y_wfs(1).npix          = 6;
    y_wfs(1).pixsize       = 0.3;
    y_wfs(1).fracsub       = 0.8;
    y_wfs(1).xpos          = 0.0;
    y_wfs(1).ypos          = 0.0;
    y_wfs(1).lambda        = 0.5;
    y_wfs(1).gsmag         = 5.;
    y_wfs(1).optthroughput = 0.5;
    y_wfs(1).zerop         = 1.e11;
    y_wfs(1).noise         = -1;
    y_wfs(1).gsalt         = 0.;
  }
  if (tconf == 3) { // one pyr 
    y_wfs                  = array(wfs_struct(),1); // clean start
    y_wfs(1).type          = "pyr";
    y_wfs(1).nxsub         = long(teldiam * 1/0.5);
    y_wfs(1).npix          = 4;
    y_wfs(1).fssize        = 1.6;
    y_wfs(1).fracsub       = 1.0;
    y_wfs(1).xpos          = 0.0;
    y_wfs(1).ypos          = 0.0;
    y_wfs(1).lambda        = 0.5;
    y_wfs(1).gsmag         = 5.;
    y_wfs(1).optthroughput = 0.5;
    y_wfs(1).zerop         = 1.e11;
    y_wfs(1).noise         = -1;
    y_wfs(1).gsalt         = 0.;
    
    y_wfs(1).fstop         = "round";
    y_wfs(1).pyr_npts      = 16;
    y_wfs(1).pyr_ampl      = 0.45;
  }
  if (tconf == 4) { // one pyr 
    y_wfs                  = array(wfs_struct(),1); // clean start
    y_wfs(1).type          = "pyr";
    y_wfs(1).nxsub         = long(teldiam * 1/0.2);
    y_wfs(1).npix          = 4;
    y_wfs(1).fssize        = 1.6;
    y_wfs(1).fracsub       = 1.0;
    y_wfs(1).xpos          = 0.0;
    y_wfs(1).ypos          = 0.0;
    y_wfs(1).lambda        = 0.5;
    y_wfs(1).gsmag         = 5.;
    y_wfs(1).optthroughput = 0.5;
    y_wfs(1).zerop         = 1.e11;
    y_wfs(1).noise         = -1;
    y_wfs(1).gsalt         = 0.;
    
    y_wfs(1).fstop         = "round";
    y_wfs(1).pyr_npts      = 16;
    y_wfs(1).pyr_ampl      = 0.45;
  }
  if (tconf == 5) { // one wfs 20cm / subap  
    y_wfs                  = array(wfs_struct(),1); // clean start
    y_wfs(1).type          = "sh";
    y_wfs(1).nxsub         = long(teldiam * 1/0.2);
    y_wfs(1).npix          = 6;
    y_wfs(1).pixsize       = 0.3;
    y_wfs(1).fracsub       = 0.8;
    y_wfs(1).xpos          = 0.0;
    y_wfs(1).ypos          = 0.0;
    y_wfs(1).lambda        = 0.5;
    y_wfs(1).gsmag         = 5.;
    y_wfs(1).optthroughput = 0.5;
    y_wfs(1).zerop         = 1.e11;
    y_wfs(1).noise         = -1;
    y_wfs(1).gsalt         = 0.;
  }
  if (tconf == 6) { // 1 wfs lgs
    y_wfs                  = array(wfs_struct(),1); // clean start
    y_wfs.type          = "sh";
    y_wfs.nxsub         = long(teldiam * 1/0.5);
    y_wfs.npix          = 6;
    y_wfs.pixsize       = 0.3;
    y_wfs.fracsub       = 0.8;
    y_wfs.lambda        = 0.5;
    y_wfs.gsmag         = 5.;
    y_wfs.optthroughput = 0.5;
    y_wfs.zerop         = 1.e11;
    y_wfs.noise         = -1;
    
    y_wfs(1).xpos          = 0.;
    y_wfs(1).ypos          = 0.;

    y_wfs(1).gsalt            = 90*1.e3;
    y_wfs(1).lltx             = 0.;
    y_wfs(1).llty             = 0.;
    y_wfs(1).laserpower       = 10;
    y_wfs(1).lgsreturnperwatt = 1.e3;
    y_wfs(1).proftype         = "Exp";
    y_wfs(1).beamsize         = 0.8;
  }
  if (tconf == 7) { // 3 wfs
    y_wfs                  = array(wfs_struct(),3); // clean start
    y_wfs.type          = "sh";
    y_wfs.nxsub         = long(teldiam * 1/0.5);
    y_wfs.npix          = 6;
    y_wfs.pixsize       = 0.3;
    y_wfs.fracsub       = 0.8;
    y_wfs.lambda        = 0.5;
    y_wfs.gsmag         = 5.;
    y_wfs.optthroughput = 0.5;
    y_wfs.zerop         = 1.e11;
    y_wfs.noise         = -1;
    
    y_wfs(1).xpos          = 40.*cos(2*pi/3.);
    y_wfs(1).ypos          = 40.*sin(2*pi/3.);

    y_wfs(2).xpos          = 40.*cos(4*pi/3.);
    y_wfs(2).ypos          = 40.*sin(4*pi/3.);

    y_wfs(3).xpos          = 40.
    y_wfs(3).ypos          = 0;

  }
  if (tconf == 8) { // canary
    y_wfs                  = array(wfs_struct(),4); // clean start
    y_wfs.type          = "sh";
    y_wfs.nxsub         = long(teldiam * 1/0.5);
    y_wfs.npix          = 16;
    y_wfs.pixsize       = 0.3;
    y_wfs.fracsub       = 0.8;
    y_wfs.lambda        = 0.5;
    y_wfs.gsmag         = 5.;
    y_wfs.optthroughput = 0.5;
    y_wfs.zerop         = 1.e11;
    y_wfs.noise         = -1;
    
    y_wfs(1).xpos          = 40.*cos(2*pi/3.);
    y_wfs(1).ypos          = 40.*sin(2*pi/3.);

    y_wfs(2).xpos          = 40.*cos(4*pi/3.);
    y_wfs(2).ypos          = 40.*sin(4*pi/3.);

    y_wfs(3).xpos          = 40.
    y_wfs(3).ypos          = 0;

    y_wfs(4).xpos          = 0.
    y_wfs(4).ypos          = 0;

  }
  if (tconf < 8) {
    
    pyk,swrite(format=wfs_disp._cmd+"y_wfs_clear(%d)",1);    
    //pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').clear()");
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",0);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('nwfs').set_value(%d)",numberof(y_wfs));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('nwfs').set_value(%d)",numberof(y_wfs));
    
  }
  
  if (y_wfs(1).type == "sh")  load_default_centro,"cog";
  else load_default_centro,"pyr";
  
  pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('default_centro').set_active(%d)",0);
}
