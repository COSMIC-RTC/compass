// Environment check
yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data/";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

require,yoga_ao_top+"/yorick/yoga_aolib.i"
require,yoga_ao_top+"/yorick/yoga_ao_ystruct.i"
require,yoga_ao_top+"/yorick/yoga_ao_utils.i"
require,yoga_ao_top+"/yorick/yoga_turbu.i"
require,yoga_ao_top+"/yorick/yoga_wfs.i"
require,yoga_ao_top+"/yorick/yoga_rtc.i"
require,yoga_ao_top+"/yorick/yoga_dm.i"
require,yoga_ao_top+"/yorick/yoga_kl.i"

func read_parfile(filename)
/* DOCUMENT read_parfile
   read_parfile,filename
     
   reads yoga parameter file filename

   SEE ALSO:
 */
{
  extern y_geom,y_atmos,y_tel,y_target,y_loop,y_wfs,y_rtc,y_dm;

  if (!fileExist(filename)) {
    exit,swrite(format="Could not find parameter file %s !\n",filename);}

  y_geom   = geom_struct();
  y_atmos  = atmos_struct();
  y_tel    = tel_struct();
  y_target = target_struct();
  y_loop   = loop_struct();
  y_rtc    = rtc_struct();
  
  require,filename;
}


func geom_init(pupdiam)
/* DOCUMENT geom_init
   geom_init,pupdiam
     
   inits simulation geometry, depending on pupdiam
   the linear number of pixels in the pupil
   
   requires 3 externals : y_atmos, y_geom and g_atmos
   y_atmos : a y_struct for the atmosphere
   y_geom  : a y_struct for the geometry of the simulation
   g_atmos : a yAtmos object on the gpu

  SEE ALSO:
 */
{
  extern y_atmos,y_geom,g_atmos;

  if (y_geom == []) y_geom   = geom_struct();
  y_geom.pupdiam = pupdiam;

  // power of 2 greater than pupdiam
  y_geom.ssize  = long(2^ceil(log(y_geom.pupdiam)/log(2)+1));

  // using images centered on 1/2 pixels
  y_geom.cent   = y_geom.ssize / 2 + 0.5;
  
  // valid pupil geometry
  _p        = y_geom.pupdiam;
  _p1       = long(ceil(y_geom.cent-_p/2.));
  _p2       = long(floor(y_geom.cent+_p/2.));
  _p        = _p2-_p1+1;

  y_geom.pupdiam = _p;
  y_geom._p1     = _p1;
  y_geom._p2     = _p2;

  // pupil with 2 pixels bands on each side
  y_geom._n      = _p+4;
  y_geom._n1     = _p1-2;
  y_geom._n2     = _p2+2;

  // Initialize pupil array:

  // large pupil (used for image formation)
  ipupil = float(make_pupil(y_geom.ssize,y_geom.pupdiam,xc=y_geom.cent,yc=y_geom.cent,\
                            cobs=y_tel.cobs));

  // useful pupil 
  spupil = float(make_pupil(y_geom.pupdiam,y_geom.pupdiam,xc=y_geom.pupdiam/2+0.5,\
                            yc=y_geom.pupdiam/2+0.5,  cobs=y_tel.cobs));

  // useful pupil + 4 pixels
  mpupil = float(make_pupil(y_geom._n,y_geom.pupdiam,xc=y_geom._n/2+0.5,yc=y_geom._n/2+0.5,\
                            cobs=y_tel.cobs));

  y_geom._ipupil = &ipupil;
  y_geom._spupil = &spupil;
  y_geom._mpupil = &mpupil;
}

func atmos_init(void)
/* DOCUMENT atmos_init
   atmos_init
     
   inits a yAtmos object on the gpu
   no input parameters
   
   requires 2 externals + 2 optional : y_atmos and y_geom + y_target and y_wfs
   y_atmos  : a y_struct for the atmosphere
   y_geom   : a y_struct for the geometry of the simulation
   y_target : a y_struct for the targets
   y_wfs    : a y_struct for the sensors
   creates 1 external :
   g_atmos : a yAtmos object on the gpu

  SEE ALSO:
 */
{
  extern g_atmos;
  
  // ajust layers alt using zenith angle
  y_atmos.alt = &(*y_atmos.alt / cos(y_geom.zenithangle*dtor));
  // pixel size in meter
  y_atmos.pupixsize = y_tel.diam / y_geom.pupdiam; 

  // compute total fov using targets and wfs gs
  if ((y_wfs != []) && (y_target != []))
    max_size = max(_(abs(*y_target.xpos,*y_target.ypos),abs(y_wfs.xpos,y_wfs.ypos)));
  else {
    if (y_target != [])
      max_size = max(abs(*y_target.xpos,*y_target.ypos));
    else if (y_wfs != [])
      max_size = max(abs(y_wfs.xpos,y_wfs.ypos));
    else max_size = 0.;
  }
  
  // compute corresponding meta-pupil diameters at each alt
  patch_diam = long(y_geom._n+2*(max_size*4.84814e-6*(*y_atmos.alt))/y_atmos.pupixsize+4);
  patch_diam += patch_diam % 2;

  y_atmos.dim_screens = &patch_diam;

  // compute phase screens speed in pixels / iteration
  deltax  = y_geom.pupdiam/y_tel.diam*(*y_atmos.windspeed)*cos(dtor*y_geom.zenithangle)*y_loop.ittime;
  deltay  = deltax*sin(dtor*(*y_atmos.winddir));
  deltax  = deltax*cos(dtor*(*y_atmos.winddir));

  y_atmos.deltax = &float(deltax);
  y_atmos.deltay = &float(deltay);

  if (*y_atmos.L0 == []) y_atmos.L0 = &((1.e5)(-:1:y_atmos.nscreens)); // infinite L0
  else {
    if (numberof(*y_atmos.L0) == 1)
      y_atmos.L0 = &((*y_atmos.L0)(-:1:y_atmos.nscreens));
    for (cc=1;cc <= y_atmos.nscreens;cc++) {
      // L0 should be gien in meters by the user
      // computing the fraction of the telescope diameter
      // translate that into pixels in the pupil for extrude
      frac_l0 = y_tel.diam / ((*y_atmos.L0)(cc));
      (*y_atmos.L0)(cc) = y_geom.pupdiam / frac_l0;
    }
  }
  
  // create atmos object on the gpu
  g_atmos = yoga_atmos_create(y_atmos.nscreens,y_atmos.r0,*y_atmos.L0,y_atmos.pupixsize,*y_atmos.dim_screens,
                              *y_atmos.frac,*y_atmos.alt,*y_atmos.windspeed,*y_atmos.winddir,
                              *y_atmos.deltax,*y_atmos.deltay,*y_geom._spupil);

}

func wfs_init(void)
/* DOCUMENT wfs_init
   wfs_init
     
   inits a ySeznsors object on the gpu
   no input parameters
   
   requires 2 externals : y_geom +and y_wfs
   y_geom   : a y_struct for the geometry of the simulation
   y_wfs    : a y_struct for the sensors
   creates 1 external :
   g_wfs    : a ySensors object on the gpu

  SEE ALSO:
 */
{
  extern y_geom;
  extern g_wfs;
  // first get the wfs with max # of subaps
  // we'll derive the geometry from the requirements in terms of sampling
  if (noneof(y_wfs.type == "pyr")) indmax = wheremax(y_wfs.nxsub)(1);
  else {
    if (anyof(y_wfs.type == "sh"))
      indmax = wheremax(y_wfs(where(y_wfs.type == "sh")).nxsub)(1);
    else
      indmax = wheremax(y_wfs.nxsub)(1);
  }
  
  // init geometry
  if (noneof(y_wfs.type == "geo"))
    init_wfs_geom,indmax,init=1;
  else
    init_wfs_geom,indmax,init=0;
  
  // do the same for other wfs
  for (i=1;i<=numberof(y_wfs);i++) {
    if (i != wheremax(y_wfs.nxsub)(1)) {
      init_wfs_geom,i;
    }
  }
  
  // create sensor object on gpu
   if (y_wfs(1).type == "sh")
     g_wfs = yoga_sensors(numberof(y_wfs),y_wfs(1).type,y_wfs.nxsub,y_wfs._nvalid,y_wfs.npix,y_wfs._pdiam,
                          y_wfs._nrebin,y_wfs._Nfft,y_wfs._Ntot,y_geom._n,y_wfs._subapd,
                          y_wfs._nphotons,y_wfs.gsalt > 0);
   if (y_wfs(1).type == "pyr")
     g_wfs = yoga_sensors(numberof(y_wfs),y_wfs(1).type,y_wfs.nxsub,y_wfs._nvalid,y_wfs.npix,y_wfs._pdiam,
                          y_wfs._nrebin,y_wfs._Nfft,y_wfs._Ntot,y_wfs(1).pyr_npts,y_wfs._subapd,
                          y_wfs._nphotons,y_wfs.gsalt > 0);
   if (y_wfs(1).type == "geo")
     g_wfs = yoga_sensors(numberof(y_wfs),y_wfs(1).type,y_wfs.nxsub,y_wfs._nvalid,y_wfs._pdiam,
                          y_geom._n,y_wfs._subapd);
  
  // init sensor gs object on gpu
   if (y_wfs(1).type == "geo") {
     sensors_initgs,g_wfs,y_wfs.xpos,y_wfs.ypos,y_wfs.lambda,[0](-::numberof(y_wfs)-1),
       (y_geom._n)(-::numberof(y_wfs)-1),[-1](-::numberof(y_wfs)-1);
   } else {
     sensors_initgs,g_wfs,y_wfs.xpos,y_wfs.ypos,y_wfs.lambda,y_wfs.gsmag,
       (y_geom._n)(-::numberof(y_wfs)-1),y_wfs.noise;
   }
  // fill sensor object with data
  for (i=1;i<=numberof(y_wfs);i++) {
    
    if (y_wfs(i).type == "sh")
      sensors_initarr,g_wfs,i-1,int(*y_wfs(i)._phasemap),int(*y_wfs(i)._hrmap),
        int(*y_wfs(i)._binmap),float(*y_wfs(i)._halfxy),float(*y_geom._mpupil),
        (*y_wfs(i)._fluxPerSub)(where(*y_wfs(i)._isvalid)),int(*y_wfs(i)._isvalid),
        int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1),int(*y_wfs(i)._istart+1),
        int(*y_wfs(i)._jstart+1),float(*y_wfs(i)._ftkernel);
    
    if (y_wfs(i).type == "pyr") {
      tmp = array(float,2,y_wfs(i)._Ntot,y_wfs(i)._Ntot);
      tmp(1,,) = (*y_wfs(i)._halfxy).re;
      tmp(2,,) = (*y_wfs(i)._halfxy).im;
      
      tmp2 = array(float,2,y_wfs(i)._Nfft,y_wfs(i)._Nfft);
      tmp2(1,,) = (*y_wfs(i)._pyr_offsets).re;
      tmp2(2,,) = (*y_wfs(i)._pyr_offsets).im;
      
      sensors_initarr,g_wfs,i-1,float(tmp),float(tmp2),float(*y_wfs(i)._submask),float(*y_geom._mpupil),
        int(*y_wfs(i)._isvalid),int(*y_wfs(i)._pyr_cx),int(*y_wfs(i)._pyr_cy),float(*y_wfs(i)._hrmap),
        int(*y_wfs(i)._phasemap),int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1);
    }
    
    if (y_wfs(i).type == "geo")
      sensors_initarr,g_wfs,i-1,int(*y_wfs(i)._phasemap),float(*y_wfs(i)._halfxy),float(*y_geom._mpupil),
        (*y_wfs(i)._fluxPerSub)(where(*y_wfs(i)._isvalid)),int(*y_wfs(i)._isvalid),
        int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1),int(*y_wfs(i)._istart+1),
        int(*y_wfs(i)._jstart+1);
  }

  // lgs case
  for (cc=1;cc<=numberof(y_wfs);cc++) {
    if (y_wfs(cc).gsalt > 0) {
      // lgs mode requested
      if (y_wfs(cc).proftype == "None") y_wfs(cc).proftype = "Gauss1";
      if (y_wfs(cc).proftype == "Gauss1")
        profilename = "allProfileNa_withAltitude_1Gaussian.fits";
      if (y_wfs(cc).proftype == "Gauss2")
        profilename = "allProfileNa_withAltitude_2Gaussians.fits";
      if (y_wfs(cc).proftype == "Gauss3")
        profilename = "allProfileNa_withAltitude_3Gaussians.fits";
      if (y_wfs(cc).proftype == "Exp")
        profilename = "allProfileNa_withAltitude.fits";

      prof=fits_read(YOGA_AO_SAVEPATH+profilename);
      h    = prof(,1);
      prof = prof(,2:)(,avg);
      y_wfs(cc)._altna  = &h;
      y_wfs(cc)._profna = &prof;
      
      // init sensor gs object with necessary data
      prep_lgs_prof,cc,prof,h,y_wfs(cc).beamsize;
      //sensors_loadkernels,g_wfs,cc-1,float(*y_wfs(cc)._lgskern);
    }
  }
}

func target_init(void)
/* DOCUMENT target_init
   target_init
     
   inits a yTarget object on the gpu
   no input parameters
   
   requires 4 externals + 1 optional : y_geom y_atmos and y_target + y_wfs
   y_geom   : a y_struct for the geometry of the simulation
   y_atmos  : a y_struct for the atmosphere
   y_target : a y_struct for the targets
   y_wfs    : a y_struct for the sensors
   creates 1 external :
   g_target    : a yTarget object on the gpu

  SEE ALSO:
 */
{
  extern g_target;
  
  type = "atmos";

  if (y_wfs != []) {
    if ((y_wfs != []) && (g_wfs != [])) {
      for (cc=1;cc<=numberof(y_wfs);cc++) {
        if (y_atmos != []) {
          for (dd=1;dd<=y_atmos.nscreens;dd++) {
            xoff = (y_wfs.xpos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
            yoff = (y_wfs.ypos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
            xoff = float(xoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
            yoff = float(yoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
            sensors_addlayer,g_wfs,cc-1,type,(*y_atmos.alt)(dd),xoff,yoff;
          }
        }
        if (y_dm != []) {
          for (dd=1;dd<=numberof(y_dm);dd++) {
            dims = y_dm(dd)._n2 - y_dm(dd)._n1 + 1;
            dim  = dimsof(*y_geom._mpupil)(2);
            dim_dm = max([dim,dims]);
            xoff = (y_wfs.xpos)(cc)*4.848e-6*(y_dm.alt)(dd)/(y_tel.diam / y_geom.pupdiam);
            yoff = (y_wfs.ypos)(cc)*4.848e-6*(y_dm.alt)(dd)/(y_tel.diam / y_geom.pupdiam);
            xoff = float(xoff+(dim_dm-y_geom._n)/2);
            yoff = float(yoff+(dim_dm-y_geom._n)/2);
            sensors_addlayer,g_wfs,cc-1,y_dm(dd).type,(y_dm.alt)(dd),xoff,yoff;
          }
        }
      }
    }
  }
  
  if (y_target != []) {
    sizes = y_geom.pupdiam;
    sizes = sizes(-::y_target.ntargets-1);
    
    g_target = yoga_target(y_target.ntargets,*y_target.xpos,*y_target.ypos,*y_target.lambda,*y_target.mag,sizes,*y_geom._spupil);

    for (cc=1;cc<=y_target.ntargets;cc++) {
      if (y_atmos != []) {
        for (dd=1;dd<=y_atmos.nscreens;dd++) {
          xoff = (*y_target.xpos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
          yoff = (*y_target.ypos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
          xoff = float(xoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
          yoff = float(yoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
          // to take into account the difference in screen size
          // for target screen is of the same size as spupil
          // for wfs screen is of same size as mpupil
          pupdiff = (y_geom._n - y_geom.pupdiam)/2
            xoff += pupdiff;
          yoff += pupdiff;
          target_addlayer,g_target,cc-1,type,(*y_atmos.alt)(dd),xoff,yoff;
        }
      }
      if (y_dm != []) {
        for (dd=1;dd<=numberof(y_dm);dd++) {
          dims = y_dm(dd)._n2 - y_dm(dd)._n1 + 1;
          dim  = dimsof(*y_geom._mpupil)(2);
          dim_dm = max([dim,dims]);
          xoff = (*y_target.xpos)(cc)*4.848e-6*(y_dm.alt)(dd)/(y_tel.diam / y_geom.pupdiam);
          yoff = (*y_target.ypos)(cc)*4.848e-6*(y_dm.alt)(dd)/(y_tel.diam / y_geom.pupdiam);
          xoff = float(xoff+(dim_dm-y_geom._n)/2);
          yoff = float(yoff+(dim_dm-y_geom._n)/2);
          pupdiff = (y_geom._n - y_geom.pupdiam)/2;
          xoff += pupdiff;
          yoff += pupdiff;
          if (y_dm(dd).type=="kl") {
            xoff += 2;
            yoff += 2;
          }
          target_addlayer,g_target,cc-1,y_dm(dd).type,(y_dm.alt)(dd),xoff,yoff;
        }
      }
      target_init_strehlmeter,g_target,cc-1;
    }
  }
}

func dm_init(void)
/* DOCUMENT dm_init
   dm_init
     
   inits a yDM object on the gpu
   no input parameters
   
   requires 3 externals : 
   y_ geom   : a y_struct for the geometry
   y_ tel    : a y_struct for the telescope
   y_dm    : a y_struct for the dm
   creates 1 external :
   g_dm    : a yDMs object on the gpu

  SEE ALSO:
 */
{
  extern y_dm,g_dm;
  
  if (y_dm != []) {
    g_dm = yoga_dms(numberof(y_dm));
    for (n=1;n<=numberof(y_dm);n++) {
      if (y_dm(n).pupoffset!=[])                                          \
        y_dm(n)._puppixoffset = long(y_dm(n).pupoffset/y_tel.diam*y_geom.pupdiam);
      
      if (y_dm(n).type == "pzt") {
      // find out the support dimension for the given mirror.
        patchDiam = long(y_geom.pupdiam+2*max(abs([y_wfs.xpos,y_wfs.ypos]))*
                         4.848e-6*abs(y_dm(n).alt)/(y_tel.diam/y_geom.pupdiam));
        patchDiam;
        y_dm(n)._pitch = patchDiam / (y_dm(n).nact -1);

        extent = y_dm(n)._pitch*(y_dm(n).nact+3.); // + 1.5 pitch each side
        y_dm(n)._n1 = long(clip(floor(y_geom.cent-extent/2.),1,));
        y_dm(n)._n2 = long(clip(ceil(y_geom.cent+extent/2.),,y_geom.ssize));
      } else {  // we are dealing with a TT, or kl
        //y_dm(n)._n1 = 1;
        //y_dm(n)._n2 = y_geom.ssize;
        extent=y_geom.pupdiam+16;
        y_dm(n)._n1 = long(clip(floor(y_geom.cent-extent/2.),1,));
        y_dm(n)._n2 = long(clip(ceil(y_geom.cent+extent/2.),,y_geom.ssize));
      }
      
      if (y_dm(n).alt == 0) {
        extent=y_geom.pupdiam+64;
        y_dm(n)._n1 = long(clip(floor(y_geom.cent-extent/2.),1,));
        y_dm(n)._n2 = long(clip(ceil(y_geom.cent+extent/2.),,y_geom.ssize));
      }

      if (y_dm(n).type == "pzt") {
        make_pzt_dm, n;
        
        dims       = long(y_dm(n)._n2-y_dm(n)._n1+1);
        dim       = dimsof(*y_geom._mpupil)(2);
        if (dims >= dim) dim = dims;
        
        ninflu    = long(y_dm(n)._ntotact);
        influsize = long(y_dm(n)._influsize);
        ninflupos = long(numberof(*y_dm(n)._influpos));
        n_npts    = long(numberof(*y_dm(n)._ninflu));
        
        yoga_addpzt,g_dm,float(y_dm(n).alt),dim,ninflu,influsize,ninflupos,n_npts,float(y_dm(n).push4imat);
        
        yoga_loadpzt,g_dm,float(y_dm(n).alt),float(*y_dm(n)._influ),
          int(*y_dm(n)._influpos),int(*y_dm(n)._ninflu),int(*y_dm(n)._influstart),
          int(*y_dm(n)._i1),int(*y_dm(n)._j1); 
       }
      if (y_dm(n).type == "tt") {
        //dim       = dimsof(*y_geom._mpupil)(2);
        dim       = long(y_dm(n)._n2-y_dm(n)._n1+1);
        make_tiptilt_dm, n;
        
        yoga_addtt,g_dm,float(y_dm(n).alt),dim,float(y_dm(n).push4imat);
        
        yoga_loadtt,g_dm,float(y_dm(n).alt),float(*y_dm(n)._influ);
      }
      if (y_dm(n).type == "kl") {
        dim       = long(y_dm(n)._n2-y_dm(n)._n1+1);
        make_kl_dm, n;
        
        dim       = long(y_dm(n)._n2-y_dm(n)._n1+1);
        ninflu    = long(y_dm(n).nkl);
        influsize = long((*y_dm(n)._klbas).ncp);
        nr        = long((*y_dm(n)._klbas).nr);
        np        = long((*y_dm(n)._klbas).np);

        yoga_addkl,g_dm,float(y_dm(n).alt),dim,ninflu,influsize,nr,np,float(y_dm(n).push4imat);

        yoga_loadkl,g_dm,float(y_dm(n).alt),float(*(*y_dm(n)._klbas).rabas)(*),
          float(*(*y_dm(n)._klbas).azbas)(*),int(*(*y_dm(n)._klbas).ord),float(*(*y_dm(n)._klbas).cr)(*),
          float(*(*y_dm(n)._klbas).cp)(*);
        //error;
        /*
        // verif :
res1 = pol2car(*y_dm(n)._klbas,gkl_sfi(*y_dm(n)._klbas, 1));
res2 = yoga_getkl(g_dm,0.,1);
        */
      }
    }
  }
}

func rtc_init(clean=)
/* DOCUMENT rtc_init
   rtc_init[,clean=]
     
   inits a yRTC object on the gpu
   use clean=0 to avoid measuring the imat
   
   requires 2 externals : 
   y_wfs    : a y_struct for the sensors
   y_rtc    : a y_struct for the rtc
   creates 1 external :
   g_rtc    : a yRTC object on the gpu

  SEE ALSO:
 */
{
  extern g_rtc;
  g_rtc = yoga_rtc(activeDevice());
  if (y_rtc != []) {
    centroiders = *y_rtc.centroiders;
    if (y_wfs != []) {
      for (i=1;i<=numberof(centroiders);i++) {
        
        nwfs = centroiders(i).nwfs;

        if (y_wfs(nwfs).type == "sh") {
          if (centroiders(i).type != "corr") {
            s_offset = (y_wfs(centroiders(i).nwfs).npix/2.+0.5);
          } else {
            if (centroiders(i).type_fct == "model") {
              if (y_wfs(nwfs).npix %2 == 0) {
                s_offset = (y_wfs(centroiders(i).nwfs).npix/2.+0.5);
              }
              else {
                s_offset = (y_wfs(centroiders(i).nwfs).npix/2.+1.);
              }
            } else {
              s_offset = (y_wfs(centroiders(i).nwfs).npix/2.+1);
            }
          }
          s_scale = y_wfs(centroiders(i).nwfs).pixsize;
        }

        if (y_wfs(nwfs).type == "pyr") {
          s_offset = s_scale = 0.0f
        }

        rtc_addcentro,g_rtc,nwfs-1,y_wfs(nwfs)._nvalid,centroiders(i).type,s_offset,s_scale;
        
        sensors_initbcube,g_wfs,nwfs-1,g_rtc,i-1;

        if (y_wfs(nwfs).type == "sh") {
          if (centroiders(i).type == "tcog") {
            rtc_setthresh,g_rtc,i-1,centroiders(i).thresh;
          }
          
          if (centroiders(i).type == "nmax") {
            rtc_setnmax,g_rtc,i-1,centroiders(i).nmax;
          }
          
          if (centroiders(i).type_fct == "model") {
            if (y_wfs(nwfs).gsalt > 0) {
              profilename = "allProfileNa_withAltitude_1Gaussian.fits";
              prof=fits_read(YOGA_AO_SAVEPATH+profilename);
              h    = prof(,1);
              prof = prof(,2:)(,avg);
              make_lgs_prof1d,nwfs,prof,h,y_wfs(nwfs).beamsize,center="image";
              tmp=(*y_wfs(nwfs)._lgskern);
              tmp2 = makegaussian(dimsof(tmp)(2),2*y_wfs(nwfs)._nrebin);
              tmp3 = roll( (fft(fft(tmp)*(fft(tmp2)),-1) ).re);
              
              offset = (y_wfs(nwfs)._Ntot-y_wfs(nwfs)._nrebin*y_wfs(nwfs).npix)/2;
              rr = offset+1 : offset + y_wfs(nwfs)._nrebin*y_wfs(nwfs).npix;
              
              centroiders(i).weights = &float(tmp3( rr, rr ,)(cum,cum,)
                                              (::y_wfs(nwfs)._nrebin,::y_wfs(nwfs)._nrebin,)(dif,dif,));
              /*
                nphase = sensors_getdata(g_wfs,nwfs-1,"phase")*0.0f;
                sensors_setphase,g_wfs,nwfs-1,nphase;
                sensors_compimg,g_wfs,nwfs-1;
                ref_fctn = sensors_getdata(g_wfs,nwfs-1,"bincube");
                //centroiders(i).weights = &float(ref_fctn);
                */
            } else {
              r0 = y_atmos.r0 / ((y_wfs(nwfs).lambda/0.5)^(6./5));
              seeing = RASC * (y_wfs(nwfs).lambda * 1.e-6) / r0;   
              npix = seeing / y_wfs(nwfs).pixsize;
              centroiders(i).width = npix;
            }
          }
          
          if (centroiders(i).type_fct == "gauss") {
            if (y_wfs(nwfs).npix % 2 == 1)
              centroiders(i).weights = &float(makegaussian(y_wfs(nwfs).npix,centroiders(i).width,
                                                           xc=int(y_wfs(nwfs).npix/2)+1,
                                                           yc=int(y_wfs(nwfs).npix/2)+1)(,,-:1:y_wfs(nwfs)._nvalid));
            else {
              if (centroiders(i).type == "corr")
                centroiders(i).weights = &float(makegaussian(y_wfs(nwfs).npix,centroiders(i).width,
                                                             xc=int(y_wfs(nwfs).npix/2),
                                                             yc=int(y_wfs(nwfs).npix/2))(,,-:1:y_wfs(nwfs)._nvalid));
              else
                centroiders(i).weights = &float(makegaussian(y_wfs(nwfs).npix,centroiders(i).width,
                                                             xc=int(y_wfs(nwfs).npix/2)+0.5,
                                                             yc=int(y_wfs(nwfs).npix/2)+0.5)(,,-:1:y_wfs(nwfs)._nvalid));
            }              
          }
          
          if (centroiders(i).type == "wcog") {
            sensors_initweights,g_wfs,nwfs-1,g_rtc,i-1,*(centroiders(i).weights);
          }
          
          if (centroiders(i).type == "corr") {
            aa = array(0.0f,2*y_wfs(nwfs).npix,2*y_wfs(nwfs).npix);
            aa(1:y_wfs(nwfs).npix,1:y_wfs(nwfs).npix) = 1.0;
            corrnorm = float(roll(fft(abs(fft(aa))^2,-1).re)*0.+1.0f)(2:,2:);//;
            centroiders(i).sizex = 3;
            centroiders(i).sizey = 3;
            centroiders(i).interpmat = &float(create_interp_mat(centroiders(i).sizex,centroiders(i).sizey));
            sensors_initcorr,g_wfs,nwfs-1,g_rtc,i-1,*(centroiders(i).weights),corrnorm,
              centroiders(i).sizex,centroiders(i).sizey,*(centroiders(i).interpmat);
          }
        }
      }
      
      y_rtc.centroiders = &centroiders;
      
      controlers = *y_rtc.controlers;
      
      if ((y_wfs != []) && (y_dm != [])) {
        for (i=1;i<=numberof(controlers);i++) {
          nwfs = *controlers(i).nwfs;
          if (numberof(y_wfs) == 1) nwfs = nwfs(1); // fixing a bug ... still not understood
          ndms = *controlers(i).ndm;
          controlers(i).nvalid = &(y_wfs(nwfs)._nvalid);
          controlers(i).nactu  = &(y_dm(ndms)._ntotact);
          rtc_addcontrol,g_rtc,sum(y_dm(ndms)._ntotact),controlers(i).delay,controlers(i).type;
          write,"doing imat and filtering unseen actuators";
          imat_init,i,clean=clean;
          write,"done";
          cmat_init,i,clean=clean;
          rtc_setgain,g_rtc,0,controlers(i).gain;
          mgain = array(1.0f,(y_dm._ntotact)(sum));
          // filtering tilt ...
          //mgain(-1:0) = 0.0f;
          rtc_loadmgain,g_rtc,0,mgain;
        }
      }
    }
  }
}

