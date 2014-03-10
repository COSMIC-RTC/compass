func imat_init(ncontrol,clean=)
{
  extern y_rtc,g_dm;

  dirsave = YOGA_AO_SAVEPATH+"mat/";
  mkdirp,dirsave;

  if (clean == []) clean = 1;

  if (simul_name == []) imat_clean = 1;
  else
    imat_clean = ((!fileExist(swrite(format=dirsave+"imat-%d-%s.fits",ncontrol,simul_name))) || clean);

  if (imat_clean) {
    // first check if wfs is using lgs
    // if so, load new lgs spot, just for imat
    for (cc=1;cc<=numberof(y_wfs);cc++) {
      if (y_wfs(cc).gsalt > 0) {
        // lgs mode requested
        profilename = "allProfileNa_withAltitude_1Gaussian.fits";
        prof=fits_read(YOGA_AO_SAVEPATH+profilename);
        h    = prof(,1);
        prof = prof(,2:)(,avg);
        prep_lgs_prof,cc,prof,h,y_wfs(cc).beamsize,imat=1;
        //y_wfs(cc)._altna  = &h;
        //y_wfs(cc)._profna = &prof;
      }
    }
    
    rtc_doimat,g_rtc,ncontrol-1,g_wfs,g_dm;
    imat = rtc_getimat(g_rtc,ncontrol-1);

  } else {
    if (simul_name != [])
      imat = fits_read(swrite(format=dirsave+"imat-%d-%s.fits",ncontrol,simul_name));
  }

  correct_dm,imat;
  
  /*
  g_dm = 0;
  g_dm = yoga_dms(numberof(y_dm));
    
  for (nm=1;nm<=numberof(y_dm);nm++) {
    // filter actuators only in stackarray mirrors:
    inds = 1;
    if (y_dm(nm).type == "pzt") {
      dmx = *y_dm(nm)._xpos;
      dmy = *y_dm(nm)._ypos;
      dmi1 = *y_dm(nm)._i1;
      dmj1 = *y_dm(nm)._j1;
        
      if (imat_clean) {
        tmp = resp(inds:inds+y_dm(nm)._ntotact-1);
        ok = where(tmp >  y_dm(nm).thresh*max(tmp));
        nok= where(tmp <= y_dm(nm).thresh*max(tmp));
        if (simul_name != []) {
          fits_write,swrite(format=dirsave+"pztok-%d-%s.fits",nm,simul_name),ok,overwrite=1;
          fits_write,swrite(format=dirsave+"pztnok-%d-%s.fits",nm,simul_name),nok,overwrite=1;
        }
      } else {
        ok = fits_read(swrite(format=dirsave+"pztok-%d-%s.fits",nm,simul_name));
        nok= fits_read(swrite(format=dirsave+"pztnok-%d-%s.fits",nm,simul_name));
      }
      
      y_dm(nm)._xpos    = &(dmx(ok));
      y_dm(nm)._ypos    = &(dmy(ok));
      y_dm(nm)._i1      = &(int(dmi1(ok)));
      y_dm(nm)._j1      = &(int(dmj1(ok)));
      y_dm(nm)._influ   = &((*(y_dm(nm)._influ))(,,ok));
      y_dm(nm)._ntotact = (dimsof(*(y_dm(nm)._influ)))(4);
        
      comp_dmgeom,nm;
      
      dims      = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
      dim       = dimsof(*y_geom._mpupil)(2);
      if (dims > dim) dim = dims;
      
      ninflu    = long(y_dm(nm)._ntotact);
      influsize = long(y_dm(nm)._influsize);
      ninflupos = long(numberof(*y_dm(nm)._influpos));
      n_npts    = long(numberof(*y_dm(nm)._ninflu));
      yoga_addpzt,g_dm,float(y_dm(nm).alt),dim,ninflu,influsize,ninflupos,n_npts,float(y_dm(nm).push4imat);
      
      yoga_loadpzt,g_dm,float(y_dm(nm).alt),float(*y_dm(nm)._influ),
        int(*y_dm(nm)._influpos),int(*y_dm(nm)._ninflu),int(*y_dm(nm)._influstart),
        int(*y_dm(nm)._i1),int(*y_dm(nm)._j1);
      
    } else if (y_dm(nm).type == "tt") {
        
      dim       = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
      yoga_addtt,g_dm,float(y_dm(nm).alt),dim,float(y_dm(nm).push4imat);
      
      yoga_loadtt,g_dm,float(y_dm(nm).alt),float(*y_dm(nm)._influ);
    } else if (y_dm(nm).type == "kl") {
      
      dim       = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
      ninflu    = long(y_dm(nm).nkl);
      influsize = long((*y_dm(nm)._klbas).ncp);
      nr        = long((*y_dm(nm)._klbas).nr);
      np        = long((*y_dm(nm)._klbas).np);
      
      yoga_addkl,g_dm,float(y_dm(nm).alt),dim,ninflu,influsize,nr,np,float(y_dm(nm).push4imat);
      
      yoga_loadkl,g_dm,float(y_dm(nm).alt),float(*(*y_dm(nm)._klbas).rabas)(*),
        float(*(*y_dm(nm)._klbas).azbas)(*),int(*(*y_dm(nm)._klbas).ord),float(*(*y_dm(nm)._klbas).cr)(*),
        float(*(*y_dm(nm)._klbas).cp)(*);
    }
    
    inds += y_dm(nm)._ntotact;
  }
  */
  ndms = *controllers(ncontrol).ndm;
  controllers(ncontrol).nactu  = &(y_dm(ndms)._ntotact);
  
  rtc_rmcontrol,g_rtc;
  
  rtc_addcontrol,g_rtc,sum(y_dm(ndms)._ntotact),controllers(ncontrol).delay,controllers(ncontrol).type;
  
  if (imat_clean) {
    tic;
    rtc_doimat,g_rtc,ncontrol-1,g_wfs,g_dm;
    write,format = "imat time : %f\n",tac();
    imat = rtc_getimat(g_rtc,ncontrol-1);
    if (simul_name != []) {
      fits_write,swrite(format=dirsave+"imat-%d-%s.fits",ncontrol,simul_name),imat,overwrite=1;
    }
  } else
    rtc_setimat,g_rtc,ncontrol-1,imat;
  
  controllers = *y_rtc.controllers;
  controllers(ncontrol).imat = &imat;
  y_rtc.controllers = &controllers;
  
  // now restore original profile in lgs spots
  for (cc=1;cc<=numberof(y_wfs);cc++) {
    if (y_wfs(cc).gsalt > 0) {
      // lgs mode requested
      h    = *y_wfs(cc)._altna;
      prof = *y_wfs(cc)._profna;
      prep_lgs_prof,cc,prof,h,y_wfs(cc).beamsize;
    }
  }
}

func cmat_init(ncontrol,clean=)
{
  extern y_rtc;

  dirsave = YOGA_AO_SAVEPATH+"mat/";
  mkdirp,dirsave;

  if (clean == []) clean = 1;

  if (simul_name == []) cmat_clean = 1;
  else
    cmat_clean = ((!fileExist(swrite(format=dirsave+"eigenv-%d-%s.fits",ncontrol,simul_name))) || clean);

  if (((*y_rtc.controllers(ncontrol)).type)(1) == "ls"){
  if (cmat_clean) {
    write,"doing svd";
    tic;
    rtc_imatsvd,g_rtc,ncontrol-1;
    write,format="svd time %f\n",tac();
    eigenv = controller_getdata(g_rtc,ncontrol-1,"eigenvals");
    if (simul_name != []) {
      U = controller_getdata(g_rtc,ncontrol-1,"U");
      fits_write,swrite(format=dirsave+"eigenv-%d-%s.fits",ncontrol,simul_name),eigenv,overwrite=1;
      fits_write,swrite(format=dirsave+"U-%d-%s.fits",ncontrol,simul_name),U,overwrite=1;
    }
  } else {
    eigenv  = fits_read(swrite(format=dirsave+"eigenv-%d-%s.fits",ncontrol,simul_name));
    U = fits_read(swrite(format=dirsave+"U-%d-%s.fits",ncontrol,simul_name));
    controller_setdata,g_rtc,ncontrol-1,"eigenvals",eigenv;
    controller_setdata,g_rtc,ncontrol-1,"U",U;
  }

  imat = rtc_getimat(g_rtc,ncontrol-1);

  maxcond = (*y_rtc.controllers)(ncontrol).maxcond;
  mfilt = where((eigenv/eigenv(-2)) < 1./maxcond);
  //mfilt = where(1./(eigenv/eigenv(3)) > maxcond);
  //nfilt = numberof(mfilt)+2;
  nfilt = numberof(mfilt);
  
  if ( (wfs_disp!=[]) && (numberof(*wfs_disp._winits) > 0)) {
    if ((*wfs_disp._winits)(5)) {
    window,(*wfs_disp._wins)(5);fma;logxy,0,1;
    plg, eigenv(::-1), marks=0;
    plmk, eigenv(::-1), msize = 0.3, marker=4;
    x0 = dimsof(imat)(3) - nfilt + 0.5;
    pldj, x0 ,min(eigenv), x0, max(eigenv), color="red";
    }
  }
    write,"building cmat";
    tic;
    rtc_buildcmat,g_rtc,ncontrol-1,nfilt;
  
    write,format="cmat time %f\n",tac();
  }
  if (((*y_rtc.controllers(ncontrol)).type)(1) == "mv"){
    rtc_buildcmatmv,g_rtc,ncontrol-1,y_dm(1).type;
  }

  cmat = rtc_getcmat(g_rtc,ncontrol-1);

  controllers = *y_rtc.controllers;
  controllers(ncontrol).cmat = &float(cmat);
  y_rtc.controllers = &controllers;
}


func cmat_update(ncontrol,maxcond)
{
  extern y_rtc;
//error;
  imat = rtc_getimat(g_rtc,ncontrol);

  eigenv = controller_getdata(g_rtc,ncontrol,"eigenvals");
  
  (*y_rtc.controllers)(ncontrol).maxcond = maxcond;
  mfilt = where((eigenv/eigenv(-2)) < 1./maxcond);
  //mfilt = where(1./(eigenv/eigenv(1)) > maxcond);
  nfilt = numberof(mfilt);
  write,format="nb modes filtered : %d",nfilt;
  if (numberof(*wfs_disp._winits) > 0) {
    if ((*wfs_disp._winits)(5)) {
    window,(*wfs_disp._wins)(5);fma;logxy,0,1;
    plg, eigenv(::-1), marks=0;
    plmk, eigenv(::-1), msize = 0.3, marker=4;
    x0 = dimsof(imat)(3) - nfilt + 0.5;
    pldj, x0 ,min(eigenv), x0, max(eigenv), color="red";
    }
  }
  write,"building cmat";
  tic;
  rtc_buildcmat,g_rtc,ncontrol,nfilt;
  write,format="cmat time %f\n",tac();

  cmat = rtc_getcmat(g_rtc,ncontrol);
  
  controllers = *y_rtc.controllers;
  controllers(ncontrol).cmat = &float(cmat);
  y_rtc.controllers = &controllers;
}


func manual_imat(void)
{
  slps = sensors_getslopes(g_wfs,0);
  nslp = numberof(slps);
  imat_cpu = slps(,-);
  
  for (nm=1;nm<=numberof(y_dm);nm++)
    yoga_resetdm,g_dm,y_dm(nm).type,y_dm(nm).alt;
  
  for (nm=1;nm<=numberof(y_dm);nm++) {
    for (i=1;i<=y_dm(nm)._ntotact;i++) {
      
        com = array(0.,y_dm(nm)._ntotact);
        com(i)= float(y_dm(nm).push4imat);
        yoga_setcomm,g_dm,y_dm(nm).type,y_dm(nm).alt,com;
        yoga_shapedm,g_dm,y_dm(nm).type,y_dm(nm).alt;
            
        //yoga_oneactu,g_dm,y_dm(nm).type,y_dm(nm).alt,i-1,float(y_dm(1).push4imat);
      
      dm_shape = yoga_getdm(g_dm,y_dm(nm).type,y_dm(nm).alt);
      
      sensors_trace,g_wfs,0,"dm",g_dm,1;
      //mscreen = sensors_getdata(g_wfs,0,"phase");
      
      sensors_compimg,g_wfs,0;
      mimg = sensors_getimg(g_wfs,0);
      
      //slopes_geom,g_wfs,0,0;
      rtc_docentroids,g_rtc,g_wfs,0;
      slps = sensors_getslopes(g_wfs,0);
      grow,imat_cpu,slps/float(y_dm(1).push4imat);
      fma;limits;
      display_slopes,slps,1,"Phase Difference";
      yoga_resetdm,g_dm,y_dm(nm).type,y_dm(nm).alt;
      //pause,500;
    } 
  }
  return imat_cpu(,2:);
}

func manual_cmat(ncontrol,nfilt)
{
  imat = rtc_getimat(g_rtc,ncontrol-1);
  eigenvalues   = SVdec(float(imat),u,vt);
  modToAct    = float(transpose(vt));
  mesToMod    = float(transpose(u));
  neigen = numberof(eigenvalues);
  mev   = array(float,neigen,neigen);
  for (i=1;i<=neigen-nfilt;i++) {
    mev(i,i)=float(1./eigenvalues(i));
  }
  // Compute the Command matrix:
  cmat = (modToAct(,+)*mev(+,))(,+) * mesToMod(+,);

  return cmat;
}

func imat_geom(void)
{
  slps = sensors_getslopes(g_wfs,0);
  nslp = numberof(slps);
  imat_cpu = slps(,-);
  
  for (nm=1;nm<=numberof(y_dm);nm++)
    yoga_resetdm,g_dm,y_dm(nm).type,y_dm(nm).alt;
  
  for (nm=1;nm<=numberof(y_dm);nm++) {
    for (i=1;i<=y_dm(nm)._ntotact;i++) {
      /*
      com = array(0.,y_dm(nm)._ntotact);
      com(i)= float(y_dm(nm).push4imat);
      yoga_setcomm,g_dm,y_dm(nm).type,y_dm(nm).alt,com;
      yoga_shapedm,g_dm,y_dm(nm).type,y_dm(nm).alt;
      */
      yoga_oneactu,g_dm,y_dm(nm).type,y_dm(nm).alt,i-1,float(y_dm(nm).push4imat);
      //dm_shape = yoga_getdm(g_dm,y_dm(nm).type,y_dm(nm).alt);
      //window,0;pli,dm_shape;
      
      sensors_trace,g_wfs,0,"dm",g_dm,1;
      
      //mscreen = sensors_getdata(g_wfs,0,"phase");
      //window,1;pli,mscreen;

      slopes_geom,g_wfs,0,1;
      
      slps = sensors_getslopes(g_wfs,0);
      
      grow,imat_cpu,slps/float(y_dm(nm).push4imat);
      
      //fma;limits;
      //display_slopes,slps,1,"Phase Difference";
      yoga_resetdm,g_dm,y_dm(nm).type,y_dm(nm).alt;
      //pause,500;
    } 
  }
  return imat_cpu(,2:);
}

func create_interp_mat(dimx,dimy)
/* DOCUMENT create_interp_mat(dimx,dimy)
     
   SEE ALSO:
 */
{
  nn = max([dimx,dimy]);
  tmp = indices(nn)(1:dimx,1:dimy,);
  tmp(,,1) -= (dimx/2+1);
  tmp(,,2) -= (dimy/2+1);
  
  mymat = [tmp(,,1)(*)^2,tmp(,,2)(*)^2,tmp(,,1)(*)*tmp(,,2)(*),tmp(,,1)(*),
           tmp(,,2)(*),tmp(,,1)(*)*0+1];
  return transpose(LUsolve(mymat(+,)*mymat(+,))(+,)*mymat(,+));
  // the transpose is to comply with original convention in fitmax
}

func fitmax(im,xmin,xmax,ymin,ymax,interp_mat=,full=)
/* DOCUMENT fitmax(im,xmin,xmax,ymin,ymax,interp_mat=,full=)
     
   SEE ALSO:
 */
{
  z = im(xmin:xmax,ymin:ymax)(*);
  nx = xmax-xmin+1;
  ny = ymax-ymin+1;

  xm = (xmin+xmax);
  ym = (ymin+ymax);
  if (xm % 2 == 0) {xm /= 2;
  } else xm = xm/2+1;
  if (ym % 2 == 0) {ym /= 2;
  } else ym = ym/2+1;

  if (interp_mat == []) {interp_mat = create_interp_mat(nx,ny);
  } else {if ((dimsof(interp_mat) != [2,nx*ny,6])(sum) != 0) {error,"interp_mat has wrong dimensions";}}
    
  p = interp_mat(+,) * z(+); // coeffs de p(1).x^2+p(2).y^2+p(3).xy+p(4).x+p(5).y+p(6)
  A = p(1);
  B = p(2);
  C = p(3);
  denom = C^2.-4*A*B;
  if( denom==0 ) {
    x0 = y0 = 0.0;
  } else {
    x0 = (2.*B*p(4)-p(5)*C)/denom;
    y0 = (2.*A*p(5)-p(4)*C)/denom;
  }
  if( full==1 ) {
    D = p(6)-(x0*y0*C+A*x0^2.+B*y0^2.);
    valmax = D;   // top of the polynom
    if( (B-A)==0 ) {
      t=pi/4;
      if( C==0 ) t=0.;
    } else {
      t = atan(C/(B-A))/2;
    }
    AA = B*sin(t)^2-C*cos(t)*sin(t)+A*cos(t)^2;
    BB = A*sin(t)^2+C*cos(t)*sin(t)+B*cos(t)^2;
    fwhmx = 1.66*sqrt( -D/2./AA );
    fwhmy = 1.66*sqrt( -D/2./BB );
    
    return [x0+xm,y0+ym,fwhmx,fwhmy,valmax,t];   // t = angle
  } else {
    return [x0+xm,y0+ym];
  }
}

func correlfft(a,b)
/* DOCUMENT correlfft(a,b)
     
   SEE ALSO:
 */
{
  extern normalisation;
  n = dimsof(a)(2);
  aa = bb = array(0.0,2*n,2*n);
  aa(1:n,1:n) = a;
  bb(1:n,1:n) = b;
  tmp1 = roll(fft(fft(aa)*conj(fft(bb)),-1).re)(2:,2:);
  if( normalisation!="divide" )
    return tmp1;
  
  aa(1:n,1:n) = 1.0;
  tmp2 = roll(fft(abs(fft(aa))^2,-1).re)(2:,2:);
  return tmp1/tmp2;
  
}


func findMaxParaboloid(corrMap,sizex,sizey,&imat)
/* DOCUMENT findMaxParaboloid(corrMap,sizex,sizey,imat)
     
   SEE ALSO:
 */
{
  if (imat == []) imat = create_interp_mat(sizex,sizey);
  iPos = where2(corrMap==max(corrMap))(,1);
  i0 = iPos(1);
  j0 = iPos(2);
  //res = fitmax(corrMap,i0,j0,size=3);
  res = fitmax(corrMap,i0-long(sizex)/2,i0+long(sizex)/2,j0-long(sizex)/2,j0+long(sizey)/2,interp_mat=imat);
  n = (1+dimsof(corrMap)(2))/2;
  res += (0.5-n/2);     // 0.5-n+n/2
  return res;
}

func correct_dm(imat)
{
  extern g_dm,y_dm;

  if (simul_name == []) imat_clean = 1;
  
  g_dm = 0;
  g_dm = yoga_dms(numberof(y_dm));
    
  resp = sqrt((imat^2.)(sum,));

  for (nm=1;nm<=numberof(y_dm);nm++) {
    // filter actuators only in stackarray mirrors:
    inds = 1;
    if (y_dm(nm).type == "pzt") {
      dmx = *y_dm(nm)._xpos;
      dmy = *y_dm(nm)._ypos;
      dmi1 = *y_dm(nm)._i1;
      dmj1 = *y_dm(nm)._j1;
        
      if (imat_clean) {
        tmp = resp(inds:inds+y_dm(nm)._ntotact-1);
        ok = where(tmp >  y_dm(nm).thresh*max(tmp));
        nok= where(tmp <= y_dm(nm).thresh*max(tmp));
        if (simul_name != []) {
          fits_write,swrite(format=dirsave+"pztok-%d-%s.fits",nm,simul_name),ok,overwrite=1;
          fits_write,swrite(format=dirsave+"pztnok-%d-%s.fits",nm,simul_name),nok,overwrite=1;
        }
      } else {
        ok = fits_read(swrite(format=dirsave+"pztok-%d-%s.fits",nm,simul_name));
        nok= fits_read(swrite(format=dirsave+"pztnok-%d-%s.fits",nm,simul_name));
      }
      
      y_dm(nm)._xpos    = &(dmx(ok));
      y_dm(nm)._ypos    = &(dmy(ok));
      y_dm(nm)._i1      = &(int(dmi1(ok)));
      y_dm(nm)._j1      = &(int(dmj1(ok)));
      y_dm(nm)._influ   = &((*(y_dm(nm)._influ))(,,ok));
      y_dm(nm)._ntotact = (dimsof(*(y_dm(nm)._influ)))(4);
        
      comp_dmgeom,nm;
      
      dims      = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
      dim       = dimsof(*y_geom._mpupil)(2);
      if (dims > dim) dim = dims;
      
      ninflu    = long(y_dm(nm)._ntotact);
      influsize = long(y_dm(nm)._influsize);
      ninflupos = long(numberof(*y_dm(nm)._influpos));
      n_npts    = long(numberof(*y_dm(nm)._ninflu));
      yoga_addpzt,g_dm,float(y_dm(nm).alt),dim,ninflu,influsize,ninflupos,n_npts,float(y_dm(nm).push4imat);
      
      yoga_loadpzt,g_dm,float(y_dm(nm).alt),float(*y_dm(nm)._influ),
        int(*y_dm(nm)._influpos),int(*y_dm(nm)._ninflu),int(*y_dm(nm)._influstart),
        int(*y_dm(nm)._i1),int(*y_dm(nm)._j1);
      
    } else if (y_dm(nm).type == "tt") {
        
      dim       = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
      yoga_addtt,g_dm,float(y_dm(nm).alt),dim,float(y_dm(nm).push4imat);
      
      yoga_loadtt,g_dm,float(y_dm(nm).alt),float(*y_dm(nm)._influ);
    } else if (y_dm(nm).type == "kl") {
      
      dim       = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
      ninflu    = long(y_dm(nm).nkl);
      influsize = long((*y_dm(nm)._klbas).ncp);
      nr        = long((*y_dm(nm)._klbas).nr);
      np        = long((*y_dm(nm)._klbas).np);
      
      yoga_addkl,g_dm,float(y_dm(nm).alt),dim,ninflu,influsize,nr,np,float(y_dm(nm).push4imat);
      
      yoga_loadkl,g_dm,float(y_dm(nm).alt),float(*(*y_dm(nm)._klbas).rabas)(*),
        float(*(*y_dm(nm)._klbas).azbas)(*),int(*(*y_dm(nm)._klbas).ord),float(*(*y_dm(nm)._klbas).cr)(*),
        float(*(*y_dm(nm)._klbas).cp)(*);
    }
    
    inds += y_dm(nm)._ntotact;
  }


}

// Florian features
func doklbasis(g_dm,nkl,n)
{ extern y_dm, y_geom;
  p4f = 1.;

  if (y_dm(n).type!="kl"){
    //dim       = long(y_dm(n)._n2-y_dm(n)._n1+1);
  dim       = dimsof(*y_geom._mpupil)(2);
  cobs  = y_tel.cobs;
  cent  = y_geom.cent;
  psize = y_geom.pupdiam;
  patchDiam = long(y_geom.pupdiam+2*max(abs([y_wfs.xpos,y_wfs.ypos]))*
                   4.848e-6*abs(y_dm(n).alt)/psize);

  klbas = make_klbas(nkl,cobs,patchDiam,funct="kolmo");
      
   dim       = long(y_dm(n)._n2-y_dm(n)._n1+1);
   ninflu    = long(nkl);
   influsize = long((klbas).ncp);
   nr        = long((klbas).nr);
   np        = long((klbas).np);

   yoga_addkl,g_dm,float(y_dm(n).alt),dim,ninflu,influsize,nr,np,float(y_dm(n).push4imat)/p4f;

   yoga_loadkl,g_dm,float(y_dm(n).alt),float(*(klbas).rabas)(*),
          float(*(klbas).azbas)(*),int(*(klbas).ord),float(*(klbas).cr)(*),
          float(*(klbas).cp)(*);
  
   //error;
   // rtc_doimatkl,g_rtc,0,g_wfs,g_dm,long(y_dm(n).alt);
   //rtc_doimat,g_rtc,0,g_wfs,g_dm;
   //error;
  }
  // KL basis recovery
  K = yoga_getkl(g_dm,y_dm(n).alt,0);
  tmp = (dimsof(K)(2)-y_geom._n)/2;
  K = K(tmp+1:-tmp,tmp+1:-tmp)(*)(indx_valid);
  K = K(,-);
  for (i=1;i<=nkl-1;i++) grow,K,yoga_getkl(g_dm,y_dm(n).alt,i)(tmp+1:-tmp,tmp+1:-tmp)(*)(indx_valid);
  if(y_dm(n).type!="kl")
    yoga_rmkl(g_dm,y_dm(n).alt);
  return K;
}

func docovmat(g_rtc,g_atmos,g_dm,Nactu,nkl,N,ndms,sig=)
{ extern y_dm, y_geom;
  nkl = Nactu;
  tmp = (dimsof(*y_geom._ipupil)(2)-y_geom._n)/2;
  pup = (*y_geom._ipupil)(tmp+1:-tmp,tmp+1:-tmp);
  indx_valid = where(pup);
  //if(y_dm(N(1)).type != "kl")
  //  nkl +=2;
  write,"Computing projection for covariance matrix computation...";
  if(y_dm(N(1)).type == "kl")
    K = doklbasis(g_dm,nkl,N(1));
  else
    K = dopztbasis(g_dm,ndms,Nactu);
  error;
  // Projecteur
  tmp = array(float,nkl,nkl);
  tmp2 = array(float,dimsof(K)(3),dimsof(K)(2));
  d_K = yoga_obj(K);
  d_tmp = yoga_obj(tmp);
  d_tmp2 = yoga_obj(tmp2);

  yoga_mm,d_tmp,d_K,d_K,'t','n'; // Kt*K
  yoga_potri,d_tmp;       // (Kt*K)⁻¹
  yoga_mm,d_tmp2,d_tmp,d_K, 'n', 't'; // (Kt*K)⁻¹K
  P = d_tmp2();
  d_K = d_tmp = d_tmp2 = [];
  write,"done";
  alpha = array(float,nkl);
  alpha = alpha(,-);
  niter = 1000;
  //looptime = array(float,niter);
  //movetime = array(float,niter);
  //alphatime = array(float,niter);
  //timeloop = timemove = timealpha = 0.;
  for (cc=1;cc<=niter;cc++) {
    // tic(1);
    //tic(2);
    for (ii=1;ii<=100;ii++) {
      move_atmos,g_atmos;
    }
    //movetime(cc) = tac(2);
    sensors_trace,g_wfs,0,"atmos",g_atmos;
    res=sensors_getdata(g_wfs,0,"phase")(*)(indx_valid);
    res -= avg(res);
    // calcul du vecteur alpha courant
    //tic(3);
    alphav = P(,+)*res(+);
    grow, alpha,alphav;
    //alphatime(cc) = tac(3);
    write,format="\r Computing covariance matrix ... %1.0f",floor((cc/float(niter))*100);
    //looptime(cc) = tac(1);
  }
  //error;
  // Calcul de la matrice de covariance
  alpha = alpha(,2:);
  cov_matrix = (alpha(,+) * alpha(,+)) / dimsof(alpha)(3);
  //error;
  // Inversion
  //if(y_dm(N(1)).type != "kl"){
   cov_matrix = LUsolve(cov_matrix);
   /*
    d_cov = yoga_obj(cov_matrix);
    yoga_potri,d_cov;
    cov_matrix = d_cov();
    d_cov = [];
  */
    //}
    //else{
  //error;
  //cov_matrix = 1/diag(cov_matrix);
  //  covmat = array(float,Nactu,Nactu);
  // for(jj=1;jj<=dimsof(cov_matrix)(2);jj++) 
  //     covmat(jj,jj) = 1/cov_matrix(jj,jj);
  // cov_matrix = covmat;
  //}
  if(sig==[]) sig = 1e-3;
  //sig=0;
  cov_matrix = cov_matrix*sig;
  // covmat = cov_matrix;
  // covmat = array(float,Nactu);
  // tmp = 0;
  // for (j = 1 ; j<=ndms ; j++){
  //   if(y_dm(j).type == "pzt")
  //     covmat(tmp+1:tmp+y_dm(j)._ntotact) = cov_matrix(3:y_dm(j)._ntotact+2);
  //   else
  //     covmat(tmp+1:tmp+y_dm(j)._ntotact) = cov_matrix(:y_dm(j)._ntotact);
  //  tmp += y_dm(j)._ntotact;
  //}
  return cov_matrix;
}

func klbasis_init(ndms,nctrl,nctrltot){
  extern y_dm,y_geom,y_rtc,y_wfs;
  extern g_rtc,g_dm,g_wfs;
  p4f = 1.;
  controllers = *y_rtc.controllers;
  imat = rtc_getimat(g_rtc,nctrl-1);
  tmp = 0;
  klbasis = array(float,dimsof(imat));
  for(i=1 ; i<= ndms ; i++){
    Nactu = y_dm(i)._ntotact;
    if(y_dm(i).type == "kl"){
      klbasis(,tmp+1:tmp+Nactu) = imat(,tmp+1:tmp+Nactu);
    }
    else{
      if(y_dm(i).type == "tt"){
      // Adding a temporary KL DM for imatkl computation
      dim       = dimsof(*y_geom._mpupil)(2);
      cobs  = y_tel.cobs;
      cent  = y_geom.cent;
      psize = y_geom.pupdiam;
      patchDiam = long(y_geom.pupdiam+2*max(abs([y_wfs.xpos,y_wfs.ypos]))*
                   4.848e-6*abs(y_dm(i).alt)/psize);
      klbas = make_klbas(Nactu,cobs,patchDiam,funct="kolmo");
      dim       = long(y_dm(i)._n2-y_dm(i)._n1+1);
      ninflu    = long(y_dm(i)._ntotact);
      influsize = long((klbas).ncp);
      nr        = long((klbas).nr);
      np        = long((klbas).np);
      yoga_addkl,g_dm,float(y_dm(i).alt),dim,ninflu,influsize,nr,np,float(y_dm(i).push4imat)/p4f;
      yoga_loadkl,g_dm,float(y_dm(i).alt),float(*(klbas).rabas)(*),
          float(*(klbas).azbas)(*),int(*(klbas).ord),float(*(klbas).cr)(*),
          float(*(klbas).cp)(*);
      }
      if(y_dm(i).type == "pzt"){
      dim       = dimsof(*y_geom._mpupil)(2);
      cobs  = y_tel.cobs;
      cent  = y_geom.cent;
      psize = y_geom.pupdiam;
      patchDiam = long(y_geom.pupdiam+2*max(abs([y_wfs.xpos,y_wfs.ypos]))*
                   4.848e-6*abs(y_dm(i).alt)/psize);
      klbas = make_klbas(Nactu+2,cobs,patchDiam,funct="kolmo");
      dim       = long(y_dm(i)._n2-y_dm(i)._n1+1);
      ninflu    = long(y_dm(i)._ntotact+2);
      influsize = long((klbas).ncp);
      nr        = long((klbas).nr);
      np        = long((klbas).np);
      yoga_addkl,g_dm,float(y_dm(i).alt),dim,ninflu,influsize,nr,np,float(y_dm(i).push4imat)/100000.;
      yoga_loadkl,g_dm,float(y_dm(i).alt),float(*(klbas).rabas)(*),
          float(*(klbas).azbas)(*),int(*(klbas).ord),float(*(klbas).cr)(*),
          float(*(klbas).cp)(*);
      }
      // Adding a temporary sensor layer
      cc = 1; // Only one wfs supported...
      dims = y_dm(i)._n2 - y_dm(i)._n1 + 1;
      dim  = dimsof(*y_geom._mpupil)(2);
      dim_dm = max([dim,dims]);
      xoff = (y_wfs.xpos)(cc)*4.848e-6*(y_dm.alt)(i)/(y_tel.diam / y_geom.pupdiam);
      yoff = (y_wfs.ypos)(cc)*4.848e-6*(y_dm.alt)(i)/(y_tel.diam / y_geom.pupdiam);
      xoff = float(xoff+(dim_dm-y_geom._n)/2);
      yoff = float(yoff+(dim_dm-y_geom._n)/2);
      sensors_addlayer,g_wfs,cc-1,"kl",(y_dm.alt)(i),xoff,yoff;
      // Adding a temporary controller
      rtc_addcontrol,g_rtc,y_dm(i)._ntotact,controllers(nctrl).delay,controllers(nctrl).type;
      // Computing imat KL
      if(y_dm(i).type == "pzt"){
	rtc_doimatkl4pzt,g_rtc,nctrltot,g_wfs,g_dm;
      }
      else{
	   rtc_doimatkl,g_rtc,nctrltot,g_wfs,g_dm;
      }
      imatkl = rtc_getimat(g_rtc,nctrltot);
      // Removing all the temporary objects
      rtc_rmcontrol,g_rtc;
      sensors_rmlayer,g_wfs,cc-1,"kl",(y_dm.alt)(i);
      yoga_rmkl,g_dm,(y_dm.alt)(i);
      //error;
      
      klbasis(,tmp+1:tmp+Nactu) = imatkl;
    }
    tmp += Nactu;
  }
  //error;
  rtc_loadklbasis,g_rtc,nctrl-1,klbasis;
}

func dopztbasis(g_dm,ndm,Nactu)
{ extern y_dm, y_geom;
  tmp = (dimsof(*y_geom._ipupil)(2)-y_geom._n)/2;
  pup = (*y_geom._ipupil)(tmp+1:-tmp,tmp+1:-tmp);
  indx_valid = where(pup);
  IFtot = array(float,numberof(indx_valid),Nactu);
  ind = 0;
  for(j=1;j<=ndm;j++){
    C = unit(y_dm(j)._ntotact);
    IF = build_dm_gpu(j,C(,1));
    tmp = (dimsof(IF)(2)-y_geom._n)/2;
    IF = IF(tmp+1:-tmp,tmp+1:-tmp)(*)(indx_valid);
    IF = IF(,-);
    for (i=2;i<=y_dm(j)._ntotact;i++) grow,IF,build_dm_gpu(j,C(,i))(tmp+1:-tmp,tmp+1:-tmp)(*)(indx_valid);
    IFtot(,ind+1:ind+y_dm(j)._ntotact) = IF;
    ind += y_dm(j)._ntotact;
    if(i == y_dm(j)._ntotact) raz = build_dm_gpu(j,array(0.0f,y_dm(j)._ntotact));
  }
  return IFtot;
}
