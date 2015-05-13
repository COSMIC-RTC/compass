
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
        
  } else {
    if (simul_name != [])
      imat = fits_read(swrite(format=dirsave+"imat-%d-%s.fits",ncontrol,simul_name));
  }

  if (imat_clean) {
    write,format="%s", "doing imat... ";
    tic;
    rtc_doimat,g_rtc,ncontrol-1,g_dm;
    write,format = "done in : %f s\n",tac();  
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

func imat_init_old(ncontrol,clean=)
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
    
    //yoga_start_profile;
    rtc_doimat,g_rtc,ncontrol-1,g_dm;
    imat = rtc_getimat(g_rtc,ncontrol-1);
    //yoga_stop_profile;
    
    //imat = imat_geom(meth=0);
    //rtc_setimat,g_rtc,ncontrol-1,imat;
    
  } else {
    imat = fits_read(swrite(format=dirsave+"imat-%d-%s.fits",ncontrol,simul_name));
    rtc_setimat,g_rtc,ncontrol-1,imat;
  }
  
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

func cmat_init(ncontrol,clean=,method=)
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
      if (simul_name != [] ) {
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
    if (eigenv(1) < eigenv(0)) mfilt = where((eigenv/eigenv(-2)) < 1./maxcond);
    else mfilt = where(1./(eigenv/eigenv(3)) > maxcond);
    //nfilt = numberof(mfilt)+2;
    nfilt = numberof(mfilt);
    //error;
    if ( (wfs_disp!=[]) && (numberof(*wfs_disp._winits) > 0)) {
      if ((*wfs_disp._winits)(5)) {
        window,(*wfs_disp._wins)(5);fma;logxy,0,1;
        if (eigenv(1) < eigenv(0)) {
          plg, eigenv(::-1), marks=0;
          plmk, eigenv(::-1), msize = 0.3, marker=4;
        } else {
          plg, eigenv, marks=0;
          plmk, eigenv, msize = 0.3, marker=4;
        }
        x0 = dimsof(imat)(3) - nfilt + 0.5;
        pldj, x0 ,min(eigenv), x0, max(eigenv), color="red";
      }
    }
    write,"building cmat";
    write,"filtering ",nfilt," modes";
    tic;
    rtc_buildcmat,g_rtc,ncontrol-1,nfilt;
  
    write,format="cmat time %f\n",tac();
  }
  if (((*y_rtc.controllers(ncontrol)).type)(1) == "mv"){
    N = array(0.0f,2*sum(y_wfs(*y_controllers(ncontrol).nwfs)._nvalid));
    ind = 0;
    for(i=1 ; i<=numberof(*y_controllers(ncontrol).nwfs) ; i++){
      k = (*y_controllers(ncontrol).nwfs)(i);
      N(ind+1:ind+2*y_wfs(k)._nvalid) = noise_cov(k);
      ind += 2*y_wfs(k)._nvalid;
    }
    rtc_loadnoisemat,g_rtc,ncontrol-1,N;
    write,"Building cmat...";
    rtc_buildcmatmv,g_rtc,ncontrol-1,y_controllers(ncontrol-1).maxcond;
    if(y_controllers(ncontrol-1).TTcond == 0) y_controllers(ncontrol-1).TTcond = y_controllers(ncontrol-1).maxcond;
    rtc_filtercmatmv,g_rtc,ncontrol-1,y_controllers(ncontrol-1).TTcond;
  }

  cmat = rtc_getcmat(g_rtc,ncontrol-1);

  controllers = *y_rtc.controllers;
  controllers(ncontrol).cmat = &float(cmat);
  y_rtc.controllers = &controllers;
}

func manual_imat(void)
{
  slps = rtc_getcentroids(g_rtc,0, g_wfs, 0);
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
      
      //dm_shape = yoga_getdm(g_dm,y_dm(nm).type,y_dm(nm).alt);
      
      sensors_trace,g_wfs,0,"dm",g_dm,1;
      //mscreen = sensors_getdata(g_wfs,0,"phase");
      
      sensors_compimg,g_wfs,0;
      //mimg = sensors_getimg(g_wfs,0);
      
      //slopes_geom,g_wfs,0,0;
      rtc_docentroids,g_rtc;
      if(y_wfs(1).type != "pyr")
        slps = rtc_getcentroids(g_rtc,0, g_wfs, 0);
      else
        slps = sensors_getslopes(g_wfs,0);
      grow,imat_cpu,slps/float(y_dm(1).push4imat);
      //fma;limits;
      //plg,slps;	
      //display_slopes,slps*100.,1,"Phase Difference";
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

func imat_geom(ncontrol,meth=)
{
  if (meth == []) meth = 0;
  nwfs = numberof(*y_controllers(ncontrol).nwfs);
  ndm = numberof(*y_controllers(ncontrol).ndm);
  nslp = 0;
  for(i=1 ; i<=nwfs ; i++){
    wfs = ((*y_controllers(ncontrol).nwfs)(i) - 1)(1);
    slps = sensors_getslopes(g_wfs,wfs);
    nslp += numberof(slps);
  }
  slps = array(0.f,nslp);
  imat_cpu =slps(,-);
  for (nmc=1 ; nmc<=ndm ; nmc++) {
    nm = (*y_controllers(ncontrol).ndm)(nmc)(1);
    yoga_resetdm,g_dm,y_dm(nm).type,y_dm(nm).alt;
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
      slps = array(0.f, sum(y_wfs._nvalid)*2);
      nslps = 1;
      for(nw = 1 ; nw<=nwfs ; nw++){
        wfs = ((*y_controllers(ncontrol).nwfs)(nw) - 1)(1);
        sensors_trace,g_wfs,wfs,"dm",g_dm,1;
      
      
        //mscreen = sensors_getdata(g_wfs,0,"phase");
        //window,1;pli,mscreen;
        //hitReturn;

        slopes_geom,g_wfs,wfs,meth;
      
        slps(nslps:nslps+y_wfs(wfs+1)._nvalid*2-1)=sensors_getslopes(g_wfs,wfs);
        nslps+=y_wfs(wfs+1)._nvalid*2;
        
      }
      
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

func correct_dm(imat,ncontrol)
{
  extern g_dm,y_dm;

  dirsave = YOGA_AO_SAVEPATH+"mat/";
  if (simul_name == []) imat_clean = 1;
  /*
    g_dm = 0;
    g_dm = yoga_dms(numberof(y_dm));
  */
  for(i=1 ; i<= numberof(*y_controllers(ncontrol).ndm) ; i++){
    nm = ((*y_controllers(ncontrol).ndm)(i))(1);
    yoga_rmdm,g_dm,y_dm(nm).type,y_dm(nm).alt;
  }
  resp = sqrt((imat^2.)(sum,));

  inds = 1;
  ndm = numberof(*y_controllers(ncontrol).ndm);
  for (nmc=1 ; nmc<=ndm ; nmc++) {
    nm = (*y_controllers(ncontrol).ndm)(nmc)(1);
    nactu_nm = y_dm(nm)._ntotact;
    // filter actuators only in stackarray mirrors:
    
    if (y_dm(nm).type == "pzt") {
      dmx = *y_dm(nm)._xpos;
      dmy = *y_dm(nm)._ypos;
      dmi1 = *y_dm(nm)._i1;
      dmj1 = *y_dm(nm)._j1;
      
      if(!imat_clean && fileExist(swrite(format=dirsave+"pztok-%d-%s.fits",nm,simul_name))){
        ok = fits_read(swrite(format=dirsave+"pztok-%d-%s.fits",nm,simul_name));
        nok= fits_read(swrite(format=dirsave+"pztnok-%d-%s.fits",nm,simul_name));
      } else {
        tmp = resp(inds:inds+y_dm(nm)._ntotact-1);
        ok = where(tmp >  y_dm(nm).thresh*max(tmp));
        nok= where(tmp <= y_dm(nm).thresh*max(tmp));
        if (simul_name != []) {
          fits_write,swrite(format=dirsave+"pztok-%d-%s.fits",nm,simul_name),ok,overwrite=1;
          fits_write,swrite(format=dirsave+"pztnok-%d-%s.fits",nm,simul_name),nok,overwrite=1;
        }
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
        int(*y_dm(nm)._i1),int(*y_dm(nm)._j1),float(*y_dm(nm)._influkernel);

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
    
    inds += nactu_nm;
  }


}

func dopztbasis(g_dm,ndm,Nactu)
{ extern y_dm, y_geom, y_wfs;
  //tmp = (dimsof(*y_geom._ipupil)(2)-y_geom._n)/2;
  tmp = (dimsof(*y_geom._ipupil)(2)-(y_dm(ndm)._n2 - y_dm(ndm)._n1 +1))/2;
  if (y_dm(ndm).alt == 0.){
    pup = (*y_geom._ipupil)(tmp+1:-tmp,tmp+1:-tmp);
  }
  else {
    psize=y_tel.diam/y_geom.pupdiam;
    patchDiam = long(y_geom.pupdiam+2*max(abs([y_wfs.xpos,y_wfs.ypos]))*4.848e-6*abs(y_dm(ndm).alt)/psize);
    pup = float(make_pupil(y_geom.ssize,patchDiam,xc=y_geom.cent,yc=y_geom.cent,cobs=y_tel.cobs))(tmp+1:-tmp,tmp+1:-tmp);
  }
  indx_valid = where(pup);
  IFtot = array(float,numberof(indx_valid),Nactu);
  ind = 0;
  //for(j=1;j<=ndm;j++){
  j = ndm;
  C =float( unit(y_dm(j)._ntotact));
  IF = build_dm_gpu(j,C(,1));
  //tmp = (dimsof(IF)(2)-y_geom._n)/2;
  //window,0; fma; pli,IF*pup;
  //hitReturn;
  IF = IF(*)(indx_valid);
  IF = IF(,-);
  for (i=2;i<=y_dm(j)._ntotact;i++){ 
    tmp = build_dm_gpu(j,C(,i));
    grow,IF,tmp(*)(indx_valid);
    //window,0; fma; pli,tmp*pup;
    //hitReturn;
  }
  IFtot(,ind+1:ind+y_dm(j)._ntotact) = IF;
  ind += y_dm(j)._ntotact;
  if(i == y_dm(j)._ntotact) raz = build_dm_gpu(j,array(0.0f,y_dm(j)._ntotact));
  //}
  return IFtot;
}
func DDiago(mat,del,&s)
/*
  Double diagonalisation. En general, mat est la matrice de cov statistique,
  et del la matrice de covariance geometrique.

  En sortie, on a une base de modes telle que

  b(+,)*(del(,+)*b(+,))(+,) = identite (modes orthonormes)

  et

  b1(,+) * (mat(,+)*b1(,+))(+,) = matrice diagonale (modes non correles)
*/
{
  s = SVdec(del,mp);
  m = mp / sqrt(s)(-,);
  // m1 = LUsolve(m);
  m1 = transpose(mp * sqrt(s)(-,));
  cp = m1(,+) * (mat(,+) * m1(,+))(+,);
  s = SVdec(cp,a);
  b = m(,+) * a(+,);
  return b
    }

func geo_cov(mat)
// Performs (transpose(mat) * mat) where mat is the matrix containing the influence functions in colums
// Matrice creuse --> code à optimiser
{
  d_mat = yoga_obj(mat);
  geocov = array(0.0f,dimsof(mat)(3),dimsof(mat)(3));
  d_geo = yoga_obj(geocov);
  yoga_mm,d_geo,d_mat,d_mat,'t','n';
  geocov = d_geo();
  d_mat = d_geo = [];  
 
  return geocov;
}

func stat_cov(n,r0)
// Compute the statistic covariance matrix on the actuators
{
  extern y_dm, y_geom,y_atmos;

  if (y_dm(n).type == "pzt"){
    // Actuators positions
    x = *y_dm(n)._xpos;
    y = *y_dm(n)._ypos;
    
    patchDiam = y_tel.diam+2*max(abs([y_wfs(n).xpos,y_wfs(n).ypos]))*4.848e-6*abs(y_dm(n).alt);
    interactp = double(x(2) - x(1));
    interactm = patchDiam/(y_dm(n).nact-1);
    p2m = interactm/interactp;
    norm = -(p2m/r0)^(5./3)/2.;
    //norm = -(0.0269484/r0)^(5./3);
    //norm = -(p2m*patchDiam/(2*r0))^(5./3);
    //norm = 1.;
    //error;
    //norm *= 1/75000.;
  }
  else if (y_dm(n).type == "kl"){
    N = int(ceil(sqrt(y_dm(n).nkl)));
    pup = make_pupil(N,N-1,cobs=y_tel.cobs);
    ind_sub = where(pup);
    while (numberof(ind_sub) < y_dm(n).nkl){
      N+=1;
      pup = make_pupil(N,N-1,cobs=y_tel.cobs);
      ind_sub = where(pup);
    }
    x = span(-1,1,N)(,-:1:N);
    y = transpose(x)(*)(ind_sub);
    x = x(*)(ind_sub);
    patchDiam = y_tel.diam+2*max(abs([y_wfs(n).xpos,y_wfs(n).ypos]))*4.848e-6*abs(y_dm(n).alt);
    norm = (patchDiam/(2*r0))^(-5./3) / y_atmos.nscreens;
  }
  // Kolmogorov statistic
  m = 6.88 * abs(x(*)(,-)-x(*)(-,),y(*)(,-)-y(*)(-,))^(5./3) ;
  
  // Filtering piston
  //F = double(unit(dimsof(m)(3)) - array(1./(dimsof(m)(3)),dimsof(m)(3),dimsof(m)(3)));
  //m = (F(,+)*m(+,))(,+)*F(,+);

  if (y_dm(n).type == "kl")
    m = m(:y_dm(n).nkl,:y_dm(n).nkl);

  return norm*m;
}

func create_sigmaTur(n){
  //Some constants
  k1 = 0.1716613621245709486;
  k2 = 1.0056349179985892838;
  //L0 = (*y_atmos.L0)(1);
  L0 = 1.e5;
  r0 = y_atmos.r0;
  // Actuators positions
  x = *y_dm(n)._xpos;
  y = *y_dm(n)._ypos;
  
  r0 = y_atmos.r0 * (y_wfs(nm).lambda/0.5)^(6/5);
  L0 = (*y_atmos.L0)(1)*y_tel.diam/y_geom.pupdiam; //conversion du L0 de pixels en metres
    
  patchDiam = y_tel.diam+2*max(abs([y_wfs.xpos,y_wfs.ypos]))*4.848e-6*abs(y_dm(n).alt);
  interactp = double(x(2) - x(1));
  interactm = patchDiam/(y_dm(n).nact-1);
  p2m = interactm/interactp;

  r = abs(x(*)(,-)-x(*)(-,),y(*)(,-)-y(*)(-,)) * p2m;
  return (k1 * k2 * L0^(5./3)  -  0.5 * rodconan(r, L0) ) * r0^(-5./3);
  //return ( -  0.5 * rodconan(r, L0) ) * r0^(-5./3);
  
    
}
  
func create_dmo(nm,nw) {
  // coefficient permettant d''obtenir les meme pentes entre l''ASO COMPASS et celle obtenues par D_MO apres injection du meme ecran de phase.

  // D_Mo[arcsec/micron] = D_Mo[px/rad] * coeff
  coeff=y_wfs(nw).pixsize*2*pi/y_wfs(nw).lambda;


  nb_p   = (y_wfs(nw)._nvalid) * 2;
  nb_act = y_dm(nm)._ntotact;
  DMo = array(0.0f,nb_p,nb_act);

  offs = 0;
  incx = (*y_dm(nm)._i1)(2) - (*y_dm(nm)._i1)(1);
  tmpx = (*y_dm(nm)._i1-offs)/incx+1;
  tmpy = (*y_dm(nm)._j1-offs)/incx;
  mask = array(0.0f,y_dm(nm).nact,y_dm(nm).nact);
  mask(tmpx + tmpy * y_dm(nm).nact) = 1;
  masq_act = where(mask);
  
  nLenslet = y_wfs(nw).nxsub;
  iD = 1;
  for (j=1;j<=nLenslet;j++) {
    for (i=1;i<=nLenslet;i++) {
      if ((*y_wfs(1)._isvalid)(i,j) == 1) {
        MD_ssp = array(0.0f,nLenslet+1,nLenslet+1);
        // à GAUCHE
        MD_ssp(i,j) = -1;
        MD_ssp(i+1,j) = 1;
        // à DROITE 
        MD_ssp(i,j+1) = -1;
        MD_ssp(i+1,j+1) = 1;
        MD_ssp = (MD_ssp/2/pi)*coeff;
        // repositionnement en colonne ( pupille CARRÉE  nActuator*nActuator  )
        V_ssp = MD_ssp(*);
        // repositionnement en colonne ( pupille DISQUE nb_act )
        VP_ssp = V_ssp(masq_act);
        // Ligne de DMo --> pente en X
        DMo(iD,) = VP_ssp;
        iD ++;
      }
    }
  }
  for (j=1;j<=nLenslet;j++) {
    for (i=1;i<=nLenslet;i++) {
      if ((*y_wfs(1)._isvalid)(i,j) == 1) {
        MD_ssp = array(0.0f,nLenslet+1,nLenslet+1);
        // à GAUCHE
        MD_ssp(i,j) = -1;
        MD_ssp(i,j+1) = 1;
        // à DROITE 
        MD_ssp(i+1,j) = -1;
        MD_ssp(i+1,j+1) = 1;
        MD_ssp = (MD_ssp/2/pi)*coeff;
        // repositionnement en colonne ( pupille CARRÉE  nActuator*nActuator  )
        V_ssp = MD_ssp(*);
        // repositionnement en colonne ( pupille DISQUE nb_act )
        VP_ssp = V_ssp(masq_act);
        // Ligne de DMo --> pente en X
        DMo(iD,) = VP_ssp;
        iD ++;
      }
    }
  }
  /*DMo2 =  array(0.0f,nb_p,nb_act);
    DMo2(1:nb_p/2,) = DMo(nb_p/2+1:nb_p,);
    DMo2(nb_p/2+1:nb_p,) = DMo(1:nb_p/2,);
    return DMo2;*/
  return DMo;
}

func create_nact(nm) {
  nb_act = y_dm(nm)._ntotact;
  nact = array(0.0f,nb_act,nb_act);

  tmpx = *y_dm(nm)._i1;
  tmpy = *y_dm(nm)._j1;
  offs = ((y_dm(1)._n2-y_dm(1)._n1+1) - (max(tmpx) - min(tmpx)))/2 - min(tmpx);
  tmpx += offs+1;
  tmpy += offs+1;
  mask = yoga_getdm(g_dm,y_dm(nm).type,y_dm(nm).alt) * 0;
  mask(tmpx + (tmpy-1) * dimsof(mask)(2)) = 1;
  masq_act = where(mask);
  for (i=1;i<=y_dm(nm)._ntotact;i++) {
    com = array(0.,y_dm(nm)._ntotact);
    com(i)= float(y_dm(nm).push4imat);
    yoga_setcomm,g_dm,y_dm(nm).type,y_dm(nm).alt,com;
    yoga_shapedm,g_dm,y_dm(nm).type,y_dm(nm).alt;
    shape = yoga_getdm(g_dm,y_dm(nm).type,y_dm(nm).alt);
    nact(i,) = shape(*)(masq_act) / float(y_dm(nm).push4imat);
  }

  return nact;
}

func create_sigmav(SigmaTur, isZonal, ordreAR, atur, btur)
{ //atur et btur sont des vecteurs (sous-diagonales de A1)
  dims = dimsof(SigmaTur);
  if ( (dims(1)!=2) || (dims(2)!=dims(3))){
    write, "Incorrect dimensions of matrix SigmaTur in function create_sigmav"
      }

  nb_az = dims(2);
  if (ordreAR == 1)
    {
      (A1_Tur = array(structof(atur), nb_az, nb_az))(1:nb_az*nb_az:nb_az+1) = atur;
      A1_TurT=transpose(A1_Tur);
      SigmaV1 = A1_Tur(,+)*SigmaTur(+,);
      SigmaV = SigmaTur - (SigmaV1(,+)*A1_TurT(+,));
    }
  else if (ordreAR == 2)
    {
      if (isZonal == 0)
        {
          (A2_Tur = array(structof(atur), nb_az, nb_az))(1:nb_az*nb_az:nb_az+1) = atur;
          (B2_Tur = array(structof(btur), nb_az, nb_az))(1:nb_az*nb_az:nb_az+1) = btur;
          A2Tur_SigmaTur = A2_Tur(,+)*SigmaTur(+,);
          Sig1 = LUsolve(unit(width_of(SigmaTur)) - B2_Tur , A2Tur_SigmaTur);
     
          B2Tur_SigmaTur = B2_Tur(,+)*SigmaTur(+,);
          A2Tur_Sig1 = A2_Tur(,+)*Sig1(+,);
          B2Tur_Sig1 = B2_Tur(,+)*Sig1(+,);


          Sig2 = A2Tur_SigmaTur(,+)*A2_Tur(+,) + B2Tur_SigmaTur(,+)*B2_Tur(+,) + A2Tur_Sig1(,+)*B2_Tur(+,) + B2Tur_Sig1(,+)*A2_Tur(+,);
          if (isZonal == 0) Sig2 = (Sig2 + transpose(Sig2))/2;
          SigmaV = SigmaTur - Sig2;
        }
    }
  //window,1;pli, SigmaTur
  //window,2;pli, SigmaV;error;
  return SigmaV;
}

//include,"dphi_rico.i";
/*
  func dphi_lowpass(r,x0,L0,rmax) {
  return (r^(5./3.)) *  Ij0t83(r*(pi/x0),L0,rmax)*(2*(2*pi)^(8/3.)*0.0228956);
  }

  func Ij0t83(x,L0,rmax) {
  extern YLARGE, XLARGE;
  if( YLARGE==[] ) {
  write,"computing array of tabulated integral";
  n = 10000;
  XLARGE = span(0,rmax,n);
  dX = (XLARGE(0)-XLARGE(1))/(n-1);
  y = (XLARGE^2 + (1./L0)^2)^(-8./6.) * unMoinsJ0(XLARGE);
  YLARGE = y(zcen)(cum) * dX;
  }
  return interp(YLARGE,XLARGE,x);
  }

  func unMoinsJ0(x,L=)
  {
  if( numberof(dimsof(x))==1 ) {   // x is a scalar
  if( x<0.1 ) {
  // J0(x) = x^2/4 - x^4/64 + ...
  //       = (x/2)^2  ( 1 - (x/2)^2 / 4 + ... ) to minimize computation errors
  x22 = (x/2.)^2; 
  return (1-x22/4.)*x22;
  } else {
  // classique
  return (1-bessj0(x));
  }
  } else {  // x is a vector
  y = double(x);
  for(i=1; i<=numberof(x); i++)
  y(i) = unMoinsJ0(x(i));
  return y;
  }
  }
*/

func selectDMforLayers(nc,Nlayers,indLayers){

  pztDM = where(y_dm(*y_controllers(nc).ndm).type == "pzt");
  
  for(i = 1 ; i<= y_atmos.nscreens ; i++){
    alt_diff = y_dm(pztDM).alt - (*y_atmos.alt)(i);
    ind = where(abs(alt_diff) == min(abs(alt_diff)));
    if(numberof(ind)>1) ind = ind(1);
    indLayers(i) = (*y_controllers(nc).ndm)(ind)-1;
    Nlayers(ind) += 1;
  }
}

func mat_cphim_gpu(nc){
  // Positions des coins inferieurs gauches des ssp valides ramenees dans ipupil (origine au centre)
  s2ipup = (dimsof(*y_geom._ipupil)(2) - dimsof(*y_geom._spupil)(2))/2.;
  
  nvalid = sum(y_wfs._nvalid(*y_controllers(nc+1).nwfs));
  ind = 0;
  X = Y = array(0.0,nvalid); // Tableaux finaux
  for(k=1 ; k<= numberof(*y_controllers(nc+1).nwfs) ; k++){
    kk = (*y_controllers(nc+1).nwfs)(k); // On ne considère que les WFS vus par le RTC
    posx = *y_wfs(kk)._istart + s2ipup ; // Positions ramenées dans ipupil
    posx = (posx * (*y_wfs(kk)._isvalid))(*); // Positions des ssp valides uniquements
    posx = posx(where(posx!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1; // et centrées sur la pupille
    posy = *y_wfs(kk)._jstart + s2ipup ;
    posy = (transpose(posy * (*y_wfs(kk)._isvalid)))(*);
    posy = posy(where(posy!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1;

    // Conversion en metres
    sspDiam = posx(2) - posx(1);
    p2m = (y_tel.diam/y_wfs(kk).nxsub) / sspDiam;
    posx *= p2m;
    posy *= p2m;
    
    //Remplis les tableaux finaux
    X(ind+1:ind+y_wfs(kk)._nvalid) = posx;
    Y(ind+1:ind+y_wfs(kk)._nvalid) = posy;
    ind += y_wfs(kk)._nvalid;
  }
  nactu = 0;
  for(k=1 ; k<=numberof(*y_controllers(nc+1).ndm) ; k++){
    kk = (*y_controllers(nc+1).ndm)(k);
    if(y_dm(kk).type == "pzt")
      nactu += y_dm(kk)._ntotact;
  }

  Xactu = Yactu = array(0.0,nactu);
  k2 = array(0.0,numberof(where(y_dm(*y_controllers(nc+1).ndm).type == "pzt")));
  pitch = k2;
  ind = 0;
  indk = 1;
  for(k=1 ; k<=numberof(*y_controllers(nc+1).ndm) ; k++){
    kk = (*y_controllers(nc+1).ndm)(k);
    if(y_dm(kk).type == "pzt"){
      //p2m = (y_tel.diam/(y_dm(kk).nact - 1)) / y_dm(kk)._pitch;
      p2m = y_tel.diam/y_geom.pupdiam;
      // Position des actuateurs et origine ramenee au centre de ipupil
      actu_x = (*y_dm(kk)._xpos - dimsof(*y_geom._ipupil)(2)/2 )*p2m;
      actu_y = (*y_dm(kk)._ypos - dimsof(*y_geom._ipupil)(2)/2 )*p2m;
      pitch(indk) = actu_x(2) - actu_x(1);
      
      //Remplis les tableaux finaux
      Xactu(ind+1 : ind+y_dm(kk)._ntotact) = actu_x;
      Yactu(ind+1 : ind+y_dm(kk)._ntotact) = actu_y;
      ind += y_dm(kk)._ntotact;
      k2(indk) = y_wfs(1).lambda / 2. / pi / y_dm(kk).unitpervolt / (y_dm(kk)._pitch * p2m) * 1.058 ; // 1.058 hardcoded for debug... need to find where it comes from
      indk ++;
    }
  }
  pztDM = where(y_dm(*y_controllers(nc).ndm).type == "pzt");
  NlayersDM = array(0,numberof(pztDM)); // Number of layers per DM
  indlayersDM = array(0,y_atmos.nscreens); // Index of the DM for each layer
  selectDMforLayers,nc+1,NlayersDM,indlayersDM;
  FoV = y_wfs.pixsize * y_wfs.npix / RASC;
  
  L0_d = &(double(*y_atmos.L0));
  alphaX = alphaY = array(0.0,numberof(y_wfs.xpos));
  alphaX(:numberof(y_wfs.xpos)) = y_wfs.xpos/RASC;
  alphaY(:numberof(y_wfs.xpos)) = y_wfs.ypos/RASC;

  write,"Computing Cphim...";
  rtc_doCphim,g_rtc,nc,g_wfs,g_atmos,g_dm,*L0_d,float(*y_atmos.frac * (y_atmos.r0)^(-5./3)),alphaX,alphaY,X,Y,Xactu,Yactu,y_tel.diam,k2,NlayersDM,indlayersDM,FoV,pitch,y_dm(pztDM).alt;
  write,"Done";

  write,"Computing Cmm ...";
  rtc_doCmm,g_rtc,nc,g_wfs,g_atmos,double(y_tel.diam),double(y_tel.cobs),*L0_d, double(*y_atmos.frac * (y_atmos.r0)^(-5./3)),alphaX,alphaY;
  write,"Done";

  Nact = F = array(0.0f,numberof(Xactu),numberof(Xactu));
  ind = 0;
  for(k=1 ; k<=numberof(*y_controllers(nc+1).ndm) ; k++){
    kk = (*y_controllers(nc+1).ndm)(k);
    if(y_dm(kk).type == "pzt"){
      Nact(ind+1:ind+y_dm(kk)._ntotact,ind+1:ind+y_dm(kk)._ntotact) = create_nact_geom(kk);
      tmp = array(-1./y_dm(kk)._ntotact,y_dm(kk)._ntotact,y_dm(kk)._ntotact);
      for(j=1 ; j<= y_dm(kk)._ntotact ; j++)
        tmp(j,j) = 1 - 1/y_dm(kk)._ntotact;
      F(ind+1:ind+y_dm(kk)._ntotact,ind+1:ind+y_dm(kk)._ntotact) = tmp;
      ind += y_dm(kk)._ntotact;
    }
  }

  rtc_filterCphim,g_rtc,0,F,Nact;

}
func mat_cphim(void)
{
  // Positions des coins inferieurs gauches des ssp valides ramenees dans ipupil (origine au centre)
  s2ipup = (dimsof(*y_geom._ipupil)(2) - dimsof(*y_geom._spupil)(2))/2.;
  posx = *y_wfs(1)._istart + s2ipup ;
  posx = (posx * (*y_wfs(1)._isvalid))(*);
  posx = posx(where(posx!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1;
  posy = *y_wfs(1)._jstart + s2ipup ;
  posy = (transpose(posy * (*y_wfs(1)._isvalid)))(*);
  posy = posy(where(posy!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1;

  // Conversion en metres
  sspDiam = posx(2) - posx(1);
  p2m = (y_tel.diam/y_wfs(1).nxsub) / sspDiam;
  posx *= p2m;
  posy *= p2m;
  sspDiam = sspDiam*p2m;
  
  //Calcul de la position des points A,B (centre des bords gauche et droit des ssp)
  Ax = posx;
  Ay = posy + sspDiam/2.;
  Bx = posx + sspDiam;
  By = Ay;
  //Calcul de la position des points C,D (centre des bords bas et haut des ssp)
  Cx = posx + sspDiam/2.;
  Cy = posy;
  Dx = Cx;
  Dy = posy + sspDiam;

  // Position des actuateurs et origine ramenee au centre de ipupil
  actu_x = (*y_dm(1)._xpos - dimsof(*y_geom._ipupil)(2)/2 )*p2m;
  actu_y = (*y_dm(1)._ypos - dimsof(*y_geom._ipupil)(2)/2 )*p2m;

  //Position de l'actuateur de reference
  actu_x0 = actu_x(numberof(actu_x)/2 + 1) * 0;
  actu_y0 = actu_x(numberof(actu_y)/2 + 1) * 0; 
  //error;

  //Taille max du support
  dmax = 0.;
  maxalt = *y_atmos.alt(y_atmos.nscreens);
  for (cc=0 ; cc<=numberof(y_wfs) ; cc++) {
    tmp = abs(y_wfs(cc).xpos/RASC,y_wfs(cc).ypos/RASC);
    if (tmp > dmax) dmax = tmp;
  }
  rmax = dmax * 2 * maxalt + y_tel.diam;

  //error;
  cphim = array(0.0f,y_dm(1)._ntotact,2*y_wfs(1)._nvalid);
  seeing = y_wfs(1).lambda * 1e-6 / y_atmos.r0;
  for (cs=1 ; cs<= y_atmos.nscreens ; cs++){
    L0 = (*y_atmos.L0)(cs);
    delta_h = abs(y_dm(1).alt - (*y_atmos.alt)(cs));
    x0 = (actu_x(2) - actu_x(1)) * 0.5/*+ (delta_h * seeing / 2.)*/;
    k1 = RASC * y_wfs(1).lambda * 1.e-6 / 2. / pi / sspDiam;
    k2 = y_wfs(1).lambda / 2. / pi / y_dm(1).unitpervolt;
    k3 = 1.;
    //error;
    for(i = 1 ; i<=dimsof(cphim)(2) ; i++){
      //Covariance selon x
      for(j=1 ; j<=dimsof(cphim)(3)/2 ; j++){
        ca = sqrt((actu_x(i) - Ax(j))^2 + (actu_y(i) - Ay(j))^2);
        cb = sqrt((actu_x(i) - Bx(j))^2 + (actu_y(i) - By(j))^2);
        bo = sqrt((Bx(j) - actu_x0)^2 + (By(j) - actu_y0)^2);
        ao = sqrt((Ax(j) - actu_x0)^2 + (Ay(j) - actu_y0)^2);
      
        //cphim(i,j) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi_lowpass(ca,x0,L0,rmax) + dphi_lowpass(bo,x0,L0,rmax) - dphi_lowpass(cb,x0,L0,rmax) - dphi_lowpass(ao,x0,L0,rmax));
        cphim(i,j) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi_lowpass(ca,x0) + dphi_lowpass(bo,x0) - dphi_lowpass(cb,x0) - dphi_lowpass(ao,x0));
        //cphim(i,j) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi_lowpass(ca,x0)  - dphi_lowpass(cb,x0));
      }
      j=1;
      //Covariance selon y
      for(jj=dimsof(cphim)(3)/2 + 1 ; jj<=dimsof(cphim)(3) ; jj++){
        cc = sqrt((actu_x(i) - Cx(j))^2 + (actu_y(i) - Cy(j))^2);
        cd = sqrt((actu_x(i) - Dx(j))^2 + (actu_y(i) - Dy(j))^2);
        Do = sqrt((Dx(j) - actu_x0)^2 + (Dy(j) - actu_y0)^2);
        co = sqrt((Cx(j) - actu_x0)^2 + (Cy(j) - actu_y0)^2);
        j++;
      
        //cphim(i,jj) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi_lowpass(cc,x0,L0,rmax) + dphi_lowpass(Do,x0,L0,rmax) - dphi_lowpass(cd,x0,L0,rmax) - dphi_lowpass(co,x0,L0,rmax));
        cphim(i,jj) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi_lowpass(cc,x0) + dphi_lowpass(Do,x0) - dphi_lowpass(cd,x0) - dphi_lowpass(co,x0));
        //cphim(i,jj) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi_lowpass(cc,x0) - dphi_lowpass(cd,x0));
      }
    }
  }
  
  return k1*k2*k3*cphim;
}

func mat_cmm(void)
{
  // Positions des coins inferieurs gauches des ssp valides ramenees dans ipupil (origine au centre)
  s2ipup = (dimsof(*y_geom._ipupil)(2) - dimsof(*y_geom._spupil)(2))/2.;
  posx = *y_wfs(1)._istart + s2ipup ;
  posx = (posx * (*y_wfs(1)._isvalid))(*);
  posx = posx(where(posx!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1;
  posy = *y_wfs(1)._jstart + s2ipup ;
  posy = (transpose(posy * (*y_wfs(1)._isvalid)))(*);
  posy = posy(where(posy!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1;

  // Conversion en metres
  sspDiam = posx(2) - posx(1);
  p2m = (y_tel.diam/y_wfs(1).nxsub) / sspDiam;
  posx *= p2m;
  posy *= p2m;
  sspDiam = sspDiam*p2m;
  
  //Calcul de la position des points A,B (centre des bords gauche et droit des ssp)
  Ax = posx;
  Ay = posy + sspDiam/2.;
  Bx = posx + sspDiam;
  By = Ay;
  //Calcul de la position des points C,D (centre des bords bas et haut des ssp)
  Cx = posx + sspDiam/2.;
  Cy = posy;
  Dx = Cx;
  Dy = posy + sspDiam;

  //error;
  cmm = array(0.0f,2*y_wfs(1)._nvalid,2*y_wfs(1)._nvalid);
  seeing = y_wfs(1).lambda * 1e-6 / y_atmos.r0 * RASC;
  for (cs=1 ; cs<= y_atmos.nscreens ; cs++){
    delta_h = abs(y_dm(1).alt - (*y_atmos.alt)(cs));
    L0 = *y_atmos.L0(cs);
    x0 = 0.;
    k1 = RASC * y_wfs(1).lambda * 1.e-6 / 2. / pi / sspDiam;
    //error;
    for(i = 1 ; i<=dimsof(cmm)(2)/2 ; i++){
      //Covariance  xx
      for(j=1 ; j<=dimsof(cmm)(3)/2 ; j++){
        aa = sqrt((Ax(i) - Ax(j))^2 + (Ay(i) - Ay(j))^2);
        bb = sqrt((Bx(i) - Bx(j))^2 + (By(i) - By(j))^2);
        ab = sqrt((Ax(i) - Bx(j))^2 + (Ay(i) - By(j))^2);
        ba = sqrt((Bx(i) - Ax(j))^2 + (By(i) - Ay(j))^2);
      
        cmm(i,j) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi(ab,x0) + dphi(ba,x0) - dphi(aa,x0) - dphi(bb,x0));
      }
      j=1;
      //Covariance  xy
      for(jj=dimsof(cmm)(3)/2 + 1 ; jj<=dimsof(cmm)(3) ; jj++){
        ac = sqrt((Ax(i) - Cx(j))^2 + (Ay(i) - Cy(j))^2);
        ad = sqrt((Ax(i) - Dx(j))^2 + (Ay(i) - Dy(j))^2);
        bc = sqrt((Bx(i) - Cx(j))^2 + (By(i) - Cy(j))^2);
        bd = sqrt((Bx(i) - Dx(j))^2 + (By(i) - Dy(j))^2);
        j++;

        cmm(i,jj) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi(bc,x0) + dphi(ad,x0) - dphi(ac,x0) - dphi(bd,x0));
      }
    }
    //Covariance yy
    ii=1;
    jj=1;
    for(i=dimsof(cmm)(3)/2 + 1 ; i<=dimsof(cmm)(3) ; i++){
      for(j=dimsof(cmm)(3)/2 + 1 ; j<=dimsof(cmm)(3) ; j++){
        cc = sqrt((Cx(ii) - Cx(jj))^2 + (Cy(ii) - Cy(jj))^2);
        dd = sqrt((Dx(ii) - Dx(jj))^2 + (Dy(ii) - Dy(jj))^2);
        cd = sqrt((Cx(ii) - Dx(jj))^2 + (Cy(ii) - Dy(jj))^2);
        dc = sqrt((Dx(ii) - Cx(jj))^2 + (Dy(ii) - Cy(jj))^2);
        jj++;
      
        cmm(i,j) += 0.5 * (y_atmos.r0)^(-5./3.) * (*y_atmos.frac)(cs) * (dphi(cd,x0) + dphi(dc,x0) - dphi(cc,x0) - dphi(dd,x0));
      }
      ii++;
      jj=1;
    }
  }
  cmm(dimsof(cmm)(3)/2+1:,:dimsof(cmm)(3)/2) = cmm(:dimsof(cmm)(3)/2,dimsof(cmm)(3)/2+1:);
  
  return k1*k1*cmm;
}

func create_nact_geom(nm) {
  
  nb_act = y_dm(nm)._ntotact;
  nact = array(0.0f,nb_act,nb_act);
  coupling = y_dm(1).coupling;

  // Actuators positions
  tmpx = *y_dm(nm)._i1;
  tmpy = *y_dm(nm)._j1;
  offs = ((y_dm(nm)._n2-y_dm(nm)._n1+1) - (max(tmpx) - min(tmpx)))/2 - min(tmpx);
  tmpx += offs+1;
  tmpy += offs+1;
  mask = yoga_getdm(g_dm,y_dm(nm).type,y_dm(nm).alt) * 0;
  mask(tmpx + (tmpy-1) * dimsof(mask)(2)) = 1;
  masq_act = where(mask);
  pitch = masq_act(2) - masq_act(1);
  //error;
  for (i=1;i<=y_dm(nm)._ntotact;i++) {
    shape=mask*0;
    //Diagonal
    shape(tmpx(i) + (tmpy(i)-1) * dimsof(mask)(2)) = 1;
    // Left, right, above and under the current actuator
    shape(tmpx(i) - pitch + (tmpy(i)-1) * dimsof(mask)(2)) = coupling;
    shape(tmpx(i) + (tmpy(i) - pitch -1) * dimsof(mask)(2)) = coupling;
    shape(tmpx(i) + pitch + (tmpy(i)-1) * dimsof(mask)(2)) = coupling;
    shape(tmpx(i) + (tmpy(i) + pitch -1) * dimsof(mask)(2)) = coupling;
    // Diagonals of the current actuators 
    shape(tmpx(i) - pitch + (tmpy(i) - pitch -1) * dimsof(mask)(2)) = coupling^2 ;
    shape(tmpx(i) + pitch + (tmpy(i) - pitch -1) * dimsof(mask)(2)) = coupling^2 ;
    shape(tmpx(i) - pitch + (tmpy(i) + pitch -1) * dimsof(mask)(2)) = coupling^2 ;
    shape(tmpx(i) + pitch + (tmpy(i) + pitch -1) * dimsof(mask)(2)) = coupling^2 ;

    nact(i,) = shape(*)(masq_act);
  }

  return nact;
}

func openLoopSlp(nrec,nctrl){
  ol_slopes = array(0.0f,sum(y_wfs._nvalid * 2),nrec);
  for(sc=1 ; sc <= nrec ; sc++){
    write,format="\r Recording %d open-loop slopes : %d %%",y_controllers(nctrl).nrec,int(sc/float(nrec)*100);
    if (g_target != []) move_atmos,g_atmos;
    if ((y_target != []) && (g_target != [])) {
      // loop on targets
      for (tt=1;tt<=y_target.ntargets;tt++) {
        target_atmostrace,g_target,tt-1,g_atmos;
      }
    }
    
    if ((y_wfs != []) && (g_wfs != [])) {
      // loop on wfs
      for (ww=1;ww<=numberof(y_wfs);ww++) {
        sensors_trace,g_wfs,ww-1,"atmos",g_atmos;

        sensors_compimg,g_wfs,ww-1;
	
        sensors_compslopes,g_rtc,nctrl-1;
        ol_slopes((ww-1)*y_wfs(ww)._nvalid*2+1:ww*y_wfs(ww)._nvalid*2,sc) = sensors_getslopes(g_wfs,ww-1);
      }
    }
  }

  return ol_slopes;
}

func compute_KL2V(nctrl){
  KL2V = array(0.0f,sum(y_dm._ntotact),sum(y_dm._ntotact));
  ndms = *y_controllers(nctrl).ndm;
  nTT = 0;
  indx_act = 0;
  for (i=1 ; i <= numberof(ndms) ; i++){
    if(y_dm(ndms(i)).type == "pzt"){
      KL2V(indx_act + 1 : indx_act + y_dm(ndms(i))._ntotact,indx_act + 1 : indx_act + y_dm(ndms(i))._ntotact) = compute_klbasis(ndms(i));
      indx_act += y_dm(ndms(i))._ntotact
        }
    if(y_dm(ndms(i)).type == "tt")
      nTT += 1; 
    
  }
  if (nTT == 0)
    KL2V = KL2V(,:y_controllers(nctrl).nmodes);
  
  else{
    KL2V = KL2V(,:y_controllers(nctrl).nmodes);
    KL2V(,:y_controllers(nctrl).nmodes-2) = KL2V(,3:);
    KL2V(,y_controllers(nctrl).nmodes-1:) = array(0.0f,sum(y_dm._ntotact),2);
    KL2V(sum(y_dm._ntotact) - 1:,y_controllers(nctrl).nmodes-1:) = unit(2);
  }
  return KL2V;
}
