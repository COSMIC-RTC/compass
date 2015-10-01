require,"yoga_ao_utils.i"

  //----------------------------------------------------
func make_pzt_dm(nm,&influ,disp=)

  /* DOCUMENT function make_pzt_dm_elt(dm_structure,disp=)
     the influence functions are in microns per volt.
     same as make_pzt_dm but returns only local IF and
     start indices
   */
{
  extern y_dm;

  coupling=y_dm(nm).coupling;

  // best parameters, as determined by a multi-dimensional fit
  // (see coupling3.i)
  a=[4.49469,7.25509,-32.1948,17.9493];
  p1 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [2.49456,-0.65952,8.78886,-6.23701];
  p2 = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  a = [1.16136,2.97422,-13.2381,20.4395];
  irc = a(1)+a(2)*coupling+a(3)*coupling^2+a(4)*coupling^3;

  dim   = y_dm(nm)._n2-y_dm(nm)._n1+1;
  size  = y_geom.ssize;
  nxact = y_dm(nm).nact;
  cent  = y_geom.cent;
  pitch = y_dm(nm)._pitch;
  ir    = irc*pitch;

  tmp=pitch/abs(ir);
  c = (coupling - 1.+ tmp^p1)/(log(tmp)*tmp^p2);

  // compute IF on partial (local) support:
  smallsize = long(ceil(2*ir+10));
  smallsize += smallsize&1;
  y_dm(nm)._influsize = smallsize;
  x     = indgen(smallsize)(,-:1:smallsize)-smallsize/2 - 0.5;
  y     = transpose(x);
  tmpx  = clip(abs(x/ir),1e-8,2.);
  tmpy  = clip(abs(y/ir),1e-8,2.);
  tmp   = (1.-tmpx^p1+c*log(tmpx)*tmpx^p2)*(1.-tmpy^p1+c*log(tmpy)*tmpy^p2);
  influ   = tmp*(tmpx <= 1.)*(tmpy <= 1.);

  // compute location (x,y and i,j) of each actuator:
  cub   = array(float,nxact,nxact,2);
  // make X and Y indices array:
  xy    = indgen(nxact)(,-:1:nxact);
  xy    = [xy,transpose(xy)];
  // express "centered" coordinate of actuator in pixels:
  xy    = (xy-1.-(nxact-1.)/2.)*pitch;
  
  // fill cub (X coord  and Y coord):
  cub(,,1) = xy(,,1); cub(,,2) = xy(,,2);

  // the following determine if an actuator is to be considered or not
  // relative to the pitchmargin parameter.
  dis      = sqrt(cub(,,1)^2.+cub(,,2)^2.);
  if (y_dm(nm).margin == 0) {
    pitchMargin = 1.44;
  } else {
    pitchMargin = y_dm(nm).margin;
  }
  rad      = ((nxact-1.)/2.+pitchMargin)*pitch;
  inbigcirc= where(dis < rad);
  // 1 if valid actuator, 0 if not:

  // converting to array coordinates:
  cub += cent;

  cub      = cub(*,);
  // cub now has two indices: first one is actuator number
  // second one is: 1:Xcoord, 2:Ycoord

  // filtering actuators outside of a disk radius = rad (see above)
  cubval   = cub(inbigcirc,);
  //cubval = cub;

  y_dm(nm)._ntotact = dimsof(cubval)(2);
  y_dm(nm)._xpos    = &(cubval(,1));
  y_dm(nm)._ypos    = &(cubval(,2));
  y_dm(nm)._i1      = &(int(long(cubval(,1)-smallsize/2+0.5)-y_dm(nm)._n1));
  y_dm(nm)._j1      = &(int(long(cubval(,2)-smallsize/2+0.5)-y_dm(nm)._n1));

  influ = influ(,,-)*array(1.f,y_dm(nm)._ntotact)(-,-,);

  if (y_dm(nm)._puppixoffset!=[]) {
    // see comment above in make_pzt_dm
    *y_dm(nm)._xpos += y_dm(nm)._puppixoffset(1);
    *y_dm(nm)._ypos += y_dm(nm)._puppixoffset(2);
    //y_dm(nm)._i1     = &(int(long(*y_dm(nm)._i1+y_dm(nm)._puppixoffset(1))));
    //y_dm(nm)._j1     = &(int(long(*y_dm(nm)._j1+y_dm(nm)._puppixoffset(2))));
  }

  fact = y_dm(nm).unitpervolt/max(influ);
  influ = float(influ*fact);
  y_dm(nm)._influ = &influ;
  comp_dmgeom,nm;
  
  // Prepare kernel convolution for comp dm_shape
  dims       = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
  dim       = dimsof(*y_geom._mpupil)(2);
  if (dims >= dim) dim = dims;
  kimg = influ(,,1);
  kernconv = array(0.,dim,dim);
  off = (dim - y_dm(nm)._influsize)/2;
  kernconv(off+1:off+y_dm(nm)._influsize,off+1:off+y_dm(nm)._influsize) = kimg;
  y_dm(nm)._influkernel = &(roll(kernconv));

  return influ;
}

func comp_dmgeom(nm)
{
  extern y_dm;

  smallsize = y_dm(nm)._influsize;
  nact      = y_dm(nm)._ntotact;
  dims      = long(y_dm(nm)._n2-y_dm(nm)._n1+1);
  dim       = dimsof(*y_geom._mpupil)(2);
  if (dims < dim) offs = (dim - dims)/2;
  else {
    offs = 0;
    dim  = dims;
  }
  mapactu = array(1,dim,dim);
  if (offs > 0) {
    mapactu(1:offs,) = 0;
    mapactu(,1:offs) = 0;
    mapactu(,-offs+1:0) = 0;
    mapactu(-offs+1:0,) = 0;
  }
  
  tmpx = indgen(smallsize)(,-:1:smallsize)-1;
  tmpy = transpose(tmpx);
  tmpx = offs + tmpx(,,-:1:nact) + (*y_dm(nm)._i1)(-:1:smallsize,-:1:smallsize,);
  tmpy = offs + tmpy(,,-:1:nact) + (*y_dm(nm)._j1)(-:1:smallsize,-:1:smallsize,);
  tmp = tmpx + dim * tmpy;
  //tmp = tmpy + dim * tmpx;
  
  if (numberof(where(tmpx < 0)) > 0) tmp(where(tmpx < 0)) = -10;
  if (numberof(where(tmpy < 0)) > 0) tmp(where(tmpy < 0)) = -10;
  if (numberof(where(tmpx > dims-1)) > 0) tmp(where(tmpx > dims-1)) = -10;
  if (numberof(where(tmpy > dims-1)) > 0) tmp(where(tmpy > dims-1)) = -10;
  
  itmps = sort(tmp(*));
  tmps = tmp(itmps);
  istart = where(tmps>-1)(1);
  tmps = tmps(istart:);
  itmps = itmps(istart:);

  npts = istart = array(0,dim*dim);
  cpt=1;
  ref=cpt;

  for (i=1;i<=dim*dim;i++) {
    if ((offs != 0) && (mapactu(*)(i) == 0)) npts(i)=0;
    else {
      while ((tmps(cpt) == i-1) && (cpt<numberof(tmps))) cpt++;
      npts(i)   = cpt-ref;
      istart(i) = ref;
      ref       = cpt;
    }
  }

  itmps = itmps(1:sum(npts));
  y_dm(nm)._influpos   = &long(itmps-1);
  y_dm(nm)._ninflu     = &int(npts);
  y_dm(nm)._influstart = &int(istart-1);
  
  tmp = *y_dm(nm)._i1;
  tmp += offs;
  y_dm(nm)._i1 = &int(tmp);
  tmp = *y_dm(nm)._j1;
  tmp += offs;
  y_dm(nm)._j1 = &int(tmp);
  
}

func build_dm(nm,command)
// algorithm to be parallelized
{
  dim   = y_dm(nm)._n2-y_dm(nm)._n1+1;
  dim2 = dim * dim;
  diminf = y_dm(nm)._influsize * y_dm(nm)._influsize;
  dmshape = array(0.0f, dim,dim);
  for (i=1;i<=dim*dim;i++) {
    for (j=0;j<(*y_dm(nm)._ninflu)(i);j++) {
      pos = (*y_dm(nm)._influpos)((*y_dm(nm)._influstart)(i)+1+j);
      ninflu = pos / diminf + 1;
      dmshape(i) += command(ninflu) * ((*y_dm(nm)._influ)(*)(pos));
    }
  }
  return dmshape;
}

func build_dm_gpu(nm,com)
// algorithm to be parallelized
{
  yoga_setcomm,g_dm,y_dm(nm).type,y_dm(nm).alt,com;
  yoga_shapedm,g_dm,y_dm(nm).type,y_dm(nm).alt;
  dm_shape = yoga_getdm(g_dm,y_dm(nm).type,y_dm(nm).alt);
  return dm_shape;
}

//----------------------------------------------------
func make_tiptilt_dm(nm,&def,disp=)

  /* DOCUMENT function make_tiptilt_dm,ndm,&def,disp=
     adapted from makeZernikeIF
     modified 2004jan22 to make it normalized at 1"
   */
{
  extern y_dm;
  
  dims      = y_dm(nm)._n2-y_dm(nm)._n1+1;
  dim       = dimsof(*y_geom._mpupil)(2);
  if (dims >= dim) dim = dims;
  
  nzer  = 2;
  cobs  = y_tel.cobs;
  cent  = y_geom.cent;
  psize = y_tel.diam/y_geom.pupdiam;
  patchDiam = long(y_geom.pupdiam+2*max(abs(y_wfs.xpos,y_wfs.ypos))*
                   4.848e-6*abs(y_dm(nm).alt)/psize);

  influ = make_zernike(nzer+1,dim,patchDiam,
                       y_geom.cent-y_dm(nm)._n1+1,y_geom.cent-y_dm(nm)._n1+1,1)(,,2:);

  // normalization factor: one unit of tilt gives 1 arcsec:
  current = influ(dim/2,dim/2,1)-influ(dim/2-1,dim/2,1);
  fact = (y_dm(nm).unitpervolt*y_tel.diam/y_geom.pupdiam)*4.848/current;

  influ = float(influ*fact);
  y_dm(nm)._ntotact = (dimsof(influ))(4);
  y_dm(nm)._influ = &influ;

  return influ;
}

//----------------------------------------------------
func make_kl_dm(nm,disp=)

  /* DOCUMENT function make_kl_dm,ndm,disp=
   */
{
  extern y_dm;
  
  //dim   = y_dm(nm)._n2-y_dm(nm)._n1+1;
  dim       = dimsof(*y_geom._mpupil)(2);
  cobs  = y_tel.cobs;
  cent  = y_geom.cent;
  psize = y_geom.pupdiam;
  patchDiam = long(y_geom.pupdiam+2*max(abs(y_wfs.xpos,y_wfs.ypos))*
                   4.848e-6*abs(y_dm(nm).alt)/psize);

  klbas = make_klbas(y_dm(nm).nkl,cobs,patchDiam,funct="kolmo");

  // normalization factor: one unit of tilt gives 1 arcsec:
  //current = influ(dim/2,dim/2,1)-influ(dim/2-1,dim/2,1);
  //fact = (y_dm(nm).unitpervolt*y_tel.diam/y_geom.pupdiam)*4.848/current;

  //influ = float(influ*fact);
  //y_dm(nm)._ntotact = (dimsof(influ))(4);

  y_dm(nm)._klbas = &klbas;
  y_dm(nm)._i1      = &(array(int((dim-patchDiam)/2.0f),y_dm(nm).nkl));
  y_dm(nm)._j1      = &(array(int((dim-patchDiam)/2.0f),y_dm(nm).nkl));
  y_dm(nm)._ntotact = y_dm(nm).nkl;
}



// Florian features -------------------------------------
func make_flo_kl_dm(nm,disp=)

  /* DOCUMENT function make_kl_dm,ndm,disp=
   */
{
  extern y_dm;
  
  dim   = y_dm(nm)._n2-y_dm(nm)._n1+1;
  dim       = dimsof(*y_geom._mpupil)(2);
  cobs  = y_tel.cobs;
  cent  = y_geom.cent;
  psize = y_geom.pupdiam;
  patchDiam = long(y_geom.pupdiam+2*max(abs(y_wfs.xpos,y_wfs.ypos))*
                   4.848e-6*abs(y_dm(nm).alt)/psize);

  klbas = make_flo_klbas(y_dm(nm).nkl,cobs);

  y_dm(nm)._klbas = &klbas;
  y_dm(nm)._i1      = &(array(int((dim-patchDiam)/2.0f),y_dm(nm).nkl));
  y_dm(nm)._j1      = &(array(int((dim-patchDiam)/2.0f),y_dm(nm).nkl));
  y_dm(nm)._ntotact = y_dm(nm).nkl;
}

func compute_klbasis(ndm){
  /* DOCUMENT function compute_klbasis,ndm

     return the KL basis on pzt actuators (only for pzt dm)
  */
  if(y_dm(ndm).type == "pzt"){
    tmp = (dimsof(*y_geom._ipupil)(2)-(y_dm(ndm)._n2 - y_dm(ndm)._n1 +1))/2;
    pup = (*y_geom._ipupil)(tmp+1:-tmp,tmp+1:-tmp);
    indx_valid = where(pup > 0)-1; 
    x = *y_dm(ndm)._xpos;
    y = *y_dm(ndm)._ypos; 
    interactp = x(2) - x(1);
    interactm = y_tel.diam/(y_dm(ndm).nact-1);
    p2m = interactm/interactp;
    norm = -(p2m*y_tel.diam/(2*y_atmos.r0))^(5./3);
    yoga_computeKLbasis,g_dm,"pzt",y_dm(ndm).alt,x,y,indx_valid,numberof(indx_valid),norm,1.0f;
    KLbasis = yoga_getKLbasis(g_dm,"pzt",y_dm(ndm).alt)(,::-1);
  }
  else{
    KLbasis = [];
    write,"DM must be pzt type";
  }
  
  return KLbasis;

}

func setcomkl(ndm,comvec,plot=){
  /*DOCUMENT function setcomkl,ndm,comvec,plot=

    Set the KL command vector comvec on the pzt DM #ndm and return the pzt command vector
applied
    Set plot=1 to apply the command and plot the DM shape

    SEE ALSO : compute_klbasis
  */
  if(y_dm(ndm).type == "pzt"){
    yoga_setcomkl,g_dm,"pzt",y_dm(ndm).alt,comvec(::-1); // Reverse vector because KL basis is reversed on GPU
    if(plot){
      yoga_shapedm,g_dm,"pzt",y_dm(ndm).alt;
      dm = yoga_getdm(g_dm,"pzt",y_dm(ndm).alt);
      window,0;
      fma;
      pli,dm;
    }
  
    return yoga_getcomm(g_dm,"pzt",y_dm(ndm).alt);
  }
  else{
    write,"DM must be pzt type";
    return 0;
  }
}

