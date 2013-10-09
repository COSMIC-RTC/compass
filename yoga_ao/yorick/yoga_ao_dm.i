require,"yoga_ao_utils.i"

  //----------------------------------------------------
func make_pzt_dm(nm,&influ,disp=)

  /* DOCUMENT function make_pzt_dm_elt(dm_structure,disp=)
     the influence functions are in microns per volt.
     same as make_pzt_dm but returns only local IF and
     start indices
   */
{
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
  pitch = y_dm(nm).pitch;
  ir = irc*pitch;

  tmp=pitch/abs(ir);
  c = (coupling - 1.+ tmp^p1)/(log(tmp)*tmp^p2);

  // compute IF on partial (local) support:
  smallsize = long(ceil(2*ir+10));
  y_dm(nm)._influsize = smallsize;
  x     = indgen(smallsize)(,:1:smallsize)-smallsize/2-0.5;
  y     = transpose(x);
  tmpx  = clip(abs(x/ir),1e-8,2.);
  tmpy  = clip(abs(y/ir),1e-8,2.);
  tmp   = (1.-tmpx^p1+c*log(tmpx)*tmpx^p2)*(1.-tmpy^p1+c*log(tmpy)*tmpy^p2);
  influ   = tmp*(tmpx <= 1.)*(tmpy <= 1.);

  // compute location (x,y and i,j) of each actuator:
  cub   = array(float,nxact,nxact,2);
  // make X and Y indices array:
  xy    = indgen(nxact)(,:1:nxact);
  xy    = [xy,transpose(xy)];
  // express "centered" coordinate of actuator in pixels:
  xy    = (xy-1.-(nxact-1.)/2.)*pitch;

  // fill cub (X coord  and Y coord):
  cub(,,1) = xy(,,1); cub(,,2) = xy(,,2);

  // the following determine if an actuator is to be considered or not
  // relative to the pitchmargin parameter.
  dis      = sqrt(cub(,,1)^2.+cub(,,2)^2.);
  if (y_dm(nm).pitchMargin == 0) {
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

  y_dm(nm)._ntotact = dimsof(cubval)(2);
  y_dm(nm)._xpos    = &(cubval(,1));
  y_dm(nm)._ypos    = &(cubval(,2));
  y_dm(nm)._i1      = &(int(long(cubval(,1)-smallsize/2+0.5)-y_dm(nm)._n1));
  y_dm(nm)._j1      = &(int(long(cubval(,2)-smallsize/2+0.5)-y_dm(nm)._n1));

  influ = influ(,,-)*array(1.f,y_dm(nm)._nact)(-,-,);

  if (y_dm(nm)._puppixoffset!=[]) {
    // see comment above in make_pzt_dm
    *y_dm(nm)._xpos += y_dm(nm)._puppixoffset(1)
    *y_dm(nm)._ypos += y_dm(nm)._puppixoffset(2)
  }

  fact = y_dm(nm).unitpervolt/max(influ);
  influ = float(influ*fact);
  y_dm(nm)._influ = &influ;

  return influ;
}


//----------------------------------------------------
func make_tiptilt_dm(nm,&def,disp=)

  /* DOCUMENT function make_tiptilt_dm,dm_structure,ActIF,disp=
     adapted from makeZernikeIF
     modified 2004jan22 to make it normalized at 1"
   */
{
  dim   = y_dm(nm)._n2-y_dm(nm)._n1+1;
  nzer  = 2;
  cobs  = y_tel.cobs;
  cent  = y_geom.cent;
  psize = y_tel.diam/y_geom.pupdiam;
  patchDiam = y_geom.pupdiam+2*max(abs([y_wfs.xpos,y_wfs.ypos]))*
    4.848e-6*abs(y_dm(nm).alt)/psize;

  prepzernike,;

  influ = make_zernike(nzer+1,dim,patchDiam,y_geom.
                       cent-y_dm(nm)._n1+1,y_geom.cent-y_dm(nm)._n1+1,1)(,,2:);

  // normalization factor: one unit of tilt gives 1 arcsec:
  current = influ(dim/2,dim/2,1)-influ(dim/2-1,dim/2,1);
  fact = (y_dm(nm).unitpervolt*y_tel.diam/y_geom.pupdiam)*4.848/current;

  influ = float(influ*fact);
  y_dm(nm)._nact = (dimsof(influ))(4);
  y_dm(nm)._influ = &influ;

  return influ;
}

