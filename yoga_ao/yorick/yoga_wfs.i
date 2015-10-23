require,yoga_ao_top+"/yorick/yoga_aolib.i";
require,yoga_ao_top+"/yorick/yoga_ao_ystruct.i";
require,yoga_ao_top+"/yorick/yoga_ao_utils.i"

func init_wfs_size(wfs,&pdiam,&Nfft,&Ntot,&nrebin,&pixsize,&qpixsize)
/* DOCUMENT init_wfs_size
 init_wfs_size,wfs,pdiam,Nfft,Ntot,nrebin,pixsize,qpixsize

 computes the quatum pixel sizes and all useful dimensions for a given wfs
 requires 2 externals : y_atmos & y_tel
 y_atmos : a atmos_struct()
 y_tel   : a tel_struct()
 wfs     : wfs_struct as input (input)
 pdiam   : pupil diam for each subap (pixels) (output)
 Nfft    : fft size for a subap (pixels) (output)
 Ntot    : hr image size for a subap (pixels) (output)
 nrebin  : rebin factor for a subap (output)
 pixsize : simulated pixel size for a subap (arcsec) (output)
 qpixsize: quantum pixel size for a subap (arcsec) (output)
 
 SEE ALSO:
 */
{
  /*
   Scheme to determine arrays sizes
   sh :
   k = 6
   p = k * d/r0
   n = long(2*d*v/lambda/RASC)+1
   N = fft_goodsize(k*n/v*lambda/r0*RASC)
   u = k * lambda / r0 * RASC / N
   n = v/u - long(v/u) > 0.5 ? long(v/u)+1 : long(v/u)
   v = n * u
   Nt = v * Npix

   pyr :
   
   */

  extern y_atmos,y_tel;

  r0 = y_atmos.r0 * ((wfs.lambda/0.5)^(6./5));
  if (r0 != 0) {
    write,format="r0 for WFS : %.2f m\n",r0;
    seeing = RASC * (wfs.lambda * 1.e-6) / r0;
    write,format="seeing for WFS : %.2f \"\n",seeing;
  }

  if (pdiam == []) {
    // this case is usualy for the wfs with max # of subaps
    // we look for the best compromise between pixsize and fov

    subapdiam = y_tel.diam / double(wfs.nxsub);// diam of subap

    k = 6;

    pdiam = long(k * subapdiam / r0);// number of phase points per subap
    if (pdiam < 16) pdiam = 16;
    //pdiam = 18;

    if (wfs.type == "sh") {
      nrebin = long(2 * subapdiam * wfs.pixsize / (wfs.lambda*1.e-6) / RASC) + 1;
      nrebin = max(2,nrebin);
      // first atempt on a rebin factor

      //Nfft = fft_goodsize(k * nrebin / wfs.pixsize * (wfs.lambda*1.e-6) / r0 * RASC);
      // since we clipped pdiam we have to be carreful in nfft computation
      Nfft = fft_goodsize(pdiam / subapdiam * nrebin / wfs.pixsize * (wfs.lambda*1.e-6) * RASC);
      // size of the support in fourier domain

      //qpixsize = k * (wfs.lambda*1.e-6) / r0 * RASC / Nfft;
      qpixsize = (pdiam * (wfs.lambda*1.e-6) / subapdiam * RASC) / Nfft;
    }

    if ((wfs.type == "pyr")||(wfs.type == "roof")) {
      //while (pdiam % wfs.npix != 0) pdiam+=1;
      padding = 2;
      nphase = pdiam * wfs.nxsub+ 2 * padding * pdiam;
      qpixsize = (pdiam * wfs.nxsub / double(nphase))*wfs.lambda/y_tel.diam/4.848;

      fssize_pixels = long(wfs.fssize / qpixsize / 2.);
      //nrebin  = pdiam / wfs.npix;

      npix = (wfs.nxsub+2*padding);
      npix_ok = npix*indgen(100);
      npix_ok = npix_ok(where(npix_ok<=(nphase/2.)));

      // now we want the quadrant image dim to be > fssize_pixels:
      w = where(npix_ok>=fssize_pixels);
      if (numberof(w)==0) {
        maxfs = npix_ok(0)*2*psize;
        error,swrite(format="wfs(%i).fssize too large (max=%.3f\")!",n,maxfs);
      }

      npix_ok = npix_ok(w(1));

      nrebin = npix_ok/npix;
      Nfft = npix_ok;
      Ntot = nphase;
      pixsize = qpixsize * nrebin;
    }

    // quantum pixel size
  } else {
    // this case is for a wfs with fixed # of phase points
    subapdiam = y_tel.diam / double(wfs.nxsub);// diam of subap

    if (wfs.type != "geo") {
      Nfft = fft_goodsize(2*pdiam);
      // size of the support in fourier domain

      qpixsize = pdiam * (wfs.lambda*1.e-6) / subapdiam * RASC / Nfft;
      // quantum pixel size
    }
  }

  if (wfs.type == "sh") {
    // actual rebin factor
    nrebin = ((wfs.pixsize/qpixsize - long(wfs.pixsize/qpixsize) > 0.5)(1) ?
        long(wfs.pixsize/qpixsize)+1 : long(wfs.pixsize/qpixsize));

    // actual pixel size
    pixsize = nrebin * qpixsize;

    if ((pixsize * wfs.npix > qpixsize * Nfft)(1))
    Ntot = fft_goodsize(long(pixsize * wfs.npix / qpixsize) + 1);
    else Ntot = Nfft;

    if (Ntot%2 != Nfft%2) Ntot+=1;
  }

  if (wfs.type != "geo") {
    write,format="quantum pixsize : %.4f \"\n",qpixsize;
    write,format="simulated FoV : %.2f\" x %.2f\"\n",Ntot * qpixsize,Ntot * qpixsize;
    write,format="actual pixsize : %.2f\"\n",pixsize;
    write,format="actual FoV : %.2f\" x %.2f\"\n",pixsize * wfs.npix,
    pixsize * wfs.npix;
    write,format="number of phase points : %d\n",pdiam;
    write,format="size of fft support : %d\n",Nfft;
    if (Ntot > Nfft) write,format="size of HR spot support : %d\n",Ntot;
  }
}

func init_wfs_geom(n,init=)
/* DOCUMENT init_wfs_geom
 init_wfs_geom,n,init=

 inits a wfs geometry (arrays used for image manipulation in the wfs model)
 requires 2 externals : y_atmos & y_wfs
 y_atmos : a atmos_struct()
 y_wfs   : a wfs_struct()
 n       : index of wfs to init
 init    : if init = 1 inits the whole simulation geometry
 
 SEE ALSO:
 */
{
  extern y_wfs;

  write,format="\n*-----------------------\nDoing inits on WFS # %d\n",n;

  if ((init == []) || (init == 0)) {
    if (y_wfs(n).type == "sh")
    pdiam = (y_geom.pupdiam % y_wfs(n).nxsub == 0 ? long(y_geom.pupdiam / y_wfs(n).nxsub) :
        long(y_geom.pupdiam / y_wfs(n).nxsub) + 1);
    if (y_wfs(n).type == "geo") {
      if (y_geom.pupdiam == 0)
      pdiam = 20;
      else
      pdiam = (y_geom.pupdiam % y_wfs(n).nxsub == 0 ? long(y_geom.pupdiam / y_wfs(n).nxsub) :
          long(y_geom.pupdiam / y_wfs(n).nxsub) + 1);
    }
  } else pdiam = [];
  
  init_wfs_size,y_wfs(n),pdiam,Nfft,Ntot,nrebin,pixsize,qpixsize;

  if (y_wfs(n).type != "geo") {
    y_wfs(n).pixsize = pixsize;
    y_wfs(n)._Nfft = Nfft;
    y_wfs(n)._Ntot = Ntot;
    y_wfs(n)._nrebin = nrebin;
    y_wfs(n)._qpixsize = qpixsize;
  }
  y_wfs(n)._subapd = y_tel.diam/y_wfs(n).nxsub;
  y_wfs(n)._pdiam = pdiam;

  if ((y_wfs(n).type == "pyr")||(y_wfs(n).type == "roof")) y_wfs(n).npix = pdiam;

  if ((init == 1) || ((y_wfs(n).type == "geo") && (n==1))) {
    //this is the wfs with largest # of subaps
    //the overall geometry is deduced from it
    geom_init,pdiam * y_wfs(n).nxsub;
  }

  if ((y_wfs(n).type == "pyr")||(y_wfs(n).type == "roof")) {
    padding = 2;
    npup = y_wfs(n)._Ntot;
    n1 = y_geom.ssize / 2 - y_geom.pupdiam / 2 - padding * y_wfs(n).npix+1;
    n2 = n1 + y_geom.pupdiam + 2 * padding * y_wfs(n).npix-1;

    y_geom._mpupil = &((*y_geom._ipupil)(n1:n2,n1:n2));
    y_geom._n1 = n1;
    y_geom._n2 = n2;
    y_geom._n = npup;

    //pup   = pup(ii,ii);
    //phase = phase(ii,ii);
    mod_ampl_pixels = y_wfs(n).pyr_ampl / y_wfs(n)._qpixsize;// modulation in pixels
    fsradius_pixels = long(y_wfs(n).fssize / y_wfs(n)._qpixsize / 2.);

    if (y_wfs(n).fstop=="round") {
      focmask = dist(npup,xc=npup/2.+0.5,yc=npup/2.+0.5)<(fsradius_pixels);
      fstop_area = pi * (y_wfs(n).fssize/2.)^2.;
    } else {
      if (y_wfs(n).fstop=="square") {
        xy = indices(npup)-(npup+1.)/2.;
        focmask = ( abs(xy(,,1)) <= (fsradius_pixels) ) *
        ( abs(xy(,,2)) <= (fsradius_pixels) );
        fstop_area = y_wfs(n).fssize^2.;
      } else error,swrite(format="wfs(%i).fstop must be round or square",n);
    }

    pyr_focmask = roll(focmask);
    /*
     if (y_wfs(n).pyr_loc!="after") {
     tmp = array(0,[3,Nfft,Nfft,4]);
     tmp(,,1) = pyr_focmask(npup-Nfft+1:,npup-Nfft+1:);
     tmp(,,2) = tmp(,,1)(::-1,);
     tmp(,,3) = tmp(,,1)(,::-1);
     tmp(,,4) = tmp(,,1)(::-1,::-1);
     y_wfs(n)._submask = &tmp;
     } else */y_wfs(n)._submask = &(pyr_focmask);

    pup = *y_geom._spupil;
    pupreb = bin2d(pup*1.,y_wfs(n).npix)/y_wfs(n).npix^2.;
    wsubok = where(pupreb>=y_wfs(n).fracsub);
    pupvalid = pupreb * 0.;
    pupvalid(wsubok) = 1;
    y_wfs(n)._isvalid = &int(pupvalid);

    pup = *y_geom._mpupil;
    pupreb = bin2d(pup*1.,y_wfs(n).npix)/y_wfs(n).npix^2.;
    wsubok = where(pupreb>=y_wfs(n).fracsub);
    pupvalid = pupreb * 0.;
    pupvalid(wsubok) = 1;
    y_wfs(n)._nvalid = numberof(wsubok);
    y_wfs(n)._validsubs = &int(where2(pupvalid));

    istart = jstart = long(span(0.5,y_geom.pupdiam + 2 * padding * y_wfs(n).npix + 0.5,y_wfs(n).nxsub+2*padding)+1)(:-1);
    y_wfs(n)._istart = &istart;
    y_wfs(n)._jstart = &jstart;

    xy = indices(npup)-(npup+1)/2.;
    phase_shift = roll( exp(1i*2*pi*(0.5*xy(,,sum))/npup) );

    y_wfs(n)._halfxy = &(phase_shift);

    xy = indices(Nfft)-(Nfft+1)/2.;
    coef1 = ( odd(y_wfs(n).nxsub*nrebin) ? 0.0:-0.5 );
    coef2 = ((odd(y_wfs(n).nxsub)&&odd(y_wfs(n).npix*nrebin)) ? 1.0:0.5 );
    pshift = exp(1i*2*pi*(coef1/Nfft+
            coef2*nrebin/y_wfs(n).npix/Nfft)*xy(,,sum));
    //roll,pshift;

    y_wfs(n)._pyr_offsets = &(pshift);

    if ((y_wfs(n).pyrtype) == "Pyramid") {

      if (*y_wfs(n).pyr_pos == []) {
        cx = lround(mod_ampl_pixels*sin(indgen(y_wfs(n).pyr_npts)*2.*pi/y_wfs(n).pyr_npts));
        cy = lround(mod_ampl_pixels*cos(indgen(y_wfs(n).pyr_npts)*2.*pi/y_wfs(n).pyr_npts));
        mod_npts = y_wfs(n).pyr_npts;
      } else {
        write,format="%s\n", "Using user-defined positions for the pyramid modulation";
        // user defined positions
        cx = lround((*y_wfs(n).pyr_pos)(:,1)/qpixsize);
        cy = lround((*y_wfs(n).pyr_pos)(:,2)/qpixsize);
        mod_npts = dimsof(cx)(2);
      }

    } else if ((y_wfs(n).pyrtype) == "RoofPrism") {
      cx = lround(2.*mod_ampl_pixels*(indgen(y_wfs(n).pyr_npts)-(y_wfs(n).pyr_npts+1)/2.)/y_wfs(n).pyr_npts);
      cy = cx;
      mod_npts = y_wfs(n).pyr_npts;

    } else {
      if (*y_wfs(n).pyr_pos == []) {
        cx = lround(mod_ampl_pixels*sin(indgen(y_wfs(n).pyr_npts)*2.*pi/y_wfs(n).pyr_npts));
        cy = lround(mod_ampl_pixels*cos(indgen(y_wfs(n).pyr_npts)*2.*pi/y_wfs(n).pyr_npts));
        mod_npts = y_wfs(n).pyr_npts;
      } else {
        write,format="%s\n", "Using user-defined positions for the pyramid modulation";
        // user defined positions
        cx = lround((*y_wfs(n).pyr_pos)(:,1)/qpixsize);
        cy = lround((*y_wfs(n).pyr_pos)(:,2)/qpixsize);
        mod_npts = dimsof(cx)(2);
      }
    }

    y_wfs(n)._pyr_cx = &(cx);
    y_wfs(n)._pyr_cy = &(cy);

    y_wfs(n)._nphotons = y_wfs(n).zerop*2.51189^(-y_wfs(n).gsmag)*
    y_loop.ittime*y_wfs(n).optthroughput;

    // spatial filtering by the pixel extent:
    // *2/2 intended. min should be 0.40 = sinc(0.5)^2.
    xy2 = xy/(Nfft-1)*2/2;
    sincar = roll(__mysinc(pi*xy2(,,1))*__mysinc(pi*xy2(,,2)));
    y_wfs(n)._hrmap = &float(sincar);

    // this defines how we cut the phase into subaps
    phasemap = array(0,pdiam,pdiam,y_wfs(n)._nvalid);
    tmp = indices(y_geom._n)-1;// we need c-like indices
    tmp = tmp(,,1)+tmp(,,2)*(y_geom._n);

    for (i=1;i<=y_wfs(n)._nvalid;i++) {
      indi = istart((*y_wfs(n)._validsubs)(1,i))+2;
      indj = jstart((*y_wfs(n)._validsubs)(2,i))+2;
      phasemap(,,i) = tmp(indi:indi+pdiam-1,indj:indj+pdiam-1);
    }
    y_wfs(n)._phasemap = &int(phasemap(*,));
  }

  if ((y_wfs(n).type == "sh") || (y_wfs(n).type == "geo")) {
    // this is the i,j index of lower left pixel of subap
    istart = jstart = long(span(0.5,y_geom.pupdiam + 0.5,y_wfs(n).nxsub+1)+1)(:-1);
    y_wfs(n)._istart = &istart;
    y_wfs(n)._jstart = &jstart;

    if(y_wfs(n).atmos_seen == 0){
      cent = y_geom.pupdiam/2 + 0.5;
      pup = float(make_pupil(y_geom.pupdiam,y_geom.pupdiam,type_ap=y_tel.type_ap,angle=y_tel.pupangle,spiders_type=y_tel.spiders_type,t_spiders=y_tel.t_spiders,nbr_miss_seg=y_tel.nbrmissing,std_ref_err=y_tel.referr,xc=cent,yc=cent,real=,cobs=));
      pup = pad_array(pup,y_geom._n);
    }
    else pup = *y_geom._mpupil;
    //pup = *y_geom._mpupil;
    // sorting out valid subaps
    fluxPerSub = array(float,y_wfs(n).nxsub,y_wfs(n).nxsub);

    for (i=1;i<=y_wfs(n).nxsub;i++) {
      for (j=1;j<=y_wfs(n).nxsub;j++) {
        indi = istart(i)+2;
        indj = jstart(j)+2;
        fluxPerSub(i,j) = sum(pup(indi:indi+pdiam-1,indj:indj+pdiam-1));
      }
    }
    fluxPerSub = fluxPerSub/pdiam^2.;
    tmp = fluxPerSub >= y_wfs(n).fracsub;
    //tmp(where(tmp == 0)) = nvalid+10;
    y_wfs(n)._validsubs = &int(where2(tmp));
    y_wfs(n)._isvalid = &int(tmp);
    y_wfs(n)._nvalid = sum(*y_wfs(n)._isvalid);
    y_wfs(n)._fluxPerSub = &float(fluxPerSub);
    // this defines how we cut the phase into subaps
    phasemap = array(0,pdiam,pdiam,y_wfs(n)._nvalid);
    tmp = indices(y_geom._n)-1;// we need c-like indices
    tmp = tmp(,,1)+tmp(,,2)*(y_geom._n);
    for (i=1;i<=y_wfs(n)._nvalid;i++) {
      indi = istart((*y_wfs(n)._validsubs)(1,i))+2;
      indj = jstart((*y_wfs(n)._validsubs)(2,i))+2;
      phasemap(,,i) = tmp(indi:indi+pdiam-1,indj:indj+pdiam-1);
    }
    //verif
    //pli,reform((*y_geom._mpupil)(*)(phasemap(*,1)+1),pdiam,pdiam)
    y_wfs(n)._phasemap = &int(phasemap(*,));

    //this is a phase shift of 1/2 pix in x and y
    if (y_wfs(n).type == "sh") {
      halfxy = (span(0,2*pi,y_wfs(n)._Nfft+1)(1:y_wfs(n)._pdiam) / 2.)(,-:1:y_wfs(n)._pdiam);
      halfxy += transpose(halfxy);

      y_wfs(n)._halfxy = &float(halfxy*0.);

      //  if (y_wfs(n).gsalt == 0.) {
      /*
       if ((y_wfs(n).npix % 2 < y_wfs(n)._Nfft % 2) ||
       (y_wfs(n).npix % 2 != y_wfs(n)._pdiam % 2)) 
       y_wfs(n)._halfxy = &float(halfxy);
       else y_wfs(n)._halfxy = &float(halfxy*0.);
       
       if (y_wfs(n).gsalt != 0.) {
       if ((y_wfs(n).npix*y_wfs(n)._nrebin) % 2 != y_wfs(n)._Nfft % 2)
       y_wfs(n)._halfxy = &float(*y_wfs(n)._halfxy-2.*halfxy);
       if (y_wfs(n)._Nfft % 2 == 0)
       y_wfs(n)._halfxy = &float(*y_wfs(n)._halfxy+halfxy);
       }
       */
      if ((y_wfs(n).npix % 2 == 1) && (y_wfs(n)._nrebin %2 == 1))
      y_wfs(n)._halfxy = &float(halfxy*0.0f);
      else y_wfs(n)._halfxy = &float(halfxy);
      //  } else {
      //if (y_wfs(n).npix % 2 < y_wfs(n)._pdiam % 2) y_wfs(n)._halfxy = &float(halfxy);
      //else y_wfs(n)._halfxy = &float(halfxy*0.);
      //}
    } else {
      halfxy = (span(0,2*pi,y_wfs(n)._pdiam+1)(1:y_wfs(n)._pdiam) / 2.)(,-:1:y_wfs(n)._pdiam);
      y_wfs(n)._halfxy = &float(halfxy*0.0f);
    }
  }

  if (y_wfs(n).type == "sh") {
    //this defines how we create a larger fov if required
    if (y_wfs(n)._Ntot != y_wfs(n)._Nfft) {
      x1=long((y_wfs(n)._Ntot-y_wfs(n)._Nfft)/2.)+1;
      x2=long(x1+y_wfs(n)._Nfft-1);

      tmp = indices(y_wfs(n)._Nfft);
      hrpix = array(0.,y_wfs(n)._Ntot,y_wfs(n)._Ntot);
      hrpix(x1:x2,x1:x2) = roll(tmp(,,1)+(tmp(,,2)-1)*y_wfs(n)._Nfft);
      hrmap =where(roll(hrpix));

      y_wfs(n)._hrmap = &int(hrmap-1);
    } else y_wfs(n)._hrmap = &int(0);

    if (y_wfs(n)._nrebin*y_wfs(n).npix % 2 != y_wfs(n)._Ntot % 2)
    x1=long((y_wfs(n)._Ntot-y_wfs(n)._nrebin*y_wfs(n).npix)/2.)+2;
    else x1=long((y_wfs(n)._Ntot-y_wfs(n)._nrebin*y_wfs(n).npix)/2.)+1;
    x2=long(x1+y_wfs(n)._nrebin*y_wfs(n).npix-1);

    binindices = array(0,y_wfs(n)._Ntot,y_wfs(n)._Ntot);
    tmp = long((indices(y_wfs(n)._nrebin*y_wfs(n).npix) -1) / y_wfs(n)._nrebin);
    binindices(x1:x2,x1:x2) = tmp(,,1)+tmp(,,2)*y_wfs(n).npix+1;

    binmap = array(0,y_wfs(n)._nrebin*y_wfs(n)._nrebin,y_wfs(n).npix*y_wfs(n).npix);
    tmp = indices(y_wfs(n)._Ntot)-1;
    tmp = tmp(,,1)+tmp(,,2)*y_wfs(n)._Ntot;
    for (i=1;i<=y_wfs(n).npix*y_wfs(n).npix;i++) {
      if (y_wfs(n).gsalt > 0) binmap(,i) = tmp(where(binindices == i));
      else
      binmap(,i) = tmp(where(roll(binindices) == i));
    }
    y_wfs(n)._binmap = &int(binmap);

    /* verif
     fim = array(0.0f,y_wfs(n).npix,y_wfs(n).npix);
     for(i=1;i<=y_wfs(n).npix*y_wfs(n).npix;i++) fim(*)(i) = res(,,7)(*)(binmap(,i)+1)(sum);
     */

    dr0 = y_tel.diam/y_atmos.r0*(0.5/y_wfs(n).lambda)^1.2/cos(y_geom.zenithangle*dtor)^0.6;
    fwhmseeing = y_wfs(n).lambda/
    (y_tel.diam/sqrt(y_wfs(n).nxsub^2.+(dr0/1.5)^2.))/4.848;
    kernelfwhm = sqrt(fwhmseeing^2.+y_wfs(n).kernel^2.);

    sdim = y_wfs(n)._Ntot;
    tmp = eclat(makegaussian(sdim,kernelfwhm/y_wfs(n)._qpixsize,
            xc=sdim/2+1,yc=sdim/2+1));

    tmp(1,1) = 1.; // this insures that even with fwhm=0, the kernel is a dirac
    tmp = tmp/sum(tmp);
    tmp2 = fft(tmp)/sdim/sdim;

    tmp = array(float,2,sdim,sdim);
    tmp(1,,) = tmp2.re;
    tmp(2,,) = tmp2.im;
    y_wfs(n)._ftkernel = &(tmp);

    // dealing with photometry
    telSurf = pi/4.*(1-y_tel.cobs^2.)*y_tel.diam^2.;

    // from the guide star 
    if (y_wfs(n).gsalt == 0) {
      if (y_wfs(n).zerop == 0) y_wfs(n).zerop = 1e11;
      y_wfs(n)._nphotons = y_wfs(n).zerop*10^(-0.4*y_wfs(n).gsmag)*
      y_wfs(n).optthroughput*                 // include throughput to WFS
      (y_tel.diam/y_wfs(n).nxsub)^2./telSurf*// for unobstructed subaperture
      y_loop.ittime;// per iteration
    } else {  // we are dealing with a LGS
      y_wfs(n)._nphotons = y_wfs(n).lgsreturnperwatt*// detected by WFS
      y_wfs(n).laserpower*// ... for given power
      y_wfs(n).optthroughput*// include throughput to WFS
      (y_tel.diam/y_wfs(n).nxsub)^2.*1e4*// for unobstructed subaperture
      y_loop.ittime;// per iteration
    }
    write,format="nphotons : %.1f\n",y_wfs(n)._nphotons;
  }
}

func prep_lgs_prof(numwfs,prof,h,beam,center=,imat=)
/* DOCUMENT prep_lgs_prof
 prep_lgs_prof,numwfs,prof,h,beam,center=

 The function returns an image array(double,n,n) of a laser beacon elongated by perpective
 effect. It is obtaind by convolution of a gaussian of width "lgsWidth" arcseconds, with the
 line of the sodium profile "prof". The altitude of the profile is the array "h".
 prof     : Na profile intensity, in arbitrary units
 h        : altitude, in meters. h MUST be an array with EQUALLY spaced elements.
 beam     : size in arcsec of the laser beam
 center   : string, either "image" or "fourier" depending on where the centre should be.
 
 Computation of LGS spot from the sodium profile:
 Everything is done here in 1D, because the Na profile is the result of the convolution of a function
 P(x,y) = profile(x) . dirac(y)
 by a gaussian function, for which variables x and y can be split :
 exp(-(x^2+y^2)/2.s^2)  =  exp(-x^2/2.s^2) * exp(-y^2/2.s^2)
 The convolution is (symbol $ denotes integral)
 C(X,Y) = $$ exp(-x^2/2.s^2) * exp(-y^2/2.s^2) * profile(x-X) * dirac(y-Y)  dx  dy
 First one performs the integration along y
 C(X,Y) = exp(-Y^2/2.s^2)  $ exp(-x^2/2.s^2) * profile(x-X)  dx
 which shows that the profile can be computed by
 - convolving the 1-D profile
 - multiplying it in the 2nd dimension by a gaussian function
 
 If one has to undersample the inital profile, then some structures may be "lost". In this case,
 it's better to try to "save" those structures by re-sampling the integral of the profile, and
 then derivating it afterwards.
 Now, if the initial profile is a coarse one, and that one has to oversample it, then a
 simple re-sampling of the profile is adequate.
 
 SEE ALSO:
 */
{
  extern y_wfs;

  y_wfs(numwfs)._prof1d = &float(prof);
  y_wfs(numwfs)._profcum = &float(prof(cum));

  subapdiam = y_tel.diam / double(y_wfs(numwfs).nxsub); // diam of subap
  if (y_wfs(numwfs).nxsub > 1)
  xsubs = span(-y_tel.diam/2+subapdiam/2,y_tel.diam/2-subapdiam/2,y_wfs(numwfs).nxsub);
  else xsubs = 0.;
  ysubs = xsubs;

  np = dimsof(prof)(2);// number of points of the profile
  hG = sum(h*prof)/sum(prof);// center of gravity of the profile
  x=(indgen(y_wfs(numwfs)._Ntot)-(y_wfs(numwfs)._Ntot/2.+1));
  // x expressed in pixels. (0,0) is in the fourier-center.
  x = x*y_wfs(numwfs)._qpixsize;// x expressed in arcseconds
  dx = x(2)-x(1);
  dh = h(2)-h(1);

  if (y_wfs(numwfs).nxsub > 1)
  dOffAxis = sqrt((xsubs((*y_wfs(numwfs)._validsubs)(1,))-y_wfs(numwfs).lltx)^2 +
      (ysubs((*y_wfs(numwfs)._validsubs)(2,))-y_wfs(numwfs).llty)^2);
  else
  dOffAxis = sqrt((xsubs-y_wfs(numwfs).lltx)^2 +
      (ysubs-y_wfs(numwfs).llty)^2);

  if (imat) dOffAxis *= 0.0f;

  w = beam / 2.35482005;//  sigma.   (2*sqrt(2*log(2)))=2.35482005
  if( w==0 ) {
    g = array(0.0,n);
    if( center=="image" )
    g(n/2:n/2+1)=0.5;
    else
    g(n/2+1)=1;
  }
  else {
    if( center=="image" )
    g = exp( -((x+y_wfs(numwfs)._qpixsize/2)^2/(2*w^2.) ) );
    else
    g = exp( -(x^2/(2*w^2.) ) );
  }

  y_wfs(numwfs)._ftbeam = &float(transpose([fft(g,1).re,fft(g,1).im]));
  y_wfs(numwfs)._beam = &float(g);
  // convolved profile in 1D.

  if (numberof(xsubs) > 1)
  azimut=atan(ysubs((*y_wfs(numwfs)._validsubs)(2,))-y_wfs(numwfs).llty,
      xsubs((*y_wfs(numwfs)._validsubs)(1,))-y_wfs(numwfs).lltx);
  else
  azimut=atan(ysubs-y_wfs(numwfs).llty,xsubs-y_wfs(numwfs).lltx);

  y_wfs(numwfs)._azimuth = &float(azimut);

  sensors_initlgs,g_wfs,numwfs-1,numberof(*y_wfs(numwfs)._prof1d),hG,h(1),dh,y_wfs(numwfs)._qpixsize,
  dOffAxis,float(*y_wfs(numwfs)._prof1d),float(*y_wfs(numwfs)._profcum),float(*y_wfs(numwfs)._beam),
  float(*y_wfs(numwfs)._ftbeam),float(*y_wfs(numwfs)._azimuth);

}

func wfs_map(arr,wfs,type=)
/* DOCUMENT wfs_map
 wfs_map,arr,wfs,type=

 maps an array of images onto a wfs
 arr     : the array to map
 wfs     : the wfs on which to map
 type    : type of mapping
 
 SEE ALSO:
 */
{
  if (type == []) type = "subaps"
  if (type == "subaps") {
    if (numberof(arr) != (*wfs._isvalid)(*)(sum))
    error,"wfs_map : wrong input dims";
    tmp = array(structof(arr),wfs.nxsub,wfs.nxsub);
    tmp(where(*wfs._isvalid)) = arr;
    return tmp;
  }
  if (type == "image") {
    if (dimsof(arr)(1) != 3)
    error,"wfs_map : wrong input dims";
    sz = dimsof(arr)(2);
    tmp = array(structof(arr),wfs.nxsub*sz,wfs.nxsub*sz);
    for (cc=1;cc<=wfs._nvalid;cc++) {
      indi = ((*wfs._validsubs)(1,cc)-1)*sz+1;
      indj = ((*wfs._validsubs)(2,cc)-1)*sz+1;
      tmp(indi:indi+sz-1,indj:indj+sz-1) = arr(,,cc);
    }
    return tmp;
  }
}

func make_lgs_prof1d(numwfs,prof,h,beam,center=)
/* DOCUMENT make_lgs_prof1d
 make_lgs_prof1d,numwfs,prof,h,beam,center=

 same as prep_lgs_prof but cpu only. original routine from rico
 
 SEE ALSO:
 */
{
  extern y_wfs;

  y_wfs(numwfs)._prof1d = &float(prof);
  y_wfs(numwfs)._profcum = &float(prof(cum));

  subapdiam = y_tel.diam / double(y_wfs(numwfs).nxsub); // diam of subap
  if (y_wfs(numwfs).nxsub > 1)
  xsubs = span(-y_tel.diam/2+subapdiam/2,y_tel.diam/2-subapdiam/2,y_wfs(numwfs).nxsub);
  else xsubs = 0.;
  ysubs = xsubs;

  np = dimsof(prof)(2);// number of points of the profile
  hG = sum(h*prof)/sum(prof);// center of gravity of the profile
  x=(indgen(y_wfs(numwfs)._Ntot)-(y_wfs(numwfs)._Ntot/2.+1));
  // x expressed in pixels. (0,0) is in the fourier-center.
  x = x*y_wfs(numwfs)._qpixsize;// x expressed in arcseconds
  dx = x(2)-x(1);
  dh = h(2)-h(1);

  if (y_wfs(numwfs).nxsub > 1)
  dOffAxis = sqrt((xsubs((*y_wfs(numwfs)._validsubs)(1,))-y_wfs(numwfs).lltx)^2 +
      (ysubs((*y_wfs(numwfs)._validsubs)(2,))-y_wfs(numwfs).llty)^2);
  else
  dOffAxis = sqrt((xsubs-y_wfs(numwfs).lltx)^2 +
      (ysubs-y_wfs(numwfs).llty)^2);

  profi = array(0.0,y_wfs(numwfs)._Ntot,y_wfs(numwfs)._nvalid);
  subsdone = array(1,y_wfs(numwfs)._nvalid);
  dif2do = array(0,y_wfs(numwfs)._nvalid);

  while (subsdone(*)(sum) > 0) {
    tmp = dOffAxis(where(subsdone)(1));
    inds = where(dOffAxis == tmp);
    // height, translated in arcsec due to perspective effect
    zhc = (h-hG)*(206265.*tmp/double(hG)^2.);

    //x=(indgen(y_wfs(numwfs)._Ntot)-(y_wfs(numwfs)._Ntot/2.+0.5));
    //x = x*y_wfs(numwfs)._qpixsize;  // x expressed in arcseconds
    //xtmp = x / (206265./double(hG)^2.)/tmp + hG;
    //xind = (xtmp-h(1))/deltah+1;

    dzhc = zhc(2)-zhc(1);
    if (y_wfs(numwfs)._qpixsize > dzhc) {
      profi(,inds) = interp( prof(cum), zhc(pcen), x(pcen) )(dif);
    } else {
      profi(,inds) = interp( prof, zhc, x );
    }
    subsdone(inds) = 0;
  }

  //profi /= profi(max,)(-,);

  w = beam / 2.35482005;//  sigma.   (2*sqrt(2*log(2)))=2.35482005
  if( w==0 ) {
    g = array(0.0,n);
    if( center=="image" )
    g(n/2:n/2+1)=0.5;
    else
    g(n/2+1)=1;
  }
  else {
    if( center=="image" ) {
      if ((y_wfs(numwfs).npix*y_wfs(numwfs)._nrebin) % 2 != y_wfs(numwfs)._Nfft % 2)
      g = exp( -((x+y_wfs(numwfs)._qpixsize)^2/(2*w^2.)));
      else
      g = exp( -((x+y_wfs(numwfs)._qpixsize/2)^2/(2*w^2.)));
    } else
    g = exp( -(x^2/(2*w^2.) ) );
  }

  y_wfs(numwfs)._ftbeam = &float(transpose([fft(g,1).re,fft(g,1).im]));
  y_wfs(numwfs)._beam = &float(g);
  // convolved profile in 1D.

  p1d = roll(fft(fft(profi,[1,0]) * fft(g(,-:1:y_wfs(numwfs)._nvalid),[1,0]),[-1,0]).re,
      [long(y_wfs(numwfs)._Ntot/2.+0.5),0]);
  // abs is just to ensure only positive values (else values ~ -1e-12 may appear)
  p1d = abs(p1d);
  im = p1d(,-,) * g(-,,-:1:y_wfs(numwfs)._nvalid);

  if (numberof(ysubs) > 1)
  azimut=atan(ysubs((*y_wfs(numwfs)._validsubs)(2,))-y_wfs(numwfs).llty,
      xsubs((*y_wfs(numwfs)._validsubs)(1,))-y_wfs(numwfs).lltx);
  else
  azimut=atan(ysubs-y_wfs(numwfs).llty,xsubs-y_wfs(numwfs).lltx);

  y_wfs(numwfs)._azimuth = &float(azimut);

  if (center == "image")
  xcent = ycent = y_wfs(numwfs)._Ntot/2+0.5;
  else
  xcent = ycent = y_wfs(numwfs)._Ntot/2+1;

  if (numberof(ysubs) > 1) {
    im = rotate3d(im,azimut*180./pi,xcent,ycent);
    tmp = im(*,)(max,);
    im /= tmp(-,-,);
  } else {
    im = rotate(im,azimut*180./pi,xcent,ycent);
    im /= max(im);
  }
  y_wfs(numwfs)._lgskern = &float(im);
}

func make_lgs_prof1d_slow(numwfs,prof,h,beam,center=)
/* DOCUMENT make_lgs_prof1d
 make_lgs_prof1d,numwfs,prof,h,beam,center=

 same as prep_lgs_prof but cpu only. original routine from rico
 
 SEE ALSO:
 */
{
  extern y_wfs;

  subapdiam = y_tel.diam / double(y_wfs(numwfs).nxsub); // diam of subap
  if (y_wfs(numwfs).nxsub > 1)
  xsubs = span(-y_tel.diam/2+subapdiam/2,y_tel.diam/2-subapdiam/2,y_wfs(numwfs).nxsub);
  else xsubs = 0.;
  ysubs = xsubs;

  np = dimsof(prof)(2);// number of points of the profile
  hG = sum(h*prof)/sum(prof);// center of gravity of the profile
  x=(indgen(y_wfs(numwfs)._Ntot)-(y_wfs(numwfs)._Ntot/2.+0.5));
  // x expressed in pixels. (0,0) is in the fourier-center.
  x = x*y_wfs(numwfs)._qpixsize;// x expressed in arcseconds
  dx = x(2)-x(1);
  dh = h(2)-h(1);

  if (y_wfs(numwfs).nxsub > 1)
  dOffAxis = sqrt((xsubs((*y_wfs(numwfs)._validsubs)(1,))-y_wfs(numwfs).lltx)^2 +
      (ysubs((*y_wfs(numwfs)._validsubs)(2,))-y_wfs(numwfs).llty)^2);
  else
  dOffAxis = sqrt((xsubs-y_wfs(numwfs).lltx)^2 +
      (ysubs-y_wfs(numwfs).llty)^2);

  azimut=atan(ysubs((*y_wfs(numwfs)._validsubs)(2,))-y_wfs(numwfs).llty,
      xsubs((*y_wfs(numwfs)._validsubs)(1,))-y_wfs(numwfs).lltx);

  y_wfs(numwfs)._azimuth = &float(azimut);

  tmp = array(float,[3,y_wfs(numwfs)._Ntot,y_wfs(numwfs)._Ntot,y_wfs(numwfs)._nvalid]);
  for (i=1;i<=y_wfs(numwfs)._nvalid;i++) tmp(,,i) = makeLgsProfile1D(y_wfs(numwfs)._Ntot, y_wfs(numwfs)._qpixsize, beam, prof, h, dOffAxis(i), hG, azimut(i)*180./pi,center=center);
  y_wfs(numwfs)._lgskern = &float(tmp);
}

func compare_spots(void) {
  //#include "yorick/yoga_wfs.i"

  profilename = "allProfileNa_withAltitude_1Gaussian.fits";
  prof = fits_read(YOGA_AO_SAVEPATH + profilename);
  h = prof(,1);
  prof = prof(,2:)(,avg);
  make_lgs_prof1d, 1, prof, h, y_wfs(1).beamsize;
  tmp = (*y_wfs(1)._lgskern);
  make_lgs_prof1d_slow, 1, prof, h, y_wfs(1).beamsize;
  tmp2 = (*y_wfs(1)._lgskern);

}

func makeLgsProfile1D(n, z, lgsWidth, prof, h, dOffAxis, H, azimut, center=)
/* DOCUMENT

 The function returns an image array(double,n,n) of a laser beacon elongated by perpective
 effect. It is obtaind by convolution of a gaussian of width "lgsWidth" arcseconds, with the
 line of the sodium profile "prof". The altitude of the profile is the array "h".
 n        : image size
 z        : pixel size of the image (arcsec)
 lgsWidth : width of gaussian, in arcsec
 prof     : Na profile intensity, in arbitrary units
 h        : altitude, in meters. h MUST be an array with EQUALLY spaced elements.
 dOffAxis : offaxis distance of the subaperture that sees the beacon
 H        : altitude at which the beacon is supposed to be (in meters)
 azimut   : rotation of the beacon
 center   : string, either "image" or "fourier" depending on where the centre should be.
 
 Computation of LGS spot from the sodium profile:
 Everything is done here in 1D, because the Na profile is the result of the convolution of a function
 P(x,y) = profile(x) . dirac(y)
 by a gaussian function, for which variables x and y can be split :
 exp(-(x^2+y^2)/2.s^2)  =  exp(-x^2/2.s^2) * exp(-y^2/2.s^2)
 The convolution is (symbol $ denotes integral)
 C(X,Y) = $$ exp(-x^2/2.s^2) * exp(-y^2/2.s^2) * profile(x-X) * dirac(y-Y)  dx  dy
 First one performs the integration along y
 C(X,Y) = exp(-Y^2/2.s^2)  $ exp(-x^2/2.s^2) * profile(x-X)  dx
 which shows that the profile can be computed by
 - convolving the 1-D profile
 - multiplying it in the 2nd dimension by a gaussian function
 
 SEE ALSO:
 */
{
  //prof =  exp(-(altitude-9.e4)^2/(2*(4.e3)^2))*135.; // mono-modal
  // exp(-(altitude-8.7e4)^2/(2*(2.e3)^2))*110+exp(-(altitude-9.2e4)^2/(2*(2.e3)^2))*130; // bi-modal
  // exp(-(altitude-8.5e4)^2/(2*(1.5e3)^2))*55+exp(-(altitude-8.8e4)^2/(2*(2.e3)^2))*80+exp(-(altitude-9.2e4)^2/(2*(2.e3)^2))*120; // tri-modal
  // exp(-(altitude-8.7e4)^2/(2*(2.e3)^2))*130+exp(-(altitude-9.2e4)^2/(2*(2.e3)^2))*130; //bi-modal sym
  // note : the fwhm is 2*sqrt(2*log(2))*sigma ...

  np = dimsof(prof)(2);// number of points of the profile
  hG = sum(h*prof)/sum(prof);// center of gravity of the profile, expressed as an index
  // transformation de h en arcsec
  zhc = (h-hG)*(206265.*dOffAxis/double(H)^2.);// height, translated in arcsec due to perspective effect

  x=(indgen(n)-(n/2+1));// x expressed in pixels. (0,0) is in the fourier-center.
  x = x*z;// x expressed in arcseconds

  /* If one has to undersample the inital profile, then some structures may be "lost". In this case,
   it's better to try to "save" those structures by re-sampling the integral of the profile, and
   then derivating it afterwards.
   Now, if the initial profile is a coarse one, and that one has to oversample it, then a
   simple re-sampling of the profile is adequate.
   */
  dzhc = zhc(2)-zhc(1);   // sampling period of initial profile
  if( z>dzhc ) {
    aprof = interp( prof(cum), zhc(pcen), x(pcen) )(dif); // resampled profile in 1D
  } else {
    aprof = interp( prof, zhc, x );   // resampled profile in 1D
  }
  //aprof=interp(aprof,x,x+sum(x*aprof)/sum(aprof));

  dec = 0.0;
  if( center=="image" ) dec = z/2;// dec will allow to shift the image by z/2 if image-centered is required

  w = lgsWidth / 2.35482005;//  sigma.   (2*sqrt(2*log(2)))=2.35482005
  if( w==0 ) {
    g = array(0.0,n);
    if( center=="image" )
    g(n/2:n/2+1)=0.5;
    else
    g(n/2+1)=1;
  }
  else {
    if( center=="image" )
    g = exp( -((x+z/2)^2/(2*w^2.) ) );
    else
    g = exp( -(x^2/(2*w^2.) ) );
    //=mygauss2(24,13,13,wx*2*sqrt(2*log(2)),wy*2*sqrt(2*log(2)),1.,0.,0.,,deriv=0);
  }
  // convolved profile in 1D.
  p1d = roll(fft(fft(aprof) * fft(g),-1).re);
  // abs is just to ensure only positive values (else values ~ -1e-12 may appear)
  p1d = abs(p1d);
  im = p1d(,-) * g(-,);

  if (center == "image") {
    im = rotate(im,azimut,n/2+0.5,n/2+0.5);
  } else {
    im = rotate(im,azimut,n/2+1,n/2+1);
  }

  return im/max(im);
}

/*
 //pyr check
 ampli=sensors_getdata(g_wfs,0,"amplifoc")
 tmp2 = array(complex,304,304)
 tmp2.re = ampli(1,,)
 tmp2.im = ampli(2,,)
 ca = roll(tmp2,[108+14,108+34])
 q1 = ca(1:108,1:108)
 ca = roll(tmp2,[0+14,108+34])
 q2 = ca(1:108,1:108)
 ca = roll(tmp2,[108+14,0+34])
 q3 = ca(1:108,1:108)
 ca = roll(tmp2,[0+14,0+34])
 q4 = ca(1:108,1:108)

 pli,q1.re
 pli,sensors_getdata(g_wfs,0,"fttotim")(1,,,1)

 */

func noise_cov(ns)
// Compute the noise covariance matrix of the WFS
    {
  extern y_wfs,y_tel,y_dm,y_atmos;
  if (y_wfs(ns).noise >= 0) {
    // Photon noise
    validsub = *y_wfs(ns)._validsubs;
    m = validsub(2,);
    n = validsub(1,);
    ind = (m - 1) * y_wfs(ns).nxsub + n; // Index of valid subaperture
    flux = (*y_wfs(ns)._fluxPerSub)(*)(ind); // Flux for each subaperture
    Nph = flux * y_wfs(ns)._nphotons;  // Nphotons on each subaperture

    //patchDiam = y_tel.diam+2*max(abs([y_wfs(ns).xpos,y_wfs(ns).ypos]))*4.848e-6*abs(y_dm(ndm).alt);
    //d = patchDiam / (y_dm(ndm).nact-1); // Interactuator distance
    //d = 1.;
    r0 = get_r0(y_atmos.r0, 0.5, y_wfs(ns).lambda);

    sig = (pi ^ 2 / 2) * (1 / Nph) * (1. / r0) ^ 2; // Photon noise in m⁻²
    sig = sig / (2 * pi / (y_wfs(ns).lambda * 1e-6)) ^ 2; // Noise variance in rad²
    sig *= RASC ^ 2;

    //error;

    //Electronic noise
    Ns = y_wfs(ns).npix; // Number of pixel
    //Nd = Ns/3. ; // FWHM of the image in pixels (Nyquist condition ?)
    Nd = (y_wfs(ns).lambda * 1e-6) * RASC / y_wfs(ns).pixsize;
    sigphi = (pi ^ 2 / 3) * (1 / Nph ^ 2) * (y_wfs(ns).noise) ^ 2 * Ns ^ 2 * (Ns / Nd) ^ 2; // Phase variance in m⁻²
    sigsh = sigphi / (2 * pi / (y_wfs(ns).lambda * 1e-6)) ^ 2; // Noise variance in rad²
    sigsh *= RASC ^ 2; // Electronic noise variance in arcsec²

    //sig *= 0; // No photon noise for debug

    cov_n = array(0.0f, 2 * numberof(sig));
  cov_n(:numberof(sig)) = sig + sigsh;
  cov_n(numberof(sig)+1:) = sig + sigsh;
} else
  cov_n = array(0.0f, 2 * dimsof(*y_wfs(ns)._validsubs)(3));

return cov_n;
}
