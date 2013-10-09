/*
  Collection of routines for building the KL of a Kolmogorov statistics
  Derived from a collection of routines developped both at ESO and ONERA

*/

struct kl_basis_struct
{
  long    nr;            
  long    ni;  
  long    np;       
  long    nfunc;   
  float   cobs;  
  long    nord;   
  pointer radp;   
  pointer evals;   
  pointer npo;   
  pointer ord;   
  pointer rabas;   
  pointer azbas;
  /////////////////////////////////////////////////////////
  long    ncp;   
  long    ncmar;   
  pointer px;            
  pointer py;  
  pointer cr;       
  pointer cp;   
  pointer pincx;  
  pointer pincy;   
  pointer pincw;   
  pointer ap;   
};

func kolstf(dvec)
  /*DOCUMENT var=kolstf(dvec)
  This routine returns the kolmogorov phase variance at spatial
  dimension (inverse of the spatial frequency) dvec
   */
{
  return 6.88 * dvec^(5./3.);
}

func karmanstf(dvec,outscl=)
  /*DOCUMENT var=kolstf(dvec)
  This routine returns the Von Karman phase variance at spatial
  dimension (inverse of the spatial frequency) dvec. Same as kolstf
  but with a correcting factor to account for the outter scale.
  The latter should be in units of telescope diameter
   */
{
  if (dimsof(outscl)==[]) outscl = 3.;
  return 6.88 * dvec^(5./3.)*(1-1.485*(dvec/outscl)^(1./3.)+\
                              5.383*(dvec/outscl)^(2)-6.281*\
                              (dvec/outscl)^(7./3.));
}

func dblindgen(n)
  /*DOCUMENT res=dblindgen(size)

  D. Gratadour Feb 2006
  
  This routine returns a size x size array containing increasing indices
  from 1 to size x size.

  SEE ALSO : indgen, indices
   */
{
  n=long(n);
  return reform(indgen(n*n),[2,n,n]);
}

func radii(nr,np,cobs)
  /*DOCUMENT res=radii(NumberOfR,NumberOfPhi,Dim)
  This routine generates an nr x np array with np copies of the
  radial coordinate array. Radial coordinate span the range from
  r=cobs to r=1 with successive annuli having equal areas (ie, the
  area between cobs and 1 is divided into nr equal rings, and the
  points are positioned at the half-area mark on each ring). There
  are no points on the border.     
   */
{
  r2 = cobs^2 +(float(indgen(nr)-1)+0.)/nr*(1.0 - cobs^2);
  rs = sqrt(r2);
  return rs*array(1.,np)(-,);
}

func polang(r)
  /*DOCUMENT res=polang(RadialCoordArray)
  This routine generates an array with the same dimensions as r,
  but containing the azimuthal values for a polar coordinate system.     
   */
{
  s =  dimsof(r);
  nr = s(2);
  np = s(3);
  phi1 = float(indgen(np)-1)/float(np)*2.*pi;
  return phi1(-,)*array(1.,nr);
}

func setpincs(ax,ay,px,py,cobs,&pincx,&pincy,&pincw)
  /*DOCUMENT setpincs(ax,ay,px,py,cobs,&pincx,&pincy,&pincw)
  This routine determines a set of squares for interpolating
  from cartesian to polar coordinates, using only those points
  with cobs < r < 1     

  SEE ALSO : pcgeom
   */
{
  s = dimsof(ax);
  nc = s(2);
  s = dimsof(px);
  nr = s(2);
  np = s(3);
  dcar = (ax(nc) - ax(1)) / (nc-1);
  ofcar = ax(1,1);

  rlx = (px - ofcar)/dcar;
  rly = (py - ofcar)/dcar;
  lx = long(rlx);
  ly = long(rly);
  shx = rlx - lx;
  shy = rly - ly;

  pincx=[lx,lx+1,lx+1,lx]+1;
  pincy=[ly,ly,ly+1,ly+1]+1;
  
  pincw=[(1-shx)*(1-shy),shx*(1-shy),shx*shy,(1-shx)*shy];
  
  axy = ax^2 + ay^2;
  axyinap = clip(axy,cobs^2.+1.e-3,0.999);
  sizeaxyinap=(dimsof(axyinap))(2);
  pincw = pincw*axyinap(pincx+(pincy-1)*sizeaxyinap);
  pincw = pincw*(1.0/pincw(,,sum))(,,-);
}

func pcgeom(&geom,ncp,ncmar)    
  /*DOCUMENT pcgeom,&geom,ncp,ncmar
  This routine builds a geom_struct. px and py are the x and y
  coordinates of points in the polar arrays.  cr and cp are the
  r and phi coordinates of points in the cartesian grids. ncmar
  allows the possibility that there is a margin of ncmar points
  in the cartesian arays outside the region of interest
   */
{
  nr   = geom.nr;
  np   = geom.np;
  cobs = geom.cobs;
  
  nused = ncp - 2*ncmar;
  ff = 0.5 * nused;
  hw =  float(ncp-1)/2;
    
  r = radii(nr,np,cobs); 
  p = polang(r);

  px0 = r * cos(p);
  py0 = r * sin(p);
  px = ff * px0 + hw;
  py = ff * py0 + hw;
  ax = float(dblindgen(ncp)-1) % ncp - 0.5 * (ncp-1);
  ax = ax / (0.5 * nused); 
  ay = transpose(ax);
        
  setpincs, ax, ay, px0, py0, cobs,pincx, pincy, pincw;
  dpi = 2 * pi;
  cr2 = (ax^2 + ay^2); 
  ap = clip(cr2,cobs^2+1.e-3,0.999);
  //cr = (cr2 - cobs^2) / (1 - cobs^2) * nr - 0.5; 
  cr = (cr2 - cobs^2) / (1 - cobs^2) * nr; 
  cp = (atan(ay, ax) + dpi) % dpi;
  cp = (np / dpi) * cp;
    
  cr = clip(cr,1.e-3,nr-1.001);
  //fudge -----, but one of the less bad ones
  cp = clip(cp,1.e-3,np -1.001);
  //fudge -----  this is the line which
  //gives that step in the cartesian grid
  //at phi = 0.
  
  geom.ncp   = ncp;
  geom.ncmar = ncmar; 
  geom.px    = &px;
  geom.py    = &py; 
  geom.cr    = &cr;
  geom.cp    = &cp;
  geom.pincx = &pincx;
  geom.pincy = &pincy;
  geom.pincw = &pincw;
  geom.ap    = &ap;
}

func set_pctr(&bas, ncp =, ncmar=)
  /*DOCUMENT geom=set_pctr(bas, ncp =, ncmar=)
  This routine calls pcgeom to build a geom_struct with the
  right initializations. bas is a gkl_basis_struct built with
  the gkl_bas routine.
  */
{
  if (!is_set(ncmar)) ncmar = 2;
  if (!is_set(ncp)) ncp = 128;
    
  pcgeom,bas,ncp,ncmar;
}

func make_radii(nr,cobs)
  /*DOCUMENT res=radii(nr,cobs)
  This routine generates an nr x np array with np copies of the
  radial coordinate array. Radial coordinate span the range from
  r=cobs to r=1 with successive annuli having equal areas (ie, the
  area between cobs and 1 is divided into nr equal rings, and the
  points are positioned at the half-area mark on each ring). There
  are no points on the border.     
   */
{
  d = (1.-cobs*cobs)/nr;
  rad2 = cobs^2 +d/16.+ d * float(indgen(nr)-1);
  rad = sqrt(rad2);
  
  return rad; 
}

func make_kernels(cobs,nr,rad,funct=,outscl=)
  /*DOCUMENT res=make_kernels(cobs,nr,rad,funct=,outscl=)
  This routine generates the kernel used to find the KL modes.
  The  kernel constructed here should be simply a discretization
  of the continuous kernel. It needs rescaling before it is treated
  as a matrix for finding  the eigen-values. The outter scale
  should be in units of the diameter of the telescope.
   */
{
  nth = 5*nr;
  kers  = array(float,[3,nr, nr, nth]);
  cth = cos(float(indgen(nth)-1)*(2.*pi/nth));
  dth = 2.*pi/nth;
  fnorm = -1./(2*pi*(1.-cobs^2))*0.5;
  //the 0.5 is to give  the r^2 kernel, not
  //the r kernel
  for (i =1;i<=nr;i++) { 
    for (j=1;j<=i;j++) {
      te = 0.5*sqrt(rad(i)^2+rad(j)^2-(2*rad(i)*rad(j))*cth);
      //te in units of the diameter, not the radius
      if (funct=="kolmo") te = kolstf(te);
      if (funct=="karman") te = karmanstf(te,outscl=outscl);
      if ((funct!="kolmo") & (funct!="karman")) {
        write,"The statistics is not known !";
        error;
      }
      kelt =  fnorm * dth * float (fft(te,-1));
      kers (i, j,) = kers (j, i,) = kelt;
    }
  }
  return kers;
}

func piston_orth(nr)
  /*DOCUMENT piston_orth(nr)
   */
{
  s = array(float,[2,nr,nr]);
  for (j=1;j<=nr-1;j++) {
    rnm = 1./sqrt (float(j*(j+1)));
    s(1:j,j) = rnm;
    s(j+1,j)= -1*j*rnm;
  }
  rnm = 1./sqrt (nr);
  s(,nr) = rnm;
  return s;
}


func gkl_fcom(kers,cobs,nf,&evals,&nord,&npo,&ord,&rabas)
  /*DOCUMENT gkl_fcom(kers,cobs,nf,&evals,&nord,&npo,&ord,&rabas)
  This routine does the work : finding the eigenvalues and
  corresponding eigenvectors. Sort them and select the right
  one. It returns the KL modes : in polar coordinates : rabas
  as well as the associated variance : evals. It also returns
  a bunch of indices used to recover the modes in cartesian
  coordinates (nord, npo and ord).
   */
{
  s = dimsof(kers);
  nr = s(2);
  nt = s(4);
  nxt = 1;
  fktom =  (1.-cobs^2)/nr;
  fevtos = sqrt(2*nr);
  
  evs = array(float,[2,nr,nt]);
  //ff isnt used - the normalisation for
  //the eigenvectors is straightforward:
  //integral of surface^2 divided by area = 1,
  //and the cos^2 term gives a factor
  //half, so multiply zero order by
  //sqrt(n) and the rest by sqrt (2n)

  //zero order is a special case...
  //need to deflate to eliminate infinite eigenvalue - actually want
  //evals/evecs of zom - b where b is big and negative
  zom = kers(,,1);
  s = piston_orth(nr);
  ts =transpose(s);
  b1 = ((ts(,+)*zom(+,))(,+)*s(+,))(1:nr-1, 1:nr-1);
 
  newev = SVdec(fktom*b1,v0,vt);

  v1 = array(float,[2,nr, nr]);
  v1(1:nr-1,1:nr-1) = v0;
  v1(nr,nr) = 1;

  vs = s(,+)*v1(+,);
  grow,newev,0;
  evs(,nxt) = float(newev);
  kers (,, nxt) = sqrt(nr)*vs;

  nxt = 2;
  do {
    newev = SVdec(fktom*kers(,,nxt),vs,vt);
    evs(,nxt) = float(newev);
    kers (,,nxt) = sqrt(2.*nr)*vs;
    mxn = max(float(newev));
    egtmxn = floor(evs(, 1:nxt)>mxn);
    nxt = nxt + 1;
  } while ((2*sum(egtmxn)-sum(egtmxn(,1))) < nf);
  nus = nxt - 1;

  kers = kers (,,1:nus);
  evs = reform (evs (, 1:nus), nr*nus);
  a = (sort(-1.*evs))(1:nf);
  //every eigenvalue occurs twice except
  //those for the zeroth order mode. This
  //could be done without the loops, but
  //it isn't the stricking point anyway...
  no = 1;
  ni = 1;
  oind = array(long,nf+1);
  do {
       if (a(ni) < nr+1) {
         oind(no) = a(ni);
         no = no + 1;
       } else {
         oind(no) = a(ni);
         oind(no+1) = a(ni);
         no = no + 2;
       }
       ni = ni + 1;
  } while (no < (nf+1));
  
  oind = oind (1:nf);
  tord = (oind-1)/nr+1;
  odd = ((long(indgen(nf)-1) % 2) == 1);
  pio = (oind-1) % nr +1;

  evals = evs(oind);
  ord = 2 *(tord-1) - floor(tord>1 & (odd))+1;

  nord = max(ord);
  rabas = array(float,[2,nr, nf]);
  sizenpo=long(max(ord));
  npo = array(long,sizenpo);
  
  for (i=1;i<=nf;i++) {
    npo(long(ord(i))) = npo(long(ord(i))) + 1;
    rabas(, i) = kers (, pio(i), tord(i));
  }
}

func make_azimuth(nord, np)
  /*DOCUMENT piston_orth(nr)
   */
{
  azi = array(float,[2,long(1+nord), np]);
  th = float(indgen(np)-1)*(2.*pi/ np);

  azi (1,) = 1.0;
  for (i = 2; i<=nord;i+=2)  azi (i,) = cos (((i-1)/2+1) * th);
  for (i = 3; i<=nord;i+=2)  azi (i,) = sin (((i-1)/2) * th);
  return azi;
}

func make_klbas(nfunc,cobs,dim,nr=,np=,funct=,outscl=)
  /*DOCUMENT make_klbas(nfunc,cobs,nr=,np=,funct=,outscl=)
  SEE ALSO : 
   */
{
  if (cobs == []) cobs = 0;
  if (nfunc == []) nfunc = 500L;
  
  if (!is_set(nr)) nr = long(5.0f*sqrt(nfunc));
  if (!is_set(np)) np = long(5*nr);
  if (dimsof(funct)==[]) funct="kolmo";
  
  radp = make_radii(nr, cobs);

  kers = make_kernels(cobs, nr, radp,funct=funct,outscl=outscl);
  
  gkl_fcom,kers,cobs,nfunc,evals,nord,npo,ord,rabas;

  azbas = make_azimuth(nord, np);

  klbasis = kl_basis_struct();
  klbasis.nr=nr;
  klbasis.np=np;
  klbasis.nfunc=nfunc; 
  klbasis.cobs=cobs;
  klbasis.radp=&radp;
  klbasis.evals=&evals;
  klbasis.nord=nord;
  klbasis.npo=&npo;
  klbasis.ord=&ord;
  klbasis.rabas=&rabas;
  klbasis.azbas=&transpose(azbas);
  
  set_pctr,klbasis, ncp= dim;
  
  return klbasis;
}

func gkl_sfi(bas, i)
  /*DOCUMENT 
  This routine returns the i'th function from the generalised KL
  basis bas. bas must be generated first with gkl_bas.
   */
{    
  if (i>bas.nfunc) { 
    write, "the basis only contains ", nfunc, "functions";
    return 0;
  }
  /*
  nr   = bas.nr;
  np   = bas.np;
  ordp = *bas.ord;
  ord  = long(ordp(i));

  rabas = (*bas.rabas)(,i);

  azbas = (*bas.azbas)(,ord);

  sf1   = rabas(,-:1:np);

  sf2   = azbas(-:1:nr,);

  sf    = sf1 * sf2;
  */
  nr = bas.nr;
  np = bas.np;
  ordp = *bas.ord;
  ord=long(ordp(i));

  rabasp=*bas.rabas;
  rabas=rabasp(,i);

  azbasp=transpose(*bas.azbas);
  azbas=azbasp(ord, );

  sf1=array(double,[2,nr,np]);
  sf1(,*)=rabas;

  sf2=array(float,[2,np,nr]);
  sf2(,*)=azbas;  

  sf = sf1*transpose(sf2);

  return sf;
}




func pol2car(cpgeom,pol,mask=)
  /*DOCUMENT cart=pol2car(cpgeom, pol, mask=)
  This routine is used for polar to cartesian conversion.
  pol is built with gkl_bas and cpgeom with pcgeom.
  However, points not in the aperture are actually treated
  as though they were at the first or last radial polar value
  -- a small fudge, but not serious  ?*******
   */
{
  cd = bilinear(pol, *cpgeom.cr+1, *cpgeom.cp+1);
  if (mask!=[]) cd = cd*(*cpgeom.ap);
  return cd;
} 

