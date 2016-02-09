#include"fits-utils.i"

dtor = pi/180.;
radeg = 1./dtor;

RASC = 180*3600/pi;

//

func rotscale(A,alpha,M)
/*
    This function reads the array A and proceed with its Fourier transform (matrix-based, not an FFT).
    It then rotate the resulting array with the rotate2 routine, and computes a second FT with a new size M.
*/
{
    OWA=256;
    E=qft(A,1.,OWA,int(4*OWA));
    E_rot=rotate2(re_part(E),alpha)+1i*rotate2(im_part(E),alpha);
    B=qft(E_rot,OWA,1.,int(M));
    return B;
}

func qft(A,D1,D2,M)
/*
 
 This function computes the Fourier transform E of the field A.
 - D1: diameter of A
 - D2: diameter of E (in lambda/D)
 - M: number of points for E (along one axis)
 
 The FT is computed using array products.
 
 It returns E, the complex electric field.
 
 */
{
	//extern N,x,X,u,U,E,I;
	
	if (int(sqrt(numberof(A)))==0) {
		exit, "A must be a square array!"
	}
	
	N=int(sqrt(numberof(A)));
	
	//o=array(1,N);
	
	//x=span(-D1*0.5+1./(4*N),D1*0.5-1./(4*N), N);
	
	x=VECT(D1,N);
	X=MESH(D1,N);
    
	Y=transpose(X);
	
	//o=array(1,M);
	//u=span(-D2*0.5+1./(4*M),D2*0.5-1./(4*M), M);
	//U=o(-,)*u(,-);
	
	u=VECT(D2,M);
	U=MESH(D2,M);
    
	V=transpose(U);
	
	XU=x(,-)*u(-,);
	YV=transpose(XU);
	
	E1=exp(2*pi*1i*XU);
	E2=exp(2*pi*1i*YV);
	
	Ex=A(,+)*E1(+,)*D1/N;
	E=E2(,+)*Ex(+,)*D1/N;
	
	//I=abs(E)^2;
	//I=log10(I/max(I(,)));
    
    return E;
}

func VECT(D,N)
{
	o=array(1,N);
	x=span(-D*0.5+1./(4*N),D*0.5-1./(4*N), N);
	return x;
}

func MESH(D,N)
{
	o=array(1,N);
	x=span(-D*0.5+1./(4*N),D*0.5-1./(4*N), N);
	X=o(-,)*x(,-);
	return X;
}

func pad_array(A,N)
{
	L=max(dimsof(A));
    D=(N-L)/2;
    O=array(0.0,N,N);
    O(D+1:N-D,D+1:N-D)=A;
    A=O;
    return A;
}

//

func regressAx(y,x)
/* DOCUMENT coef = regressAx(y,x)
   returns the regression coefficient A, when fitting
   the data y by the function y=A.x
   This function expects that (y=Ax+B with B=0). This makes the linear fit more
   robust, since the regression line is constrained to go through (0,0).
 */
{
  return sum(y*x)/sum(x*x);
}

func fft_goodsize(size)
  /* DOCUMENT func fft_goodsize(size)
     find best size for a fft (radix 2, 3, 5 or 7) from size
   */
{
  aradix=[3,5,7];
  mradix = 2;
  /*
  for (i=1;i<=3;i++) {
    tmpf = log(size)/log(aradix(i));
    tmpl = long(tmpf);
    tmpf -= tmpl;
    mradix = ((tmpf > (log(size)/log(mradix) - long(log(size)/log(mradix))))(1) ? aradix(i) : mradix);
  }
  */
  return mradix^(long(log(size)/log(mradix))+1);
}

func rotate3d(im,ang,cx,cy,zoom)
/* DOCUMENT imrot = rotate(im,ang,cx,cy,zoom)
            imrot = rotate(im,ang)

Rotates an image of an angle "ang" (in DEGREES).
The center of rotation is cx,cy.
A zoom factor can be applied.

(cx,cy) can be omitted :one will assume one rotates around the
center of the image.
If zoom is not specified, the default value of 1.0 is taken.

modif dg : allow to rotate a cube of images with one angle per image

 */
{
  if( is_void(zoom) ) zoom=1.0;
  if( zoom==0.0 ) zoom=1.0;
  if (numberof(ang) == 1)
    if(ang==0 & zoom==1.0) return im;
  
  ang *= pi/180.0;
  s = dimsof(im);
  nx = s(2);
  ny = s(3);
  if( is_void(cx) )
    cx = nx/2 + 1;
  if( is_void(cy) )
    cy = ny/2 + 1;
  x = indgen(1:nx)(,-:1:ny) - cx;
  y = indgen(1:ny)(-:1:nx,) - cy;
  x /= zoom;
  y /= zoom;
  
  if (numberof(ang)>1) {
    rxy = array(0.,nx,ny,2,numberof(ang));
    for (i=1;i<=numberof(ang);i++) {
      matrot = [[cos(ang(i)),-sin(ang(i))],[sin(ang(i)),cos(ang(i))]];
      rxy(,,,i) = [x,y](,,+) * matrot(,+) + [cx,cy](-,-,);
    }
  } else {
    matrot = [[cos(ang),-sin(ang)],[sin(ang),cos(ang)]];
    rxy = [x,y](,,+) * matrot(,+) + [cx,cy](-,-,);
  }
  
  nn = where(rxy<1);
  if( is_array(nn) )
    rxy(nn)=1.;
  
  if (numberof(ang) > 1) {
    rx = rxy(,,1,);
    ry = rxy(,,2,);
  } else {
    rx = rxy(,,1);
    ry = rxy(,,2);
  }
  nn = where(rx>(nx-1));
  if( is_array(nn) )
    rx(nn)=nx-1;
  nn = where(ry>(ny-1));
  if( is_array(nn) )
    ry(nn)=ny-1;

  wx = rx;
  wy = ry;
  rx = long(rx);   // partie entiere
  ry = long(ry);
  wx -= rx;        // partie fractionnaire
  wy -= ry;

  ind = rx + (ry-1)*nx;
  
  if (numberof(ang) > 1) {
    nim = indgen(numberof(ang))-1;
    nim *= (nx*ny);
    nim = nim(-:1:nx,-:1:ny,);
    ind += nim;
  }
  
  imr = (im(ind)*(1-wx) + im(ind+1)*wx)*(1-wy) + (im(ind+nx)*(1-wx)+im(ind+nx+1)*wx)*wy;
  return imr;
}

func rotate(im,ang,cx,cy,zoom)
/* DOCUMENT imrot = rotate(im,ang,cx,cy,zoom)
            imrot = rotate(im,ang)

Rotates an image of an angle "ang" (in DEGREES).
The center of rotation is cx,cy.
A zoom factor can be applied.

(cx,cy) can be omitted :one will assume one rotates around the
center of the image.
If zoom is not specified, the default value of 1.0 is taken.

 */
{
  if( is_void(zoom) ) zoom=1.0;
  if( zoom==0.0 ) zoom=1.0;
  if(ang==0 & zoom==1.0) return im;
  ang *= pi/180.0;
  s = dimsof(im);
  nx = s(2);
  ny = s(3);
  if( is_void(cx) )
    cx = nx/2 + 1;
  if( is_void(cy) )
    cy = ny/2 + 1;
  x = indgen(1:nx)(,-:1:ny) - cx;
  y = indgen(1:ny)(-:1:nx,) - cy;
  x /= zoom;
  y /= zoom;
  
  matrot = [[cos(ang),-sin(ang)],[sin(ang),cos(ang)]];
  rxy = [x,y](,,+) * matrot(,+) + [cx,cy](-,-,);
  nn = where(rxy<1);
  if( is_array(nn) )
    rxy(nn)=1.;
  rx = rxy(,,1);
  ry = rxy(,,2);
  nn = where(rx>(nx-1));
  if( is_array(nn) )
    rx(nn)=nx-1;
  nn = where(ry>(ny-1));
  if( is_array(nn) )
    ry(nn)=ny-1;

  wx = rx;
  wy = ry;
  rx = long(rx);   // partie entiere
  ry = long(ry);
  wx -= rx;        // partie fractionnaire
  wy -= ry;

  ind = rx + (ry-1)*nx;
  imr = (im(ind)*(1-wx) + im(ind+1)*wx)*(1-wy) + (im(ind+nx)*(1-wx)+im(ind+nx+1)*wx)*wy;
  return imr;
}


func fft_rotate(im,angle,xc=,yc=,gband=)
/* DOCUMENT fft_rotate(im,angle)
   im    : square image
   angle : rotation angle in degrees (clockwise)
   xc,yc : x,y positions of the rotation center
   
   high precision image rotation using fft
   no information loss if : image is shannon sampled
                            image has sufficiently large guard band
   using the algorithm from :
   "Fast Fourier method for the accurate rotation of sampled images"
   Kieran G. Larkin, Michael A. Oldfield, Hanno Klemm
   Optics Communications 139 (1997) 99-106

   routine by d. gratadour 31/05/2011
   SEE ALSO:
 */

{
  if (angle == 0.) return im;
  size = dimsof(im);
  if (size(2) != size(3)) error,"works only on square images";
  nx = size(2);

  if (angle >= 360) angle = angle % 360;
  if (angle > 180) return fft_rotate(im,angle-360,xc=xc,yc=yc,gband=gband);
  if (angle < -180) return fft_rotate(im,360+angle,xc=xc,yc=yc,gband=gband);

  if (gband != 0) {
    im2=array(double,2*nx,2*nx);
    if (nx %2 == 0)
      im2(nx-nx/2+1:nx+nx/2,nx-nx/2+1:nx+nx/2) = im;
    else
      im2(nx-nx/2:nx+nx/2,nx-nx/2:nx+nx/2) = im;      
    im = im2;
    nx *= 2;
  }
  
  if (xc == []) xc = ((nx/2)%2 == 0 ? nx/2+1 : nx/2);
  if (yc == []) yc = ((nx/2)%2 == 0 ? nx/2+1 : nx/2);

  theta = angle * pi/180.;
  
  if (angle > 90) theta = pi/2;
  if (angle < -90) theta = -pi/2;

  stepx = tan(theta/2);
  stepy = -1.*sin(theta);

  if ((nx/2)%2 == 0) {
    tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-0.5));
    tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-0.5));
  } else {
    tiltx = -2.*pi/nx*(float(stepx)*(indgen(nx)-nx/2-1));
    tilty = -2.*pi/nx*(float(stepy)*(indgen(nx)-nx/2-1));
  }
  
  compltiltx=array(complex,nx);
  compltiltx.im=roll(tiltx);
  compltiltx = compltiltx(,-::nx-1);
  
  compltilty=array(complex,nx);
  compltilty.im=roll(tilty);
  compltilty = compltilty(-::nx-1,);
  
  col  = span(1,nx,nx)(,-:1:nx);
  lig  = transpose(col);
    
  tmpc=array(complex,nx,nx);

  tmpc = fft(exp(compltiltx*(lig-xc))*fft(im,[1,0]),[-1,0]);
  tmpc = fft(exp(compltilty*(col-yc))*fft(tmpc,[0,1]),[0,-1]);
  tmpc = fft(exp(compltiltx*(lig-xc))*fft(tmpc,[1,0]),[-1,0]);

  if (angle > 90) {
    return fft_rotate(tmpc.re/nx/nx/nx,angle-90.,gband=0);
  } else {
    if (angle < -90) {
      return fft_rotate(tmpc.re/nx/nx/nx,angle+90.,gband=0);
    } else {
      if ((nx/2)%2 == 0)
        return tmpc(nx/2-nx/4+1:nx/2+nx/4,nx/2-nx/4+1:nx/2+nx/4).re/nx/nx/nx;
      else
        return tmpc(nx/2-nx/4:nx/2+nx/4,nx/2-nx/4:nx/2+nx/4).re/nx/nx/nx;
    }
  }
}

func calc_dphi(phase,pup,den)
/* DOCUMENT
 * 
 */ 
{
  npix = dimsof(phase)(2);
  mi = p = dphi = array(float,[2,2*npix,2*npix]);
  mi(1:npix,1:npix) = phase;
  p(1:npix,1:npix)  = pup;

  //den  = (fft(fft(p,1)*conj(fft(p,1)),-1)).re;
  mask = den > max(den)*1.e-7;
  pix  = where(mask);
  //tmp = fft( abs(fft(mi*p))^2. ).re;
  //tmp(pix) /= (den)(pix);
  //dphi(pix) = tmp(1,1) - tmp(pix);
  dphi(pix) = fft(2*((fft(mi^2*p,1)*conj(fft(p,1))).re - abs(fft(mi*p,1))^2),-1).re(pix)/den(pix);

  return dphi; 
}

func calc_dphis(phase)
/* DOCUMENT
 * 
 */ 
{
  npix = dimsof(phase)(2);
  mi = p = dphi = array(float,[2,2*npix,2*npix]);
  mi(1:npix,1:npix) = phase;
  dphi = (mi - mi(1,1))^2;

  return dphi; 
}

func calc_dphif(phase,pup,den,conjftpup)
/* DOCUMENT 
 *  
 * KEYWORDS :  
 */ 
{
	npix = dimsof(phase)(2);
	mi = p = dphi = array(float,[2,2*npix,2*npix]);
	mi(1:npix,1:npix) = phase;
	p(1:npix,1:npix)  = pup;

	mask = den > max(den)*1.e-7;
	pix  = where(mask);
    //dphi(pix) = fft(2.*((fft(mi^2*p,1)*conjftpup) - abs(fft(mi*p,1))^2).re,-1).re(pix)/den(pix);
	dphi(pix) = fft(fft(mi^2*p, 1)*conjftpup + fft(p, 1)*conj(fft(mi^2*p, 1)) -2.*fft(mi*p, 1)*conj(fft(mi*p, 1)), -1).re(pix)/den(pix)
  
	return dphi;   
}

func circavg(a,center=,middle=)
/* DOCUMENT circavg
 *  
 * average=circavg(array[,center=,middle=])
 *
 * This routine returns the circular mean of an array. The routine dist is
 * used to compute an array of coordinate of points to be averaged.
 *
 * a      (input) : The array of which we want the circular mean. It can be
 *                  a long, a float or a double. Complex are not supported
 * center (input) : An array of the form [x,y] which give the position of
 *                  the origin for the circular mean calculation
 * middle (input) : A flag to indicate that the origin for the circular mean
 *                  calculation is the middle of the array a
 *
 * SEE ALSO: 
 */ 
{
  s=dimsof(a);

  if (!is_set(middle)) middle=0;
  if (s(1) != 2) write,"error - invalid dimensions";
  if (s(3) != s(2)) write,"error - invalid dimensions";

  dim=s(2);
  
  if (center!=[]) {
    s=dimsof(center);
    if ((s(1) != 1) | (s(1) != 2)) \
       write,"error - center has invalid dimensions";

    center=long(center);

    if (middle) {
      center=long([0,0]);
      write,"error - center and middle are not compatible keywords";
    }
  } else { 
    if (middle) center=long([0,0]);
    else center=[dim/2,dim/2];
  }
  
  r=long(roll(long(dist(dim)+.5)+1,[center(1),center(2)]));
  j=long(max(r));
  n=array(long,j);
  sx=array(double,j);
  dim2=long(dim)*long(dim);
  
  for (i=1;i<=dim2;i++) {
    j=r(i);
    sx(j)=sx(j)+a(i);
    n(j)=n(j)+1;
  }
  
  return sx/n;
}


func circavg_quad(a)
/* DOCUMENT circavg_quad
 *  
 * average=circavg_quad(array)
 *
 *
 * SEE ALSO: 
 */ 
{
  s=dimsof(a);

  if (s(1) != 2) write,"error - invalid dimensions";
  if (s(3) != s(2)) write,"error - invalid dimensions";

  dim=s(2);
  
  
  r=long(roll(dist(2*dim)+.5)+1)(1:dim,1:dim);
  j=long(max(r));
  n=array(long,j);
  sx=array(double,j);
  dim2=long(dim)*long(dim);
  
  for (i=1;i<=dim2;i++) {
    j=r(i);
    sx(j)=sx(j)+a(i);
    n(j)=n(j)+1;
  }
  
  return sx/n;
}

/*
 _____                     __   __ _    ___  
|  ___| __ ___  _ __ ___   \ \ / // \  / _ \ 
| |_ | '__/ _ \| '_ ` _ \   \ V // _ \| | | |
|  _|| | | (_) | | | | | |   | |/ ___ \ |_| |
|_|  |_|  \___/|_| |_| |_|   |_/_/   \_\___/ 
                                             
 */
func __mysinc(ar)
/* DOCUMENT func sinc(ar)
 * Return the sinus cardinal of the input array
 * F.Rigaut, 2002/04/03
 * SEE ALSO: Eric Thiebault wrote a sinc which is probably better.
 */
{
  local ar;
  ar = double(ar);
  w  = where(abs(ar) < 1e-10);
  if (exist(w)) {ar(w) = 1e-10;}
  return sin(ar)/ar;
}

func make_pupil_generic(dim,pupd,xc=,yc=,real=,cobs=,t_spiders=,spiders_type=,angle=)
  /* DOCUMENT func make_pupil(dim,pupd,xc=,yc=,real=)
   */
{
  if (real == 1) {
    pup = exp(-(dist(dim,xc=xc,yc=yc)/(pupd/2.))^60.)^0.69314;
  } else {
    //    xc;yc;info,xc;
    //    tv,dist(dim,xc=xc,yc=yc);pause,2000;
    pup = dist(dim,xc=xc,yc=yc) < (pupd+1.)/2.;
  }
  if (is_set(cobs)) {
    if (real == 1) {
      pup -= exp(-(dist(dim,xc=xc,yc=yc)/(pupd*cobs/2.))^60.)^0.69314;
    } else {
      pup -= dist(dim,xc=xc,yc=yc) < (pupd*cobs+1.)/2.;

      X=span(-0.5*(1-1./dim),0.5*(1-1./dim),dim)(,-:1:dim);
      Y=transpose(X);
        
      if (!is_set(spiders_type)) spiders_type="six";
      if (!is_set(t_spiders)) t_spiders=0.01;

      t_spiders=t_spiders*pupd/dim;
        
      if (spiders_type=="four") {

          spiders_map_1 = ((Y > (X+t_spiders/(2*sin(pi/4)))*tan(pi/4))+(Y < (X-t_spiders/(2*sin(pi/4)))*tan(pi/4)));
          spiders_map_2 = ((Y > (-X+t_spiders/(2*sin(pi/4)))*tan(pi/4))+(Y < (-X-t_spiders/(2*sin(pi/4)))*tan(pi/4)));
          pup = pup*spiders_map_1*spiders_map_2;

      } else if (spiders_type=="six") {
      
          angle = 180/15.;
          spiders_map_1 = ((Y > (-X+t_spiders/(2*sin(pi/2-pi/angle)))*tan(pi/2-pi/angle))+(Y < (-X-t_spiders/(2*sin(pi/2-pi/angle)))*tan(pi/2-pi/angle)));
          spiders_map_2 = ((Y > (X+t_spiders/(2*sin(pi/6+pi/angle)))*tan(pi/6+pi/angle))+(Y < (X-t_spiders/(2*sin(pi/6+pi/angle)))*tan(pi/6+pi/angle)));
          spiders_map_3 = ((Y > (-X+t_spiders/(2*sin(pi/6-pi/angle)))*tan(pi/6-pi/angle))+(Y < (-X-t_spiders/(2*sin(pi/6-pi/angle)))*tan(pi/6-pi/angle)));
          pup = pup*spiders_map_1*spiders_map_2*spiders_map_3;

      }
    }
  }
    
  return pup;
}

func make_VLT(dim,pupd,cobs,t_spiders)
/* DOCUMENT func make_VLT(dim,pupd,cobs,t_spiders)
 */
{
    
    diam=8.0;
    y_tel.diam=diam;
    
    cobs=0.14;
    y_tel.cobs=cobs;
    
    t_spiders=0.005;
    y_tel.t_spiders=t_spiders;
    
    X=MESH(1.,dim);
    Y=transpose(X);
    R=sqrt(X^2.+Y^2.);
    
    off=2.52/36;
    t_spiders=0.09/18;
    angle = 50.5*pi/180;
    
    pup=(R<0.5)*(R>cobs/2);
    
    spiders_map = ((Y > (X-cobs/2+t_spiders/sin(angle))*tan(angle))+(Y < (X-cobs/2)*tan(angle)))*(X>0)*(Y>0);
    spiders_map = spiders_map+spiders_map(dim:1:-1,);
    spiders_map = spiders_map+spiders_map(,dim:1:-1);
    
    pup = pup*spiders_map;
    
    return pup;
}

// EELT pupil building ; to be called by main pupil building routine (make_pupil).
func make_EELT(dim,pupd,type_ap,cobs,N_seg,nbr_miss_seg,std_ref_err,t_spiders,angle)
{
    if (fileExist(swrite(format=YOGA_AO_SAVEPATH+"apertures/"+"%s_N%d_COBS%d_CLOCKED%d_TSPIDERS%2.1f_MS%d_REFERR%d.fits",type_ap,dim,long(100*cobs),long(angle),100*t_spiders,nbr_miss_seg,long(100*std_ref_err))))
    {
    
        pup = fits_read(swrite(format=YOGA_AO_SAVEPATH+"apertures/"+"%s_N%d_COBS%d_CLOCKED%d_TSPIDERS%2.1f_MS%d_REFERR%d.fits",type_ap,dim,long(100*cobs),long(angle),100*t_spiders,nbr_miss_seg,long(100*std_ref_err)));
    
    } else {
    
        f=open(YOGA_AO_SAVEPATH+"apertures/Coord_"+type_ap+".dat","r");
        x_seg= y_seg= array(0.0,N_seg);
        read, f, x_seg, y_seg;
        close, f;
        
        f2=open(YOGA_AO_SAVEPATH+"apertures/EELT_MISSING_"+type_ap+".dat","r");
        k_seg = array(0.0,1,20);
        read, f2, k_seg;
        close, f2;
        
        W=1.45*cos(pi/6);
        D=40.;
        y_tel.diam=D;
        X=MESH(1.*D*dim/pupd,dim);
        Y=transpose(X);
        
        if (nbr_miss_seg>0) k_seg=k_seg(1:nbr_miss_seg);
        k_seg=int(k_seg(sort(k_seg)));
        
        f3=open(YOGA_AO_SAVEPATH+"apertures/EELT_REF_ERROR.dat","r");
        ref_err = array(0.0,798);
        read, f3, ref_err;
        close, f3;
        mean_ref = sum(ref_err)/798;
        std_ref = sqrt(1./798.*sum((ref_err-mean_ref)^2));
        ref_err = ref_err * std_ref_err/std_ref;
        if (nbr_miss_seg>0) ref_err(k_seg) = 1.0;
        
        pup=array(0.0,dim,dim);
        for (i=1 ; i<= N_seg ; i++)
        {
            Xt=X+x_seg(i);
            Yt=Y+y_seg(i);
            pup=pup+(1-ref_err(i))*(Yt<0.5*W)*(Yt>=-0.5*W)*(0.5*(Yt+tan(pi/3.)*Xt)<0.5*W)*(0.5*(Yt+tan(pi/3.)*Xt)>=-0.5*W)*(0.5*(Yt-tan(pi/3.)*Xt)<0.5*W)*(0.5*(Yt-tan(pi/3.)*Xt)>=-0.5*W);
        }
        t_spiders_ori=t_spiders;
        t_spiders = t_spiders*D*dim/pupd;
        spiders_map_1 = abs(X) > t_spiders/2;
        spiders_map_2 = ((Y > (X+t_spiders/(2*sin(pi/6)))*tan(pi/6))+(Y < (X-t_spiders/(2*sin(pi/6)))*tan(pi/6)));
        spiders_map_3 = ((Y > (-X+t_spiders/(2*sin(pi/6)))*tan(pi/6))+(Y < (-X-t_spiders/(2*sin(pi/6)))*tan(pi/6)));
        t_spiders = t_spiders_ori;
        pup = pup*spiders_map_1*spiders_map_2*spiders_map_3;
        if (angle != 0) pup=rotate2(pup,angle);
        "EELT pupil has been created.";
        
        fits_write, swrite(format=YOGA_AO_SAVEPATH+"apertures/"+"%s_N%d_COBS%d_CLOCKED%d_TSPIDERS%2.1f_MS%d_REFERR%d.fits",type_ap,dim,long(100*cobs),long(angle),100*t_spiders_ori,nbr_miss_seg,long(100*std_ref_err)), pup, overwrite=1;
        "New aperture saved";
        
    }
    
    return pup;
}

// Main pupil building routine.
func make_pupil(dim,pupd,type_ap=,angle=,spiders_type=,t_spiders=,nbr_miss_seg=,std_ref_err=,xc=,yc=,real=,cobs=)
{
    if (!is_set(angle)) angle=0.;
    if (!is_set(nbr_miss_seg)) nbr_miss_seg=0;
    if (!is_set(std_ref_err)) std_ref_err=0.;

    if (type_ap=="EELT-Nominal") {
        
        N_seg=798;
        cobs = 0.3;
        y_tel.cobs = cobs;
        pup=make_EELT(dim,pupd,type_ap,cobs,N_seg,nbr_miss_seg,std_ref_err,t_spiders,angle);
        
    } else if (type_ap=="EELT-BP1") {
        
        N_seg=768;
        cobs = 0.369;
        y_tel.cobs = cobs;
        pup=make_EELT(dim,pupd,type_ap,cobs,N_seg,nbr_miss_seg,std_ref_err,t_spiders,angle);
        
    } else if (type_ap=="EELT-BP3") {

        N_seg=672;
        cobs = 0.503;
        y_tel.cobs = cobs;
        pup=make_EELT(dim,pupd,type_ap,cobs,N_seg,nbr_miss_seg,std_ref_err,t_spiders,angle);
        
    } else if (type_ap=="EELT-BP5") {

        N_seg=558;
        cobs = 0.632;
        y_tel.cobs = cobs;
        pup=make_EELT(dim,pupd,type_ap,cobs,N_seg,nbr_miss_seg,std_ref_err,t_spiders,angle);
        
    } else if (type_ap=="VLT") {
        
        pup=make_VLT(dim,pupd);
        
    } else {
        
        type_ap="Generic";
        y_tel.type_ap=type_ap;
        pup=make_pupil_generic(dim,pupd,xc=xc,yc=yc,real=real,cobs=cobs,t_spiders=t_spiders,spiders_type=spiders_type,angle=angle);
    }
    
    return pup;
}

func make_apodizer(dim,pupd,filename,angle)
/* DOCUMENT func make_apodizer(dim,pupd,filename)
 */
{
    "Opening apodizer";
    pup=fits_read(filename);
    
    A=dimsof(pup);
    A=A(2);
    
    if (A>dim) {
        error, "Apodizer dimensions must be smaller.";
    }
    
    if (A != pupd) {
        pup=bilinear(pup,pupd,pupd);
    }
    
    if (angle != 0) {
        pup=rotate2(pup,angle);
    }
    
    //pup=rotscale(pup,angle,pupd);
    
    reg=where(dist(pupd)>pupd/2.);
    pup(reg)=0.;
    
    pupf=array(0.,dim,dim);
    
    if (dim != pupd) {
        if ((dim-pupd)%2 != 0) {
            pupf((dim-pupd+1)/2+1:(dim+pupd+1)/2,(dim-pupd+1)/2+1:(dim+pupd+1)/2)=pup;
        
        } else {
            pupf((dim-pupd)/2+1:(dim+pupd)/2,(dim-pupd)/2+1:(dim+pupd)/2)=pup;
        }
    } else {
        pupf=pup;
    }
    
    pupf=abs(pupf);
    
    return pupf;
}

func plvf(vy,vx,y,x,autoscale=,scale=,width=,hsize=,hang=,color=,type=,prop=)
/* DOCUMENT plvf,vy,vx,y,x,scale=,width=,hsize=,hang=,color=,type=,prop=
   Plots the vector field defined by (vx,vy) at positions (x,y)
   vx,vy,x,y must have the same size, but can be of arbitrary dimension.
   KEYWORDS:
   autoscale: set to 1 for the vector length to be autoscaled
   scale:     multiplicative factor applied to the autoscale results
              (for fine tweaking)
   width, color, type: same as in plg.
   hsize, hang: size and opening angle of the arrow head (default
       hsize=0.4, hang=20 degrees)
   prop:      set to zero if you want the same head size for *all* vector.
              Otherwise, the head size is proportionnal to the size of
              the vector (which results in something nicer to the eye).
   SEE ALSO: pldj
 */
{
  if (!scale) scale=1.;
  if (!width) width=2;
  if (hsize==[]) hsize=0.4;
  if (hang==[]) hang = 20;
  if (prop==[]) prop = 1;

  if (autoscale) {  
    if (prop) {
      sc=abs(vx,vy);
      if (max(sc)==0) sc=1.;
      //      else sc=sc/max(sc);
    } else {sc=1.;}

    // vector body autoscaling:
    xdif = abs(x(dif));
    w = where(xdif != 0);
    if (numberof(w)!=0) {
      minspace = min(xdif(w));
    }
    ydif = abs(y(dif));
    w = where(ydif != 0);
    if (numberof(w)!=0) {
      minspace = (minspace==[]? min(ydif(w)) : min([minspace,min(ydif(w))]) );
    }
    if (minspace==[]) minspace=1.;
    // autoscale normalization factor: max vector length / min space between location
    norm = max([vy,vx])/minspace*1.2;
    if (norm==0) norm=1.;
    vx = vx/norm*scale;
    vy = vy/norm*scale;
    //    hsize = hsize/norm*scale;
  } else {
  }
  sc = abs(vx,vy);

  pldj,(x+vx)(*),(y+vy)(*),x(*),y(*),width=width,color=color,type=type;
  x1=(x+vx)(*);  y1=(y+vy)(*);
  ang=atan(vy(*),vx(*))-(180-hang)*pi/180.;
  x2=x1+sc*hsize*cos(ang);
  y2=y1+sc*hsize*sin(ang);
  pldj,x2,y2,x1,y1,width=width,color=color,type=type;

  ang=atan(vy,vx)-(180+hang)*pi/180.;
  x2=x1+sc*hsize*cos(ang);
  y2=y1+sc*hsize*sin(ang);
  pldj,x2,y2,x1,y1,width=width,color=color,type=type;
}

func mydist(dim,xc=,yc=)
/* DOCUMENT func dist(dim,xc=,yc=)
 * Return an array which elements are the distance to (xc,yc). xc and
 * yc can be omitted, in which case they are defaulted to size/2+1.
 * F.Rigaut, 2003/12/10.
 * SEE ALSO:
*/

{
  dim = long(dim);

  if (is_scalar(dim)) dim=[2,dim,dim];
  if ((is_vector(dim)) && (dim(1)!=2))
    error,"Dist only deals with 2D square arrays";

  d = array(float,dim);

  if (xc!=[]) {xc = float(xc-1.);} else {xc = float(dim(2)/2);}
  if (yc!=[]) {yc = float(yc-1.);} else {yc = float(dim(3)/2);}

  res = _dist(&d,dim(2),dim(3),xc,yc);

  return d;
}

//---------------------------------------------------------

func factoriel(arg)
/* DOCUMENT factoriel(arg)
 * Return factoriel of the argument    
 * SEE ALSO:
 */
{
  if (arg == 0) {
    return 1.;
  } else {
    res = 1.;
    for (i=1;i<=arg;i++) res = res*i;
    return res;
   }
}

//---------------------------------------------------------

func zernumero(zn)
/* DOCUMENT zernumero(zn)
 * Returns the radial degree and the azimuthal number of zernike
 * number zn, according to Noll numbering (Noll, JOSA, 1976)
 * SEE ALSO: prepzernike, zernike
 */
{
  j	= 0;
  for (n=0;n<=100;n++)
   {
    for (m=0;m<=n;m++)
     {
      if (even(n-m))
       {
        j	= j+1;
        if (j == zn) {return [n,m];}
        if (m != 0)
         {
          j	= j+1;
          if (j == zn) {return [n,m];}
         }
       }
     }
   }
}

//---------------------------------------------------------

func make_zernike(nzer,size,diameter,xc,yc,ext)
/* DOCUMENT make_zernike(nzer,size,diameter,xc,yc)
 * Call this function to set up an array of zernike polynomials.
 * size : size of the 2d array on which future "zernike" will be returned
 * diameter : diameter of the pupil in pixel in the array
 * xc, yc (optional) : Coordinates (in pixels of the center of the pupil)
 * These zernikes follow the Noll (JOSA, 1976) numbering and
 * definition (rms of 1 over the pupil)
 * SEE ALSO: zernike,zernike_ext,zernumero
 */
{

  if (xc == []) xc = size/2+1;
  if (yc == []) yc = size/2+1;
  if (ext == []) ext = 0;

  radius= (diameter+1.)/2.;
  zdim	= size;
  zr	= mydist(zdim,xc=xc,yc=yc)/radius;
  zmask	= (zr <= 1.);
  zmaskmod = (zr <= 1.2);
  zrmod	= zr*zmaskmod;
  zr	= zr*zmask;
  x	= float(span(1,zdim,zdim)(,-:1:zdim));
  y	= transpose(x);
  zteta	= atan(y-yc,x-xc);

  z	= array(float,zdim,zdim,nzer);
  for (zn=1;zn<=nzer;zn++) {
    znm	= zernumero(zn) ; n=znm(1) ; m=znm(2);

    if (ext) {
      for (i=0;i<=(n-m)/2;i++) {
        z(,,zn) = z(,,zn) + (-1.)^i*zrmod^(n-2.*i)*factoriel(n-i)/
          (factoriel(i)*factoriel((n+m)/2-i)*factoriel((n-m)/2-i));
      }
    } else {
      for (i=0;i<=(n-m)/2;i++) {
        z(,,zn) = z(,,zn) + (-1.)^i*zr^(n-2.*i)*factoriel(n-i)/
          (factoriel(i)*factoriel((n+m)/2-i)*factoriel((n-m)/2-i));
      }
    }

    if (odd(zn)) {
      if (m == 0) {
        z(,,zn) = z(,,zn)*sqrt(n+1.);
      } else {
        z(,,zn) = z(,,zn)*sqrt(2*(n+1.))*sin(m*zteta);
      }
    } else {
      if (m == 0) {
        z(,,zn) = z(,,zn)*sqrt(n+1.);
      } else {
        z(,,zn) = z(,,zn)*sqrt(2*(n+1.))*cos(m*zteta);
      }
    }

  }
  
  if (ext) return z*zmaskmod(,,-:1:nezr);
  else return z*zmask(,,-:1:nezr);
}

func get_pyrimg(mimg)
{
  npix = dimsof(mimg)(2);
  tmp = array(0.,[2,2*npix+3,2*npix+3]);
  ii1 = 2:npix+1;
  ii2 = npix+3:2*npix+2;
  tmp(ii1,ii1) = mimg(,,1);
  tmp(ii2,ii1) = mimg(,,2);
  tmp(ii1,ii2) = mimg(,,3);
  tmp(ii2,ii2) = mimg(,,4);
  return tmp;
}

func get_r0(r0_at_lambda1, lambda1, lambda2)
/* DOCUMENT get_r0(r0_at_lambda1, lambda1, lambda2)
 * Give the value of r0 at lambda2 from its value at lambda1
*/
{
  return (lambda2/lambda1)^(6./5) * r0_at_lambda1;
}

func plsh(slp,nwfs,scale){
  slpx = slp(:numberof(slp)/2);
  slpy = slp(numberof(slp)/2 + 1:);

  s2ipup = (dimsof(*y_geom._ipupil)(2) - dimsof(*y_geom._spupil)(2))/2.;
  posx = *y_wfs(1)._istart + s2ipup ;
  posx = (posx * (*y_wfs(1)._isvalid))(*);
  posx = posx(where(posx!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1;
  posy = *y_wfs(1)._jstart + s2ipup ;
  posy = (transpose(posy * (*y_wfs(1)._isvalid)))(*);
  posy = posy(where(posy!=0)) - dimsof(*y_geom._ipupil)(2)/2 - 1;

  fma;
  plv,slpy(,-:1:2),slpx(,-:1:2),posy(,-:1:2),posx(,-:1:2),scale=scale;

}

func pldm(com){

  yoga_setcomm,g_dm,"pzt",0.,com;
  yoga_shapedm,g_dm,"pzt",0.;
  dm = yoga_getdm(g_dm,"pzt",0.);
  yoga_resetdm,g_dm,"pzt",0.;

  return dm;
}