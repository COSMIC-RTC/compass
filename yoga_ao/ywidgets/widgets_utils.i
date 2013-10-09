/*
yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";
*/

struct display_struct
{
  long    _lut;
  long    _ncolors;
  long    _itt;
  long    _log_itt_dex;
  long    _gui_realized;
  long    _invertlut;
  long    _defaultdpi;
  float   _zoom;
  float   _cmin;
  float   _cmax;
  string  _cmd;
  pointer _wins;
  pointer _winits;
  pointer _rlut;
  pointer _glut;
  pointer _blut;
  pointer _imdnum;
  pointer _xid;
  pointer _xytitles_adjust;
};

func widget_set_lut(_lut,_win,&widget_lut,widget_ncolors,widget_itt,widget_log_itt_dex,rlut,glut,blut,invertlut=)
// change the LookUp Table and Intensity Transfer Table
{
  require,"idl-colors.i";

  if (is_void(invertlut)) invertlut = 0;
  local r,g,b;

  if (_lut!=[]) widget_lut = _lut;
  
  if (_lut!=[]) {  // then read and set new lut
    if (_lut==0) palette,"earth.gp";
    else loadct,_lut;
    palette,query=1,rlut,glut,blut;  // store
  }

  // invert?
  if (invertlut) {
    r=rlut(::-1); g=glut(::-1); b=blut(::-1);
  } else {
    r=rlut; g=glut; b=blut;
  }

  // itt:
  if (widget_itt<=1) { // linear
    ind = span(0.,1.,widget_ncolors);
  } else if (widget_itt==2) { // sqrt
    ind = sqrt(span(0.,1.,widget_ncolors));
  } else if (widget_itt==3) { // square
    ind = (span(0.,1.,widget_ncolors))^2.;
  } else if (widget_itt==4) { // log
    ind = log10(span(10.^(-widget_log_itt_dex),1.,widget_ncolors)); // 8 dex
    ind -= min(ind);
    ind /= max(ind);
  } else if (widget_itt>=5) { // histeq
    ind = span(0.,1.,widget_ncolors);
  }
  ind = long(round(ind*(widget_ncolors-1)+1));
  r = r(ind); g = g(ind); b = b(ind);

  // and finally, load the palette:
  window,_win;
  palette,r,g,b;
}

func widget_disp(_win,_im,&widget_imdnum,widget_itt,widget_zoom,widget_defaultdpi,disp_name,disp_xtitle,disp_ytitle)
// pli display, main window
{
  _imd = _im*0.;
  
  for (i=1;i<=2;i++) {
    if (widget_itt==5) {
      if (nallof(widget_imdnum==[2,widget_itt]))     \
        _imd = widget_histeq_scale(_im);
    } else _imd = bytscl(_im,cmin=widget_cmin,cmax=widget_cmax);

    window,_win;
    fma;
    if (widget_zoom != 1) {
      npix = dimsof(_im)(2);
      newpix = long(npix/widget_zoom);
      mygap = long((npix-newpix)/2.+1);
      pli,_imd(mygap:mygap+newpix-1,mygap:mygap+newpix-1);
    } else pli,_imd;

    disp_xtitle = "pixels";
    disp_ytitle = "pixels";
    
    disp_pltitle,disp_name,defaultdpi=widget_defaultdpi;
    disp_xytitles,disp_xtitle,disp_ytitle,defaultdpi=widget_defaultdpi;
    colorbar,adjust=-0.024,levs=10;
  }
  
  widget_imdnum = [2,widget_itt];
}

func widget_disp_cpc(e,_im,&widget_cmin,&widget_cmax,all=)
{
  if (all) eq_nocopy,subim,_im;
  else subim=widget_get_subim(x1,x2,y1,y2,_im);
  
  if (subim==[]) subim = _im; // init, not defined
  
  if (e==0) {
    cmin = min(subim);
    cmax = max(subim);
  } else {
    cmin = min(widget_cpc(subim,0.1,0.999));
    cmax = max(widget_cpc(subim,0.1,0.999));
  }
  
  widget_cmin = cmin;
  widget_cmax = cmax;
  /*
  pyk,swrite(format=cmd_pray+"widget_set_cmincmax(%f,%f,%f)",float(cmin),float(cmax),float(cmax-cmin)/100.);
  */
}


func widget_cpc(im,fmin,fmax)
/* DOCUMENT func cpc(im,fmin,fmax)
   return clipped image at
   from cmin=fraction fmin of pixel intensity
   to   cmax=fraction fmax of pixel intensity
   0 <= fmin < fmax <= 1
   example: pli,cpc(im)
   SEE ALSO:
 */
{
  s = sedgesort(im(*));
  n = numberof(im);
  if (fmin==[]) fmin = 0.10;
  if (fmax==[]) fmax = 0.995;
  if ((fmin<0)||(fmax<0)) error,"fmin and fmax should be > 0";
  if ((fmin>1)||(fmax>1)) error,"fmin and fmax should be < 1";
  if (fmin>fmax) error,"fmin should be < fmax";
  x1=s(long(round(n*fmin)));
  x2=s(long(round(n*fmax)));
  return clip(im,x1,x2);
}

func widget_set_cmin(pycmin,&widget_cmin,widget_cmax)
{
  if ((widget_cmax-widget_cmin)!=0) if ((abs(pycmin-widget_cmin)/(_cmax-widget_cmin))<1e-3) return;
  if (pycmin>widget_cmax) {
    write,"cmin > cmax, ignoring";
    return;
  }
  widget_cmin = pycmin;
}

func widget_set_cmax(pycmax,widget_cmin,&widget_cmax)
{
  if ((widget_cmax-widget_cmin)!=0) if ((abs(pycmax-widget_cmax)/(_cmax-widget_cmin))<1e-3) return;
  if (pycmax<widget_cmin) {
    write,"cmax < cmin, ignoring";
    return;
  }
  widget_cmax = pycmax;
}

func widget_get_subim(_win,_im,&x1,&x2,&y1,&y2)
{
  curw = current_window();
  window,_wins;
  lim=limits();
  dims = dimsof(_im);
  x1=long(floor(clip(lim(1),1,dims(2))));
  x2=long(ceil(clip(lim(2),1,dims(2))));
  y1=long(floor(clip(lim(3),1,dims(3))));
  y2=long(ceil(clip(lim(4),1,dims(3))));
  if ( (x1==x2) || (y1==y2) ) {
    write,"WARNING: (get_subim) Nothing to show";
    return;
  }
  window,curw;
  return _im(x1:x2,y1:y2)
}

func widget_unzoom(_win,&widget_zoom)
{
  window,_wins;
  unzoom;
  limits;

  widget_zoom = 1;
}

func widget_setcuts
{
  widget_disp_cpc,all=1;
}

func widget_histeq_scale(z, top=, cmin=, cmax=)
/* DOCUMENT histeq_scale(z, top=top_value, cmin=cmin, cmax=cmax)
     returns a byte-scaled version of the array Z having the property
     that each byte occurs with equal frequency (Z is histogram
     equalized).  The result bytes range from 0 to TOP_VALUE, which
     defaults to one less than the size of the current palette (or
     255 if no pli, plf, or palette command has yet been issued).

     If non-nil CMIN and/or CMAX is supplied, values of Z beyond these
     cutoffs are not included in the frequency counts.

     Identical to histeq_scale except it uses sedgesort instead of sort.
     faster for arrays for which many elements are repeated (e.g.
     CCD arrays where pixels values are integers.
   SEE ALSO: bytscl, plf, pli
 */
{
  if (is_void(top)) top= bytscl([0.,1.])(2);  /* palette size - 1 */
  top= long(top);
  if (top<0 | top>255) error, "top value out of range 0-255";
  y= z(*);
  if (!is_void(cmin)) y= y(where(y>=cmin));
  if (!is_void(cmax)) y= y(where(y<=cmax));
  y= sedgesort(y);
  x= span(0.,1., numberof(y));
  xp= span(0.,1., top+2);
  bins= interp(y, x, xp);
  list= where(bins(dif)<=0.0);
  if (numberof(list)) {
    /* some value (or values) of z are repeated many times --
       try to handle this by adding a small slope to the sorted y */
    dy= y(0)-y(1);
    if (!dy) dy= 1.0;
    for (eps=1.e-10 ; eps<1000.1 ; eps*=10.) {
      bins= interp(y+eps*dy*x, x, xp);
      list= where(bins(dif)<=0.0);
      if (!numberof(list)) break;
    }
    if (eps>1000.) error, "impossible error??";
  }
  return char(max(min(digitize(z,bins)-2,top),0));
}

func widget_pltitle(title,&title_height,title_font,defaultdpi=)
{
  if (defaultdpi == []) defaultdpi = 70;
  plth_save = title_height;
  title_height = long(title_height*defaultdpi/83.);
  
  port= viewport();
  plth=title_height;

  plt, escapechar(title), port(zcen:1:2)(1), port(4)+0.005,
    font=title_font, justify="CB", height=plth;

  title_height = plth_save;
}


func widget_xytitles(xtitle,ytitle,adjust,&title_height,defaultdpi=)
{
  if (defaultdpi == []) defaultdpi = 70;
  plth_save = title_height;
  title_height = long(title_height*defaultdpi/83.);
  
  curw=current_window();
  if (adjust==[]) {
    adjust = [0.012,0.019];
  }
  xytitles,escapechar(xtitle),escapechar(ytitle),adjust;

  title_height = plth_save;
}

func escapechar(s)
{
  if (s==[]) return;
  s=streplace(s,strfind("_",s,n=20),"!_");
  s=streplace(s,strfind("^",s,n=20),"!^");
  return s;
}

func yoga_pltitle(title,adjust)
/* DOCUMENT pltitle, title
     Plot TITLE centered above the coordinate system for any of the
     standard Gist styles.  You may want to customize this for other
     plot styles.
   SEE ALSO: plt, xytitles
 */
{
  port= viewport();
  plt, title, port(zcen:1:2)(1)+adjust(1), port(4)+0.02+adjust(2),
    font=pltitle_font, justify="CB", height=pltitle_height;
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

func display_slopes(slopes,nwfs,type)
{

  slX = slopes(1:y_wfs(nwfs)._nvalid);
  slY = slopes(y_wfs(nwfs)._nvalid+1:);
  uu = where(*y_wfs(nwfs)._isvalid);
  
  vx=vy=array(0.,[2,y_wfs(nwfs).nxsub,y_wfs(nwfs).nxsub]);
  vx(uu)=slX;
  vy(uu)=slY;
  
  xy = indices(y_wfs(nwfs).nxsub);
  xs = xy(,,1)-0.5;
  ys = xy(,,2)-0.5;
  
  //window,(*wfs_disp._wins)(1);fma;
  plvf,vy(uu),vx(uu),ys(uu),xs(uu),color="fg";
  
  xs = xs(uu);
  ys = ys(uu);
  
  for (i=1;i<=numberof(xs);i++) {
    plg,[ys(i)-0.5,ys(i)+0.5],[xs(i)-0.5,xs(i)-0.5],marks=0,color="red";
    plg,[ys(i)-0.5,ys(i)+0.5],[xs(i)+0.5,xs(i)+0.5],marks=0,color="red";
    plg,[ys(i)-0.5,ys(i)-0.5],[xs(i)-0.5,xs(i)+0.5],marks=0,color="red";
    plg,[ys(i)+0.5,ys(i)+0.5],[xs(i)-0.5,xs(i)+0.5],marks=0,color="red";
  }
  avgX = vx(uu)(avg); avgY = vy(uu)(avg);
  rmsX = vx(uu)(rms); rmsY = vy(uu)(rms);

  myxtitle = swrite(format="%s on WFS # %d",type,nwfs);
  port= viewport();
  plt, myxtitle, port(zcen:1:2)(1), port(4)+0.005,
    font="helvetica", justify="CB", height=long(pltitle_height*wfs_disp._defaultdpi/200.);

  adjust = 0.019;
  myxtitle = swrite(format="avg: x=%.2f\" y=%.2f\" rms: x=%.2f\" y=%.2f\"",avgX,avgY,rmsX,rmsY);
  plt,myxtitle , port(zcen:1:2)(1), port(3)-0.050+adjust,
      font="helvetica", justify="CT", height=long(pltitle_height*wfs_disp._defaultdpi/200.);
}

