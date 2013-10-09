/*
 * util_fr.i
 * A collection of routines for general purpose.
 *
 * $Id: util_fr.i,v 1.3 2008/10/29 15:58:13 paumard Exp $
 *
 * Author: Francois Rigaut.
 * Written 2002
 * last revision/addition: 2004Oct15
 *
 * Copyright (c) 2003, Francois RIGAUT (frigaut@gemini.edu, Gemini
 * Observatory, 670 N A'Ohoku Place, HILO HI-96720).
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
 * $Log: util_fr.i,v $
 * Revision 1.3  2008/10/29 15:58:13  paumard
 * utils.i: reform would not work with empty dimlist. Fixed.
 * plot.i, util_fr.i, utils.i: rename functions now standard in Yorick (color_bar, rdfile, reform)
 *
 * Revision 1.2  2007/12/27 15:22:07  frigaut
 * nothing. commit before tagging.
 *
 */

require,"style.i";

/*************************/
/* CONVENIENCE functions */
/*************************/

func ls 
/* DOCUMENT ls: system command ls 
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: pwd, system, $.
 */
{system,"ls";}

func exist(arg)
/* DOCUMENT exist(arg)
 * Returns 0 if element is not set or is a <nuller>, 1 otherwise
 * F.Rigaut 2002/04/03
 * SEE ALSO: is_void, where
 */
{
  if (numberof(arg) == 0) {return 0;}
  else {return 1;}
}



func is_set(arg)
/* DOCUMENT is_set(arg)
 * Returns 0 if element is void or equal to zero, 1 otherwise
 * F.Rigaut 2002/06/03
 * SEE ALSO: is_void, where, exist
 */
{
  if (is_void(arg) | (arg == 0)) {return 0;}
  else {return 1;}
}



func tv(im,square=) 
/* DOCUMENT tv(im,square=) 
 * This routines does a frame advance, display the image
 * and set the limits to have the image full display.
 * Inspired from the IDL tvscl
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: fma, pli, plot
 */
{
  fma; 
  pli,im; 
  limits,"e","e","e","e",square=square;
}



func plot(vect,x,square=,histo=) 
/* DOCUMENT plot(vect,x,square=,histo=)
 * Short cut for a fma + plg
 * Set histo to get plh type plot
 * F.Rigaut 2001/10
 * SEE ALSO: plg, fma, tv, plh
 */
{
  fma;
  if (is_set(histo)) {
    plh,vect,x;
  } else {
    plg,vect,x;
  }
  limits,,,,,square=square;
}



func nprint(var,sep=,format=)
/* DOCUMENT func nprint(var,sep=,format=)
   Neat print of a 2d array.
   example:
   > nprint,optpos*pi/3.14e9,sep=", "
      +0, -5.003e-07,       +0, -9.005e-08,       +0,       +0 
      +0, -4.002e-07,       +0, +9.005e-08,       +0,       +0 
      +0, -3.002e-07,       +0, +9.005e-08,       +0,       +0 
      +0, -2.001e-07,       +0, +9.005e-08,       +0,       +0 
      +0, -1.801e-07,       +0, +9.005e-08,       +0,       +0 
      +0, -1.001e-07,       +0, +9.005e-08,       +0,       +0 
      +0, +1.001e-07,       +0, +9.005e-08,       +0,       +0
   sep= separator string. The default separator is two blanks ("  ").
   format= swrite format
   Restricted to 2D arrays
   SEE ALSO: pm
 */
{
  if (!is_set(sep)) sep = " ";
  if (!is_set(format)) format = "%+8.4g";
  dim = dimsof(var);
  if (dim(1) != 2) {error,"only implemented for 2D arrays";}
  for (i=1;i<=dim(3);i++) {
    for (j=1;j<=dim(2)-1;j++) {
      write,format=format+sep,var(j,i);
    }
    write,format=format,var(0,i);
    write,"";
  }
}



func typeReturn(void)
/* DOCUMENT typeReturn(void)
 * A simple convenient function that does what is name says.
 * SEE ALSO:
 */
{
  rep = rdline(prompt="type return to continue..."); 
  return rep;
}
hitReturn=typeReturn;



/****************************/
/* ADDITIONAL MATH function */
/****************************/

func round(arg) 
/* DOCUMENT round(arg)
 * Returns the rounded version of a floating point argument
 * modified 2007dec06 to fix problem with negative numbers
 * F.Rigaut 2001/10
 * SEE ALSO: ceil, floor
 */
{return long(arg+0.5)-(arg<0);}



func even(arg) 
/* DOCUMENT even(arg)
 * returns 1 is argument is even, zero otherwise. 
 * The argument should be an integer.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: odd
 */
{return ((arg % 2) == 0);}



func odd(arg) 
/* DOCUMENT odd(arg)
 * Returns 1 is argument is odd, zero otherwise. 
 * The argument should be an integer.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: even
 */
{return ((arg % 2) == 1);}



func minmax(arg) 
/* DOCUMENT minmax(arg)
 * Returns a vector containing the min and the max of the argument
 * F.Rigaut 2001/09
 * SEE ALSO:
 */
{return [min(arg),max(arg)];}



local clip;
func __clip(arg,lt,ht)
/* DOCUMENT func clip(arg, mini, maxi);
 * Returns the argument, which has been "clipped" to mini
 * and maxi, i.e. in which all elements lower than "mini"
 * have been replaced by "mini" and all elements greater
 * than "maxi" by "maxi". Array is converted to float.
 * Either "mini" and "maxi" can be ommited, in which case
 * the corresponding mini or maxi is not clipped.
 * Equivalent to the IDL ">" and "<" operators.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO:
 */
{
  if (lt != []) arg = max(arg,lt);
  if (ht != []) arg = min(arg,ht);
  return arg;
}
if (!is_func(clip)) clip = __clip;



local sinc;
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
if (!is_func(sinc)) sinc = __mysync;



/************************************/
/* SYSTEM AND PERFORMANCE functions */
/************************************/

func tic(counterNumber)
/* DOCUMENT tic(counter_number)
 * Marks the beginning of a time lapse
 * ex: tic ; do_something ; tac()
 * will print out the time ellapsed between tic and tac
 * a counter number can optionaly be specified if several
 * counters have to be used in parallel.
 * F.Rigaut 2001/10
 * SEE ALSO: tac
 */
{
  if (counterNumber == []) counterNumber = 1;
  if (counterNumber > 10) error,"tic and tac are limited to 10 time counters !";

  el = array(double,3);
  timer,el;
  _nowtime(counterNumber) = el(3);
}
extern _nowtime;
 _nowtime = array(double,10);



func tac(counterNumber)
/* DOCUMENT tac(counter_number)
 * Marks the end of a time lapse
 * ex: tic ; do_something ; tac()
 * will print out the time ellapsed between tic and tac
 * a counter number can optionaly be specified if several
 * counters have to be used in parallel.
 * F.Rigaut 2001/10
 * SEE ALSO: tic
 */
{
  if (counterNumber == []) counterNumber = 1;

  el = array(double,3);
  timer,el;
  elapsed = el(3)-_nowtime(counterNumber);

  return elapsed;
}



func spawn_fr(command) 
/* DOCUMENT spawn_fr(command)
 * This function tries to group in one call the :
 * - call to system
 * - read the file created by the system call
 * - returns it
 * Inspired from the IDL function of the same name
 * F.Rigaut 2002/04/04
 * SEE ALSO: system, popen, exec in Eric/system.i
 * 
 * DEPRECATED: Use spawn/sys instead.
 */
{
  f   = popen(command,0);
  ans = rdline(f);
  l   = ans;
  while (l)
    {
      l = rdline(f);
      if (l) {ans = grow(ans,l);}
    }
  return ans;
}  



/*******************/
/* ARRAY functions */
/*******************/

func wheremin(ar) { return  where(ar == min(ar)); }
func wheremax(ar) { return  where(ar == max(ar)); }
/* DOCUMENT func wheremin(ar)
        and func wheremax(ar)
   Short hand for where(array == min(array) or max(array)
   SEE ALSO: where, where2, min, max
 */



func indices(dim)
/* DOCUMENT indices(dim)
 * Return a dimxdimx2 array. First plane is the X indices of the pixels
 * in the dimxdim array. Second plane contains the Y indices.
 * Inspired by the Python scipy routine of the same name.
 * New (June 12 2002): dim can either be :
 *   - a single number N (e.g. 128) in which case the returned array are
 *     square (NxN)
 *   - a Yorick array size, e.g. [#dimension,N1,N2], in which case
 *     the returned array are N1xN2
 *   - a vector [N1,N2], same result as previous case
 * F.Rigaut 2002/04/03
 * SEE ALSO: span
 */
{
  if (numberof(dim) == 1)
  {
    x  = span(1,dim,dim)(,-:1:dim);
    y  = transpose(x);
    return [x,y];
  } else
  {
    if (numberof(dim) == 3) {dim = dim(2:3);}
    x  = span(1,dim(1),dim(1))(,-:1:dim(2));
    y  = span(1,dim(2),dim(2))(,-:1:dim(1));
    y  = transpose(y);
    return [x,y];
  }
    
}



local dist;
func __dist(dim,xc=,yc=)
/* DOCUMENT func dist(size,xc=,yc=)
 * Returns an array which elements contains the distance to (xc,yc). xc
 * and yc can be omitted, in which case they are defaulted to size/2+1.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: indices, radial_distance
 */
{
  dim	= long(dim);
  if (xc == []) xc = int(dim/2)+1;
  if (yc == []) yc = int(dim/2)+1;
  x	= float(span(1,dim,dim)(,-:1:dim));
  y	= transpose(x);
  d	= float(sqrt((x-xc)^2.+(y-yc)^2.));
  d	= clip(d,1e-5,);
  return d;
}
if (!is_func(dist)) {dist = __dist;}



local eclat;
func __eclat(image)
/* DOCUMENT func eclat(image)
 * Equivalent, but slightly faster (?) than roll. Transpose the four main
 * quadrants of a 2D array. Mostly used for FFT applications.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: roll.
 */
{
  d	= dimsof(image);
  dx	= d(2);
  dy	= d(3);
  x1=1; x2=dx/2 ; x3=x2+1 ; x4=dx;
  y1=1; y2=dy/2 ; y3=y2+1 ; y4=dy;
  out	= image*0.;
  out(x1:x2,y1:y2) = image(x3:x4,y3:y4);
  out(x3:x4,y1:y2) = image(x1:x2,y3:y4);
  out(x1:x2,y3:y4) = image(x3:x4,y1:y2);
  out(x3:x4,y3:y4) = image(x1:x2,y1:y2);
  return out;
}
if (!is_func(eclat)) {eclat = __eclat;}



func makegaussian(size,fwhm,xc=,yc=,norm=) 
/* DOCUMENT makegaussian(size,fwhm,xc=,yc=)
 * Returns a centered gaussian of specified size and fwhm.
 * F.Rigaut 2001/09
 * norm returns normalized 2d gaussian
 * SEE ALSO:
 */
{
  tmp = exp(-(dist(size,xc=xc,yc=yc)/(fwhm/1.66))^2.);
  if (is_set(norm)) tmp = tmp/fwhm^2./1.140075;
  return tmp;
}



func bin2(image)
/* DOCUMENT bin2(image)
 * Returns the input image, binned by a factor of 2.
 * That is, a 512x512 image is transformed in a 256x256 image.
 * one output pixel is the average of the 4 corresponding input ones,
 * so that it conserves the total intensity.
 * SEE ALSO: undersample
 */
{
  d = dimsof(image);
  if (d(1) != 2) {
    error,"Bin only accepts images";
  }
  if (((d(2) % 2) != 0) || ((d(3) % 2) != 0)) {
    error,"Bin only accepts dimensions with even # of pixels";
  }
  sim= image+roll(image,[-1,-1])+roll(image,[-1,0])+roll(image,[0,-1]);
  return sim(::2,::2);
}



/*****************************/
/* STRING AND FILE functions */
/*****************************/

func fileExist(filename)
/* DOCUMENT fileExist(filename)
 *  Returns "1" if the file(s) exist(s), "0" if it does not.
 *  filename can be an array, in which case the results is an array
 *  of 0s and 1s.
 *  F.Rigaut 2002/01
 *  SEE ALSO:
 */
{
  exist = [];
  for (i=1;i<=numberof(filename);i++) {
    grow,exist,(open(filename(i),"r",1) ? 1:0);
  }
  if (numberof(filename) == 1) {return exist(1);}
  return exist;
}



func findfiles(files)
/* DOCUMENT findfiles(filter)
 * This routines returns a list of files which satisfy the filter
 * argument. The list is a string vector. If no files were found,
 * the results is the empty string.
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: spawn
 */
{
  // parse input parameter into path and filename:
  tmp = strtok(files,"/",20);
  tmp = tmp(where(tmp));
  
  filereg = tmp(0);

  if (numberof(tmp)>1) {
    path = sum(tmp(1:-1)+"/");
    if (strpart(files,1:1)=="/") path = "/"+path;
  } else path=".";

  l = lsdir(path);

  // process result list:
  if (noneof(l)) return;

  w = where(strglob(filereg,l));
  if (numberof(w)==0) return;
  
  res = l(w);
  if (path==".") return res;

  return (path+res);
}



func __rdfile(file)
/* DOCUMENT func rdfile(file)
   Open, read, close and return the whole content of ascii file "file".
   AUTHOR : F.Rigaut, Oct 2004.
   SEE ALSO: load_text, dump_text
 */
{
  f = open(file,"r");
  fcontent = [];
  while (line=rdline(f)) grow,fcontent,line;
  return fcontent;
}
if (!rdfile) rdfile=__rdfile;




func parsedate_fr(strdate,format,dec=)
/* DOCUMENT parsedate(strdate,format,dec=)
 * Returns the date in integers as an array [[year],[month],[day]],
 * or in decimal if "dec" is set. The format has to be specified.
 * Does not yet accept dates in the form "2003jun05".
 * strdate = string array containing the date e.g. "2003/10/25"
 * format  = format in the form "yyyy/mm/dd"
 * Examples:
 * > parsedate(["2003/05/21","2003/02/15"],"yyyy/mm/dd",dec=1)
 * [2003.38,2003.12]
 * > parsedate(["2003/05/21","2003/02/15"],"yyyy/mm/dd")
 * [[2003,2003],[5,2],[21,15]]
 * SEE ALSO: parsetime
 ***** NEED TO UPGRADE THIS FUNCTION FOR THE 1.6 STR FUNCTIONS ****
 */
{
  exit,"not upgraded for use with the yorick-1.6.01 str functions";

  a = b = c = array(long,numberof(strdate));
  strdate = strtrim(strdate);

  format = strtolower(format);
  wy = strfind(format,"y"); ys = strjoin(array("y",numberof(wy)));
  wm = strfind(format,"m"); ms = strjoin(array("m",numberof(wm)));
  wd = strfind(format,"d"); ds = strjoin(array("d",numberof(wd)));
  format = strreplace(format,ys,"%"+swrite(numberof(wy),format="%d")+"d");
  format = strreplace(format,ms,"%"+swrite(numberof(wm),format="%d")+"d");
  format = strreplace(format,ds,"%"+swrite(numberof(wd),format="%d")+"d");

  sread,strdate,format=format,a,b,c;

  if (wy(1) < wm(1) & wm(1) < wd(1)) ymd= [a,b,c];
  if (wd(1) < wm(1) & wm(1) < wy(1)) ymd= [c,b,a];
  if (wm(1) < wd(1) & wd(1) < wy(1)) ymd= [c,a,b];
  if (wy(1) < wd(1) & wd(1) < wm(1)) ymd= [a,c,b];
  if (wd(1) < wy(1) & wy(1) < wm(1)) ymd= [b,c,a];
  if (wm(1) < wy(1) & wy(1) < wd(1)) ymd= [b,a,c];

  if (!dec) {return ymd;}

  mlength = [31.,28,31,30,31,30,31,31,30,31,30,31];
  // rought cut at bisextile years:
  bisext = (ymd(..,1)/4. == ymd(..,1)/4);
  // none of the date are bisextile:
  if (noneof(bisext)) {
    return ymd(..,1)+(mlength(cum)(ymd(..,2))+ymd(..,3)-1.)/mlength(sum);
  }

  // all of the date are bisextile:
  if (allof(bisext)) {
    mlength(2) = 29;
    return ymd(..,1)+(mlength(cum)(ymd(..,2))+ymd(..,3)-1.)/mlength(sum);
  }

  // mix of bisextile and not:
  tnotbi = ymd(..,1)+(mlength(cum)(ymd(..,2))+ymd(..,3)-1.)/mlength(sum);
  mlength(2) = 29;
  tbi    = ymd(..,1)+(mlength(cum)(ymd(..,2))+ymd(..,3)-1.)/mlength(sum);
  return tnotbi*(1-bisext) + tbi*bisext;
}



func parsetime_fr(strtime,dec=)
/* DOCUMENT parsetime(strtime,dec=)
 * Parse the input string array "strtime" and returns the time as a [hh,mm,ss]
 * vector (or array if strtime is an array of string) or returns the time
 * as a single decimal number (or vector if strtime is an array) if the
 * "dec" keyword is set.
 * Examples:
 * > parsetime(["22:30:25.792","22:20:33.852"],dec=1)
 * [22.5072,22.3427]
 * > parsetime(["22:30:25.792","22:20:33.852"])
 * [[22,22],[30,20],[25.792,33.852]]
 ***** NEED TO UPGRADE THIS FUNCTION FOR THE 1.6 STR FUNCTIONS ****
 * SEE ALSO: parsedate
 */
{
  exit,"not upgraded for use with the yorick-1.6.01 str functions";
  n = numberof(strtime);
  h = m = array(long,n);
  s = array(float,n);
  strtime = strtrim(strtime);
  if (allof(strmatch(strtime,":"))) {delim=":";} \
  else if (allof(strmatch(strtime," "))) {delim=" ";} \
  else if (allof(strmatch(strtime,";"))) {delim=";";} \
  else {error,"Can't figure out the delimiter";}
  f = "%2d"+delim+"%2d"+delim+"%f";
  sread,strtime,format=f,h,m,s;
  if (dec) return h+m/60.+s/3600.;
  return [h,m,s];
}


func secToHMS(time)
/* DOCUMENT secToHMS(time)
 * Convert from time (float in sec) to string in the form hh:mm:ss.s
 * AUTHOR : F.Rigaut, June 13 2002.
 * SEE ALSO: 
 */
{
  lt    = long(time);
  hh    = lt/3600;
  lt    = lt-hh*3600;
  mm    = lt/60;
  sec   = float(time)-hh*3600-mm*60;
  ts    = swrite(hh,mm,sec,format="%02i:%02i:%04.1f");
  return ts;
} 


/*********************/
/* GRAPHIC functions */
/*********************/

func colorbar(cmin, cmax,adjust=,levs=)
/* DOCUMENT colorbar
            colorbar, cmin, cmax
     draw a color bar to the right of the plot.  If CMIN and CMAX
     are specified, label the top and bottom of the bar with those
     numbers.
     adjust: x adjust. typically +/- 0.01 
     levs: number of ticks in color bar (plus one)
     upgraded 2007june02
   SEE ALSO: color_bar
 */
{
  if (adjust==[]) adjust=0.;
  cursys = plsys();

  get_style,landscape, systems, legends, clegends;

  left   = systems(cursys).viewport(2)+0.03+adjust;
  right  = systems(cursys).viewport(2)+0.05+adjust;
  bottom = systems(cursys).viewport(3);
  top    = systems(cursys).viewport(4);
  middle = (left+right)/2.;
  
  plsys, 0;
  pli, span(0,1,200)(-,), left, bottom, right, top, legend="";
  plg, [bottom,top,top,bottom],[right,right,left,left], closed=1,
    marks=0,color="fg",width=1,type=1,legend="";
  if (levs) {
    for (i=1;i<=(levs-1);i++) {
      y=bottom+(top-bottom)/levs*i;
      plg, [y,y],[left,right], marks=0,color="fg",width=1,type=1,legend="";
    }
  }
  plsys, cursys; 
  if (!is_void(cmin)) {
    plt, pr1(cmin), middle, bottom-0.005, justify="CT";
    plt, pr1(cmax), middle, top+0.005, justify="CB";
  }
}

// general equivalent (remove if you don't like them)
man=help;
