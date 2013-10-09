require,"yoga.i";
require,"util_fr.i";
  
func check_yoga_fft(size)
{
  if (size == []) size = 256;
  a = yoga_random("scomplex",[2,size,size]);
  res = a();
  b=res(1,,)+1i*res(2,,);
  
  yoga_fft,a;
  res = a();
  c = res(1,,)+1i*res(2,,);
  res2 = fft(b,1);
  
  "Error Re: ";
  minmax(c.re-float(res2.re))/max([max(abs(res2.re)),max(abs(res2.im))]);
  "Error Im: ";
  minmax(c.im-float(res2.im))/max([max(abs(res2.re)),max(abs(res2.im))]);

  "";
  "checkFFT";
  "";
  checkFFT,size;

  "";
  "checkFFTMulti";
  "";
  checkFFTMulti,1000,64,64;
  
  "";
  "checkConvFFT";
  "";
  checkConvFFT,size;
}

/*
       _               _    _____ _____ _____ 
   ___| |__   ___  ___| | _|  ___|  ___|_   _|
  / __| '_ \ / _ \/ __| |/ / |_  | |_    | |  
 | (__| | | |  __/ (__|   <|  _| |  _|   | |  
  \___|_| |_|\___|\___|_|\_\_|   |_|     |_|  
                                             
*/
func checkFFT(sizex, sizey,ctype=) {

  if (is_void(ctype)) ctype = "all";
  if(is_void(sizex)) sizex=512;
  if(is_void(sizey)) sizey=sizex;
  
  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    a = yoga_random("scomplex",[2,sizex,sizey]);
    grow, timeProfile, tac(timeTic);
    res=a();ref=res(1,,)+1i*res(2,,);
    grow, timeProfile, tac(timeTic);
    yoga_fft,a;
    //grow, timeProfile, tac(timeTic);
    //yoga_fft,a;
    grow, timeProfile, tac(timeTic);
    res2 = a(); res=res2(1,,)+1i*res2(2,,);
    grow, timeProfile, tac(timeTic); 
    a = [];
    grow, timeProfile, tac(timeTic);
    
    "";
    "FFT gpu time profile: ";
    "Alloc     h2dm     comp    d2hm    free" ;
    timeProfile; 
    "";
    "FFT gpu time individual: ";
    "h2dm        comp        d2hm       free";
    timeProfile(dif); 
    "";
  }
  if ((ctype == "all") || (ctype == "cpu")) {
    //myTime=tic(); for(i=0; i<nbsamp; i++) res=fft(imd, 1); tac(myTime)/nbsamp;
    timeTic = tic();
    res2=fft(ref, 1);
    timecpu=tac(timeTic);
    if (ctype != "all") timeProfile = tac(timeTic);
    write,format="FFT cpu time : %f\n",tac(timeTic);
    if (ctype == "all") {
      "Error Re: ";
      minmax(res.re-float(res2.re))/max([max(abs(res2.re)),max(abs(res2.im))]);
      "Error Im: ";
      minmax(res.im-float(res2.im))/max([max(abs(res2.re)),max(abs(res2.im))]);
      tmp = timeProfile(dif);
      "";
      write,format="Acceleration factor w memcopy : %.1f\n",timecpu/tmp(1:3)(sum);
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(2);
    }
  }
 
  return timeProfile;
}

func bench_fft(nmin=,nmax=)
{
  device=setDeviceId(0);
  
  if (is_void(nmin)) nmin = 16;
  if (is_void(nmax)) nmax = 4096;

  memcpy_in = compute_fft = memcpy_out = cpures = [];
  
  for (cc=nmin;cc<=nmax;cc++) {
    gpures = checkFFT(cc,ctype="gpu");
    grow,memcpy_in,gpures(dif)(1);
    grow,compute_fft,gpures(dif)(2);
    grow,memcpy_out,gpures(dif)(3);
    grow,cpures,checkFFT(cc,ctype="cpu");
  }

  mysize = indgen(nmax - (nmin-1)) + (nmin-1);
  window,dpi=150;
  //plg,,mysize,marks=0,width = 4;
  plg,cpures/compute_fft,mysize,marks=0,width = 4;
  plg,cpures/(compute_fft+memcpy_in+memcpy_out),mysize,marks=0,width = 4,color="red";
  plg,cpures/(compute_fft+memcpy_out),mysize,marks=0,width = 4,color="green";
  xytitles,"width in pixels","Accel. factor";
  pltitle,"FFT on "+device;
}

/*
       _               _    _____ _____ _____ __  __       _ _   _ 
   ___| |__   ___  ___| | _|  ___|  ___|_   _|  \/  |_   _| | |_(_)
  / __| '_ \ / _ \/ __| |/ / |_  | |_    | | | |\/| | | | | | __| |
 | (__| | | |  __/ (__|   <|  _| |  _|   | | | |  | | |_| | | |_| |
  \___|_| |_|\___|\___|_|\_\_|   |_|     |_| |_|  |_|\__,_|_|\__|_|
                                                                  
*/
func checkFFTMulti(nbsamp, sizex,sizey,ctype=) {
   
  if (is_void(ctype)) ctype = "all";
  
  if(is_void(nbsamp)) nbsamp=1000;

  if(is_void(sizex)) sizex=int(64);
  if(is_void(sizey)) sizey=sizex;

  //calcul en single precision
  a = yoga_random("scomplex",[3,sizex,sizey,nbsamp]);
  imd = a();
  imz=imd(1,,,)+1i*imd(2,,,); //float

  //calcul en double precision
  //res=imd=imz=array(complex,dims_data);
  //for (i=0;i<nbsamp;i+=2) imz_c(,,i) =random(size,size)+1i*random(size,size);
  //typeFFT="Z2Z";

  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    yoga_fft,a;
    grow, timeProfile, tac(timeTic);
    yoga_random,a; 
    grow, timeProfile, tac(timeTic);
    yoga_fft,a;
    grow, timeProfile, tac(timeTic);
    res = a();
    grow, timeProfile, tac(timeTic);
    a=[];
    grow, timeProfile, tac(timeTic);
    
    res=res(1,,,)+1i*res(2,,,);

    "";
    "FFT gpu time profile: ";
    "warm     rand     comp    d2hm    free"  
    timeProfile; 
    "";
    "FFT gpu time individual: ";
    "rand        comp        d2hm       free"  
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    res2 = imz * 0.;
  
    timeTic = tic();
    for (i=1;i<=nbsamp;i++)
      res2(,,i) =fft(imz(,,i), -1);
    timecpu=tac(timeTic);
    if (ctype != "all") timeProfile = tac(timeTic);
    write,format="FFTMulti cpu time : %f\n",tac(timeTic);
  }

  if (ctype == "all") {
    "Error Re: ";
    minmax((res-res2).re)/max([max(abs(res2.re)),max(abs(res2.im))]);
    "Error Im: ";
    minmax((res-res2).im)/max([max(abs(res2.re)),max(abs(res2.im))]);
    tmp = timeProfile(dif);
    "";
    write,format="Acceleration factor w memcopy : %.1f\n",timecpu/tmp(1:3)(sum);
    write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(2);
  }

  return timeProfile;
}

func bench_fftmulti(nmin=,nmax=)
{
  device=setDeviceId(0);

  if (is_void(nmin)) nmin = 10;
  if (is_void(nmax)) nmax = 1000;


  sizex = [8,16,32,64,128];
  for (dd=1;dd<=numberof(sizex);dd++) {
    memcpy_in = compute_fft = memcpy_out = cpures = [];
    for (cc=nmin;cc<=nmax;cc++) {
      gpures = checkFFTMulti(cc,int(sizex(dd)),ctype="gpu");
      grow,memcpy_in,gpures(dif)(1);
      grow,compute_fft,gpures(dif)(2);
      grow,memcpy_out,gpures(dif)(3);
      grow,cpures,checkFFTMulti(cc,int(sizex(dd)),ctype="cpu");
    }
    mysize = indgen(nmax - (nmin-1)) + (nmin-1);
    window,dd,dpi=150;
    //plg,,mysize,marks=0,width = 4;
    plg,cpures/compute_fft,mysize,marks=0,width = 4;
    plg,cpures/(compute_fft+memcpy_in+memcpy_out),mysize,marks=0,width = 4,color="red";
    plg,cpures/(compute_fft+memcpy_out),mysize,marks=0,width = 4,color="green";
    xytitles,"# of images","Accel. factor";
    pltitle,swrite(format="FFT Multi on %dx%d images on "+device,sizex(dd),sizex(dd));
  }
}


func bench_fftmulti2(nmin=,nmax=)
{
  device=setDeviceId(0);

  if (is_void(nmin)) nmin = 10;
  if (is_void(nmax)) nmax = 1000;


  sizex = [8,16,32,64,128];
  for (dd=1;dd<=numberof(sizex);dd++) {
    memcpy_in = compute_fft = memcpy_out = cpures = [];
    for (cc=nmin;cc<=nmax;cc++) {
      gpures = checkFFTMulti(cc,int(sizex(dd)),ctype="gpu");
      grow,memcpy_in,gpures(dif)(1);
      grow,compute_fft,gpures(dif)(2);
      grow,memcpy_out,gpures(dif)(3);
    }
    mysize = indgen(nmax - (nmin-1)) + (nmin-1);
    window,dd;fma;limits;
    //plg,,mysize,marks=0,width = 4;
    plg,compute_fft,mysize,marks=0,width = 4;
    xytitles,"# of images","Exec Time (s)";
    pltitle,swrite(format="FFT Multi on %dx%d images on "+device,sizex(dd),sizex(dd));
  }
}

/*
      _               _     ____                 _____ _____ _____ 
  ___| |__   ___  ___| | __/ ___|___  _ ____   _|  ___|  ___|_   _|
 / __| '_ \ / _ \/ __| |/ / |   / _ \| '_ \ \ / / |_  | |_    | |  
| (__| | | |  __/ (__|   <| |__| (_) | | | \ V /|  _| |  _|   | |  
 \___|_| |_|\___|\___|_|\_\\____\___/|_| |_|\_/ |_|   |_|     |_|  
                                                                   
*/  
func checkConvFFT(size,sizeker,width, method,ctype=)
{
  extern ref, hres;
  if (is_void(ctype)) ctype = "all";
  if(is_void(size)) size=1024;
  if(is_void(sizeker)) sizeker=512;
  if(is_void(width)) width=20;
  
  hdata = float(makegaussian(size,10));
  hkern = float(makegaussian(size,width));
  hkern2 = float(makegaussian(sizeker,width));
  //sizefft = snapTransformSize(size + sizeker -1);


  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    a=yoga_obj(hdata);
    grow, timeProfile, tac(timeTic);
    k=yoga_obj(hkern2);
    grow, timeProfile, tac(timeTic);
    pad_a=yoga_fftconv_init(a,k,"real");
    grow, timeProfile, tac(timeTic);
    pad_spec=yoga_fftconv_init(a,k,"complex");
    grow, timeProfile, tac(timeTic);
    c=yoga_obj(hdata*0.0f);
    grow, timeProfile, tac(timeTic);
    yoga_fftconv,c,a,k,pad_a,pad_spec,int(sizeker/2),int(sizeker/2);
    grow, timeProfile, tac(timeTic);
    //write,format="gpu time for convFFT (only) %f\n",tac(myTime);
    hres=c();
    grow, timeProfile, tac(timeTic);
    a=[];k=[];pad_a=[];pad_spec=[];c=[];
    grow, timeProfile, tac(timeTic);
    "convFFT gpu time profile: "; timeProfile; 
    "convFFT gpu time individual: "; timeProfile(dif); 
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic=tic();
    ref = float(fft(fft(hdata,1)*fft(hkern,1),-1)/(size)^2);
    write,format="convFFT cpu time : %f\n",tac(timeTic);  
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      "minmax(ref - eclat(hres))";
      minmax(ref - eclat(hres))/max(ref);
    }
  }

  return timeProfile;
}

func checkConvFFT3d(size,sizeker,nim,width, method,ctype=)
{
  extern ref, hres;
  if (is_void(ctype)) ctype = "all";
  if(is_void(size)) size=128;
  if(is_void(sizeker)) sizeker=64;
  if(is_void(width)) width=20;
  if(is_void(nim)) nim=20;
  
  hdata = float(makegaussian(size,10))(,,-::nim-1);
  hkern = float(makegaussian(size,width))(,,-::nim-1);
  hkern2 = float(makegaussian(sizeker,width))(,,-::nim-1);
  //sizefft = snapTransformSize(size + sizeker -1);


  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    a=yoga_obj(hdata);
    grow, timeProfile, tac(timeTic);
    k=yoga_obj(hkern2);
    grow, timeProfile, tac(timeTic);
    pad_a=yoga_fftconv_init(a,k,"real");
    grow, timeProfile, tac(timeTic);
    pad_spec=yoga_fftconv_init(a,k,"complex");
    grow, timeProfile, tac(timeTic);
    c=yoga_obj(hdata*0.0f);
    grow, timeProfile, tac(timeTic);
    yoga_fftconv,c,a,k,pad_a,pad_spec,int(sizeker/2),int(sizeker/2);
    grow, timeProfile, tac(timeTic);
    //write,format="gpu time for convFFT (only) %f\n",tac(myTime);
    hres=c();
    grow, timeProfile, tac(timeTic);
    a=[];k=[];pad_a=[];pad_spec=[];c=[];
    grow, timeProfile, tac(timeTic);
    "convFFT gpu time profile: "; timeProfile; 
    "convFFT gpu time individual: "; timeProfile(dif); 
  }
  /*
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic=tic();
    ref = float(fft(fft(hdata,1)*fft(hkern,1),-1)/(size)^2);
    write,format="convFFT cpu time : %f\n",tac(timeTic);  
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      "minmax(ref - eclat(hres))";
      minmax(ref - eclat(hres))/max(ref);
    }
  }
  */
  return timeProfile;
}

