require,"yoga.i";
require,"util_fr.i";

/*
  ____ _   _ _____ ____ _  ______  
 / ___| | | | ____/ ___| |/ / ___| 
| |   | |_| |  _|| |   | ' /\___ \ 
| |___|  _  | |__| |___| . \ ___) |
 \____|_| |_|_____\____|_|\_\____/ 
                                   
 */

func checkall_yoga_prng(size)
{
  if (size == []) size = 1024;

  // check random
  seed = 12345;
  a = yoga_random("float",[2,size,size]);  // create object
  res1 = a();
  res2 =random(size,size);

  "uniform result GPU min / max / rms / avg";
  write,format="%f %f %f %f\n",min(res1),max(res1),res1(*)(rms),res1(*)(avg);
  "uniform result CPU min / max / rms / avg";
  write,format="%f %f %f %f\n",min(res2),max(res2),res2(*)(rms),res2(*)(avg);

  // check randomn
  a = yoga_random_n("float",[2,size,size]);  // create object
  res1 = a();

  res2 =random_n(size,size);

  "normal result GPU min / max / rms / avg";
  write,format="%f %f %f %f\n",min(res1),max(res1),res1(*)(rms),res1(*)(avg);
  "normal result CPU min / max / rms / avg";
  write,format="%f %f %f %f\n",min(res2),max(res2),res2(*)(rms),res2(*)(avg);

  a=[];
  
  checkRNG,size,'U';

  checkRNG,size,'N';

}

func checkRNG(size,type,ctype=)
{
  if (is_void(ctype)) ctype = "all";
  if(is_void(type)) type='U';
  
  if (type == 'U') write,"Uniform Distribution";
  if (type == 'N') write,"Normal Distribution";
  
  if(is_void(size)) size=1024;
  tmp=array(float,size,size);
  
  seed=12345;
  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    if (type == 'U') a = yoga_random("float",[2,size,size]);  // create object
    if (type == 'N') a = yoga_random_n("float",[2,size,size]);  // create object
    grow, timeProfile, tac(timeTic);
    if (type == 'P') {
      xx = indgen(size*size) - size*size/2.;
      tmp2 = float(100.*exp(-(xx^2/(size*sqrt(size/2))^2)));
      //plg,tmp2,marks=0;
      a = yoga_obj(tmp2);            // if poisson copy image
      grow, timeProfile, tac(timeTic);
      tmp = tmp2*0.0f;
    }  
    if (type == 'U') yoga_random,a;             // if uniform randomCU
    if (type == 'N') yoga_random_n,a;            // if gauss random_nCU
    if (type == 'P') yoga_poisson,a;            // if gauss random_nCU

    grow, timeProfile, tac(timeTic);
    
    tmp = a();
    grow, timeProfile, tac(timeTic);
    
    a=[];

    grow, timeProfile, tac(timeTic);
    "";
    "RNG gpu time profile: ";
    "Alloc     comp    d2hm    free"  ;
    timeProfile; 
    "";
    "RNG gpu time individual: ";
    "comp        d2hm       free"  ;
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic = tic();
    if (type == 'U') tmp2=random(size,size);
    if (type == 'N') tmp2=random_n(size,size);
    if (type == 'P') {
      xx = indgen(size*size) - size*size/2.;
      tmp2 = float(100.*exp(-(xx^2/(size*sqrt(size/2))^2)));
      _poidev,tmp2,size*size;
    }
    write,format="RNG cpu time : %f\n",tac(timeTic);  
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      timecpu=tac(timeTic);
      "result GPU min / max / rms / avg";
      write,format="%f %f %f %f\n",min(tmp),max(tmp),tmp(*)(rms),tmp(*)(avg);
      "result CPU min / max / rms / avg";
      write,format="%f %f %f %f\n",min(tmp2),max(tmp2),tmp2(*)(rms),tmp2(*)(avg);
      tmp = timeProfile(dif);
      "";
      write,format="Acceleration factor w memcopy : %.1f\n",timecpu/tmp(1:2)(sum);
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(1);
    }
  }
  return timeProfile;
}
