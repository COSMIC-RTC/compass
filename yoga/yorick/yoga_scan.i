require,"yoga.i";

func check_yoga_scan(size)
{
  if (size == []) size = 128;

  // check sumCU
  "";
  write,"Check Scan";
  "";
  test = float(random([1,size*size]));
  obj1 = yoga_obj(test);
  res2=yoga_add(obj1);
  //res2 = scanCU(obj1);
  res1 = test(*)(sum);
  write,format="scan add error %f\n",res2-res1;

  res2 = yoga_max(obj1);
  res1 = max(test(*));
  write,format="scan max error %f\n",res2-res1;
  
  res2 = yoga_min(obj1);
  res1 = min(test(*));
  write,format="scan min error %f\n",res2-res1;

  checkScan,size;
}


/*
  ____ _               _     ____                  
 / ___| |__   ___  ___| | __/ ___|  ___ __ _ _ __  
| |   | '_ \ / _ \/ __| |/ /\___ \ / __/ _` | '_ \ 
| |___| | | |  __/ (__|   <  ___) | (_| (_| | | | |
 \____|_| |_|\___|\___|_|\_\|____/ \___\__,_|_| |_|
                                                   
*/  
func checkScan(size,type=,ctype=) {
  
  if (is_void(ctype)) ctype = "all";
  if(is_void(size)) size=1024;
  if(is_void(type)) type = "sum";
  
  A = float(random_n(size*size));
  
  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    obj1 = yoga_obj(A);
    grow, timeProfile, tac(timeTic);
    if (type == "max") res1 = yoga_max(obj1);
    if (type == "min") res1 = yoga_min(obj1);
    if (type == "sum") res1 = yoga_add(obj1);
    grow, timeProfile, tac(timeTic);
    obj1 = [];
    grow, timeProfile, tac(timeTic);
    "";
    "Scan gpu time profile: ";
    "Alloc     comp    free"  
    timeProfile; 
    "";
    "Scan gpu time individual: ";
    "comp        free"  
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic=tic();
    if (type == "sum") sout2=A(sum);
    if (type == "max") sout2=max(A);
    if (type == "min") sout2=min(A);
    write,format="cpu time for scan %f\n",tac(timeTic);
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      timecpu=tac(timeTic);
      write,format="scan error : %f\n",(res1-sout2)/sout2;
      tmp = timeProfile(dif);
      "";
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(1);
    }
  }

  return timeProfile;
}
