require,"yoga.i";

func check_yoga_sort(size)
{
  if (size == []) size = 128;

  // check sumCU
  "";
  write,"Check Sort";
  "";
  test =float(random(size,size));
  res2 = sort(test(*));
  obj = yoga_obj(test);
  res1 = (yoga_sort(obj))()+1;
  write,format="sorted keys error %f\n",max(abs(test(res2)-obj()(*)));
  write,format="sorted values error %d\n",max(abs(res2-res1(*)));

  checkSort,size*size;
}


/*
      _               _      ____             _   
  ___| |__   ___  ___| | __ / ___|  ___  _ __| |_ 
 / __| '_ \ / _ \/ __| |/ / \___ \ / _ \| '__| __|
| (__| | | |  __/ (__|   <   ___) | (_) | |  | |_ 
 \___|_| |_|\___|\___|_|\_\ |____/ \___/|_|   \__|
                                                  
*/  
func checkSort(size,ctype=) {
  
  if (is_void(ctype)) ctype = "all";
  if(is_void(size)) size=1024;
  
  A = float(random_n(size));
  
  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    obj = yoga_obj(A);
    grow, timeProfile, tac(timeTic);
    ind1=yoga_sort(obj)()(*);
    grow, timeProfile, tac(timeTic);
    sout1=obj()(*);
    grow, timeProfile, tac(timeTic); 
    obj = [];
    grow, timeProfile, tac(timeTic); 
    "";
    "Sort gpu time profile: ";
    "Alloc     comp    d2hm    free"  
    timeProfile; 
    "";
    "Sort gpu time individual: ";
    "comp        d2hm       free"  
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic=tic();
    ind2 = sort(A);
    sout2=A(ind2);
    write,format="cpu time for sort %f\n",tac(timeTic);
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      timecpu=tac(timeTic);
      ind1+=1;
      write,format="%f %f %d %d\n",min(sout1-sout2),max(sout1-sout2),min(ind1-ind2),max(ind1-ind2);
      tmp = timeProfile(dif);
      "";
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(1);
    }
  }

  return timeProfile;
}
