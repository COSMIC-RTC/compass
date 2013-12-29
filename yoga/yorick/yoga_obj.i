require,"yoga.i";
require,"util_fr.i";

/*
  ____ _   _ _____ ____ _  ______  
 / ___| | | | ____/ ___| |/ / ___| 
| |   | |_| |  _|| |   | ' /\___ \ 
| |___|  _  | |__| |___| . \ ___) |
 \____|_| |_|_____\____|_|\_\____/ 
                                   
 */

func check_yoga_obj(size)
{
  if (is_void(size)) size= 512;

  checkall_yoga_cublas;
  
  check_transpose,size;  
}


func checkall_yoga_cublas(size)
{
  write,"BLAS1 checks";
  check_cublas1;

  write,"BLAS2 checks";
  check_cublas2;
  
  write,"BLAS3 checks";
  check_cublas3;

  "";
  write,"checkMV";
  "";
  checkMV,size;
  "";
  write,"checkMM";
  "";
  checkMM,size;
}

func check_cublas1(size)
{
/* DOCUMENT check_cublas1(size)
     this templates shows how to exercise yoga cublas1 fonctions
   SEE ALSO:
 */
  if (is_void(size)) size= 128;
  
// create object
  a = yoga_random("float",[1,size*size]);
// cpu version of the data  
  test = a();
  
// check max
  ind1 = yoga_imax(a);
  ind2 = wheremax(test);

  write,format="error for max :%d\n",ind1-ind2;
     
// check min
  ind1 = yoga_imin(a);
  ind2 = wheremin(test);

  write,format="error for min :%d\n",ind1-ind2;

// check asum
  res1 = yoga_asum(a);
  res2 = abs(test)(*)(sum);

  write,format="error for asum :%f\n",res1-res2;

// check nrm2
  res1 = yoga_nrm2(a);
  res2 = sqrt(test(*)(+)*test(*)(+));
  write,format="error for nrm2 :%f\n",res1-res2;

// check scale
  test1 = test;

  yoga_scale,a,10.f;
  test = a();
  
  test1*=10.f;
  write,format="error for scale: %f\n",max(abs(test-test1));

// check swap
  b = yoga_obj("float",[2,size,size]);
  yoga_swap,b,a;
  test=b();
  write,format="error for swap: %f\n",max(abs(test(*)-test1));

// check copy
  c = yoga_obj("float",[2,size,size]);
  yoga_copy,c,b;
  test=c();
  write,format="error for copy: %f\n",max(abs(test(*)-test1));

//check axpy
  a = yoga_random("float",[1,size*size]);
  b = yoga_obj(a);
  yoga_scale,b,1.4f;
  
  data1 = a();
  data2 = 1.4f*a();

  yoga_axpy,b,10.0f,a;
  res = data1*10.+data2;
  test=b();
  write,format="error for axpy : %f\n",max(abs(test-res))/max(res);
  
//check dot
  res1 = yoga_dot(b,b);
  res2 = float(test(*)(+)*test(*)(+));
  write,format="error for dot : %f\n",(res1-res2)/res2;
  
// delete objects
  a = b = [];
}

func check_cublas2(sizem,sizen)
{
/* DOCUMENT check_cublas1(size)
     this templates shows how to exercise yoga cublas2 fonctions
   SEE ALSO:
 */
  if (sizem==[]) sizem = 512;
  if (sizen==[]) sizen = 1024;

// input data
  matA = float(random(sizem,sizen));
  vectx = float(random(sizen));
  vecty = float(random(sizem));

// create objects
  matA_gpu = yoga_obj(matA);
  vectx_gpu = yoga_obj(vectx);
  vecty_gpu = yoga_obj(vecty);
  
// compute gemv
  yoga_mv,vecty_gpu,matA_gpu,vectx_gpu;
  
// transfer back for checks

  res2 = matA(,+)*vectx(+);
  res2 += vecty;

  write,format="error for gemv : %f\n",max(abs(vecty_gpu()-res2))/max(abs(res2));

// compute ger
  matA_gpu = yoga_obj(matA);
  vectx_gpu = yoga_obj(vectx);
  vecty_gpu = yoga_obj(vecty);
  
  yoga_rank1,matA_gpu,vecty_gpu,vectx_gpu;
  
  res2 = vecty(,-::sizen-1)(,+)*vectx(-::sizen-1,)(+,);
  res2 /= sizen;
  res2 += matA;
  write,format="error for ger : %f\n",max(abs(matA_gpu()(*)-res2(*)))/max(abs(res2));
  
// delete objects
  matA_gpu =vectx_gpu = vecty_gpu = [];

// create symmetric matrix
  matA = float(dist(sizem));
  vectx = float(random(sizem));
  vecty = array(0.0f,sizem);

// create objects
  matA_gpu = yoga_obj(matA);
  vectx_gpu = yoga_obj(vectx);
  vecty_gpu = yoga_obj(vecty);

  yoga_symv,vecty_gpu,matA_gpu,vectx_gpu;
  
// transfer back for checks

  res2 = matA(,+)*vectx(+);
  res2 += vecty;

  write,format="error for symv : %f\n",max(abs(vecty_gpu()-res2))/max(abs(res2));
  
  
}

func check_cublas3(sizem,sizen,sizek)
{
/* DOCUMENT check_cublas3(sizem,sizen,sizek)
     this templates shows how to exercise yoga cublas3 fonctions
   SEE ALSO:
 */

  if (sizem==[]) sizem = 256;
  if (sizen==[]) sizen = 512;
  if (sizek==[]) sizek = 1024;

  //check gemm
  matA = float(random(sizem,sizek));
  matB = float(random(sizek,sizen));
  matC = float(random(sizem,sizen));

// create objects
/*
  matA_gpu = yoga_obj(matA);
  matB_gpu = yoga_obj(matB);
  matC_gpu = yoga_obj(matC);
*/
  matA_gpu = yoga_setm(matA);
  matB_gpu = yoga_setm(matB);
  matC_gpu = yoga_setm(matC);
  
// compute gemm
  yoga_mm,matC_gpu,matA_gpu,matB_gpu;
  
  res2 = matA(,+)*matB(+,);
  res2 += matC*0.0f;
  
  write,format="gemm diff : %f\n",max(abs(matC_gpu()(*)-res2(*)))/max(res2);
  
// delete objects
  matA_gpu = matB_gpu = matC_gpu = [];
  
//check symm
  matA = float(dist(sizek));
  matB = float(random(sizek,sizen));
  matC = array(0.0f,sizek,sizen);
  
  matA_gpu = yoga_setm(matA);
  matB_gpu = yoga_setm(matB);
  matC_gpu = yoga_setm(matC);
  
// compute symm
  yoga_symm,matC_gpu,matA_gpu,matB_gpu;
  
  res2 = matA(,+)*matB(+,);
  res2 += matC;
  
  write,format="symm diff : %f\n",max(abs(matC_gpu()(*)-res2(*)))/max(res2);
  
// delete objects
  matA_gpu = matB_gpu = matC_gpu = [];

//check dmm
  matA = float(random(sizek,sizen));
  vectx = float(random(sizek))
  matC = array(0.0f,sizek,sizen);
  
  //matA_gpu = yoga_setm(float(transpose(matA)));
  matA_gpu = yoga_setm(matA);
  vectx_gpu = yoga_setv(vectx);
  matC_gpu = yoga_setm(matC);
  
// compute dmm
  yoga_dmm,matC_gpu,matA_gpu,vectx_gpu;
  tmp = array(0.0f,sizek,sizek);
  tmp(*)(1::sizek+1) = vectx;
  res2 = tmp(,+)*matA(+,);
  //res2 = matA(+,)*matA(+,); // op == 't'
  res2 += matC;
  
  write,format="left dmm diff : %f\n",max(abs(matC_gpu()(*)-res2(*)))/max(res2(*));
  
  vecty = float(random(sizen))
  matC = array(0.0f,sizek,sizen);
  
  vecty_gpu = yoga_setv(vecty);
  matC_gpu = yoga_setm(matC);
  
  yoga_dmm,matC_gpu,matA_gpu,vecty_gpu,'r';

  tmp = array(0.0f,sizen,sizen);
  tmp(*)(1::sizen+1) = vecty;
  res2 = matA(,+)*tmp(+,);
  //res2 = matA(+,)*matA(+,); // op == 't'
  res2 += matC;

  write,format="right dmm diff : %f\n",max(abs(matC_gpu()(*)-res2(*)))/max(res2(*));
// delete objects
  matA_gpu = matC_gpu = vectx_gpu = ecty_gpu = [];

//check syrk
  //matA = fits_read("imat.fits");
  //sizek = dimsof(matA)(2);
  //sizen = dimsof(matA)(3);
  matA = float(random(sizen,sizen));
  matC = array(0.0f,sizen,sizen);
  
  //matA_gpu = yoga_setm(float(transpose(matA)));
  matA_gpu = yoga_setm(matA);
  matC_gpu = yoga_setm(matC);
  
// compute syrk
  yoga_syrk,matC_gpu,matA_gpu;

  msk = where(matC_gpu() != 0);
  
  res2 = matA(,+)*matA(,+);
  //res2 = matA(+,)*matA(+,); // op == 't'
  res2 += matC;
  
  write,format="syrk diff : %f\n",max(abs(matC_gpu()(msk)-res2(msk)))/max(res2(msk));
  
// delete objects
  matA_gpu = matC_gpu = [];

//check syrk
  matA = fits_read("imat.fits");
  sizek = dimsof(matA)(2);
  sizen = dimsof(matA)(3);
  //matA = float(random(sizen,sizen));
  matC = array(0.0f,sizek,sizek); 
  
  //matA_gpu = yoga_setm(float(transpose(matA)));
  matA_gpu = yoga_obj(float(matA));
  matC_gpu = yoga_obj(matC);
  
// compute syrk
  //yoga_syrk,matC_gpu,matA_gpu;

  yoga_mm,matC_gpu,matA_gpu,matA_gpu,'n','t';
  
  msk = where(matC_gpu() != 0);
  
  res2 = matA(,+)*matA(,+);
  //res2 = matA(+,)*matA(+,); // op == 't'
  res2 += matC;
  
  write,format="syrk (with mm) diff : %f\n",max(abs(matC_gpu()(*)-res2(*)))/max(res2(*));
// delete objects
  matA_gpu = matC_gpu = [];

}


func checkMV(sizem, sizen,ctype=)
/* DOCUMENT checkMV(sizem, sizen,ctype=)
     this templates does the profiling of matrix-vector multiply
   SEE ALSO:
 */
{
  if (is_void(ctype)) ctype = "all";
  if(is_void(sizem)) sizem=512;
  if(is_void(sizen)) sizen=sizem;
  
  matA=float(random(sizem, sizen));  vectx=float(random(sizen));
  vecty=array(0.f, sizem);

  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    
    matA_gpu = yoga_obj(matA);
    grow, timeProfile, tac(timeTic);
    vectx_gpu = yoga_obj(vectx);
    grow, timeProfile, tac(timeTic);
    vecty_gpu = yoga_obj(vecty);
    grow, timeProfile, tac(timeTic);
    
    yoga_mv,vecty_gpu,matA_gpu,vectx_gpu;
    grow, timeProfile, tac(timeTic);
  
    vectyb=vecty_gpu();
    grow, timeProfile, tac(timeTic);

    matA_gpu = vectx_gpu = vecty_gpu = [];
    grow, timeProfile, tac(timeTic);
    "";
    "MV mult gpu time profile: ";
    "h2dm     h2dv1    h2dv2    comp    d2hv    free"  
    timeProfile; 
    "";
    "MV mult gpu time individual: ";
    "h2dv1       h2dv2       comp        d2hv       free"  
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic=tic();
    vecty=matA(,+)*vectx(+);
    write,format="MV mult cpu time : %f\n",tac(timeTic);
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      timecpu=tac(timeTic);
      write,format="error for mv : %f\n",max(abs(vecty-vectyb))/max(abs(vecty));
      tmp = timeProfile;
      "";
      write,format="Acceleration factor w memcopy : %.1f\n",timecpu/tmp(1:4)(sum);
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(dif)(3);
    }
  }
  return timeProfile;
}

func checkMM(sizem, sizen, sizek,ctype=)
/* DOCUMENT checkMM(sizem, sizen, sizek,ctype=)
     this templates does the profiling of matrix-matrix multiply
   SEE ALSO:
 */
{
  if (is_void(ctype)) ctype = "all";
  if(is_void(sizem)) sizem=512;
  if(is_void(sizen)) sizen=sizem;
  if(is_void(sizek)) sizek=sizem;
  
  matA=float(random(sizem, sizek)); matB=float(random(sizek, sizen));
  matC=matCb=array(float, dimsof(matA)(2), dimsof(matB)(3));

  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    
    matA_gpu = yoga_obj(matA);
    grow, timeProfile, tac(timeTic);
    matB_gpu = yoga_obj(matB);
    grow, timeProfile, tac(timeTic);
    matC_gpu = yoga_obj(matC);
    grow, timeProfile, tac(timeTic);
  
    yoga_mm,matC_gpu,matA_gpu,matB_gpu;
    grow, timeProfile, tac(timeTic);
  
    matCb=matC_gpu();
    grow, timeProfile, tac(timeTic);
    
    matA_gpu = matB_gpu = matC_gpu;
    grow, timeProfile, tac(timeTic);
    
    "";
    "MM mult gpu time profile: ";
    "h2dm     h2dv1    h2dv2    comp    d2hv    free"  
    timeProfile; 
    "";
    "MM mult gpu time individual: ";
    "h2dm2       h2dm3       comp        d2hm       free"  
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    "matD=matA(,+)*matB(+,)";
    timeTic=tic();
    matD=matA(,+)*matB(+,); 
    write,format="MM cpu time : %f\n",tac(timeTic);

    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      write,format="error for mm : %f\n",max(abs(matCb-matD))/max(abs(matD));
      timecpu=tac(timeTic);
      tmp = timeProfile;
      "";
      write,format="Acceleration factor w memcopy : %.1f\n",timecpu/tmp(1:4)(sum);
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(dif)(3);
    }
  }
  return timeProfile;
}

func check_transpose(size)
{
  if (size == []) size = 512;

  // check sumCU
  test =float(random(size,size));
  obj1 = yoga_obj(test);
  yoga_host2device,obj1,test;
  obj2 = yoga_obj("float",[2,size,size]);
  timeTic = tic();
  for (i=1;i<=100;i++)
    yoga_transpose,obj2,obj1;
  tmpGPU = tac(timeTic)/100.;
  res1 = obj2();
  res2 = transpose(test);

  write,format="error for transpose %f\n",max(abs(res2-res1))/max(abs(res2));

  obj1 = obj2 = [];
  
  write,format="gpu time for transpose %f\n",tmpGPU;
  
  checkTranspose,size;
}

func checkTranspose(size,ctype=) {
  if (is_void(type)) ctype = "all";
  if(is_void(size)) size=1024;
  A = float(random_n(size,size));
  
  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    obj1 = yoga_obj(A);
    obj2 = yoga_obj("float",[2,size,size]);
    grow, timeProfile, tac(timeTic);
    yoga_transpose,obj2,obj1;
    grow, timeProfile, tac(timeTic);
    sout1 = obj2();
    grow, timeProfile, tac(timeTic); 
    obj1 = obj2 = [];
    grow, timeProfile, tac(timeTic);
    "";
    "Transpose gpu time profile: ";
    "h2dm       comp      d2hm      free"  
    timeProfile; 
    "";
    "Transpose gpu time individual: ";
    "comp        d2hm       free"  
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic=tic();
    sout2=transpose(A);
    write,format="cpu time for transpose %f\n",tac(timeTic);
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      timecpu=tac(timeTic);
      write,format="error for transpose %f\n",max(abs(sout1-sout2))/max(abs(sout2));
      tmp = timeProfile;
      "";
      write,format="Acceleration factor w memcopy : %.1f\n",timecpu/tmp(1:3)(sum);
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(dif)(1);
    }
  }

  return timeProfile;
}

/*

func check_generic(size) {
  
  if (size == []) size = 1024;
  tmp = float(random(size*size));
  
  obj1 = yoga_create("float",[2,size,size]);
  obj2 = yoga_create("float",[2,size,size]);
  
  yoga_host2device,obj1,tmp;
  _generic2CUS,*obj1.handler,*obj2.handler; // warmup !

  res1 = tmp*0.0f;
  yoga_host2device,obj1,tmp;
  //_generic2CUS,*obj1.handler,*obj2.handler; // warmup !

  timeTic = tic();
  _generic2CUS,*obj1.handler,*obj2.handler;
  tmpGPU = tac(timeTic);

  yoga_device2host,obj2,res1;
  write,format="kernel error : %.8f\n",max(abs(res1-sin(2.0f*tmp)));
  write,format="gpu time for kernel2d %f\n",tmpGPU;

  yoga_destroy,obj1;
  yoga_destroy,obj2;

  obj1 = yoga_create("float",[1,size*size]);
  obj2 = yoga_create("float",[1,size*size]);
  yoga_host2device,obj1,tmp;
  _genericCUS,*obj1.handler,*obj2.handler; // warmup !
  res1 = tmp*0.0f;
  
  yoga_host2device,obj1,tmp;
  
  timeTic = tic();
  _genericCUS,*obj1.handler,*obj2.handler;
  tmpGPU = tac(timeTic);
  
  yoga_device2host,obj2,res1;
  write,format="kernel error : %.8f\n",max(abs(res1-sin(2.0f*tmp)));
  write,format="gpu time for kernel1d %f\n",tmpGPU;
  
  yoga_destroy,obj1;
  yoga_destroy,obj2;
  
  timeTic = tic();
  tmp2 = sin(2.0f*tmp);
  write,format="cpu ref time  %f\n",tac(timeTic);

}


func check_sum(size)
{
  if (size == []) size = 512;

  // check sumCU
  test =float(random(size,size));
  obj = yoga_create("float",[2,size,size]);
  yoga_host2device,obj,test;
  res1 = sumCU(obj);
  res2 = test(*)(sum);

  write,format="sum result GPU %f\n",res1;
  write,format="sum result CPU %f\n",res2;
  write,format="error %f\n",abs(res2-res1)/res2;

  yoga_destroy,obj;
  
  checkSum,size;
}


func checkSum(size,ctype=)
{
  if (is_void(ctype)) ctype = "all";
  if(is_void(size)) size=1024;
  test = float(random_n(size,size));
  res1=0.f;

  if ((ctype == "all") || (ctype == "gpu")) {
    timeTic = tic();
    timeProfile=[];
    obj = yoga_create("float",[2,size,size]);
    grow, timeProfile, tac(timeTic);
    yoga_host2device,obj,test;
    grow, timeProfile, tac(timeTic);
    res1 = sumCU(obj);
    grow, timeProfile, tac(timeTic);
    yoga_destroy,obj;
    grow, timeProfile, tac(timeTic);
    "";
    "Sum gpu time profile: ";
    "Alloc     h2dm     comp    free"  
    timeProfile; 
    "";
    "Sum gpu time individual: ";
    "h2dm        comp        free"  
    timeProfile(dif); 
    "";
  }
  
  if ((ctype == "all") || (ctype == "cpu")) {
    timeTic = tic();
    res2 = test(*)(sum);
    write,format="sum cpu time : %f\n",tac(timeTic);  
    if (ctype != "all") timeProfile = tac(timeTic);
    else {
      timecpu=tac(timeTic);
      write,format="error for sum %f\n",abs(res1-res2)/res2;
      tmp = timeProfile(dif);
      "";
      write,format="Acceleration factor w memcopy : %.1f\n",timecpu/tmp(1:2)(sum);
      write,format="Acceleration factor w/ memcopy : %.1f\n",timecpu/tmp(2);
    }
  }

  return timeProfile;
}

func check_transfer(size)
{
  if (is_void(size)) size= 100;
  
  test = float(indgen(size)); // the input data
  
// create first object
  obj1 = yoga_create("float",[1,size]);

// fill with data  
  yoga_host2device,obj1,test;
     
// create second object
  obj2 = yoga_create("float",[1,size]);
  
// transfer data !
  yoga_device2device,obj2,obj1,size;
  
// get data in second object (for checks)
  test2 = array(0.f,size);
  yoga_device2host,obj2,test2;
  
  write,format="error for transfer :%f\n",max(abs(test2-test));
  
// delete objects
  yoga_destroy,obj1;
  yoga_destroy,obj2;
}

*/
