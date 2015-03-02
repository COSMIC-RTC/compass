require,"yoga.i";
require,"util_fr.i";

if(is_func(diag)==0) { //if yutils is not installed
  include, ["func diag(a) {return 1;}"]
}

func bench_evd(n, niter)
{
  tmp = yoga_obj(random(n, n)-0.5 + diag(random(n)));
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);

  d_U = yoga_obj(mat*0.);
  h_EV = array(0., n);
  write, "doing yoga_syevd, d_mat, h_EV, d_U... ";
  tps = array(0., niter);
  for(i=1; i<=niter; i++) {
    tic; yoga_syevd, d_mat, h_EV, d_U; tps(i)=tac();
    write, format="%d/%d\r", i, niter;
  }
  write, "done";

  write, "[min, avg, max]";
  [min(tps),avg(tps),max(tps)];
  fits_write, swrite(format="bench_evd_gcc_%d.fits", n), tps;
  return tps;
}

func check_getri(n, compare_yorick=)
{
  if(is_void(compare_yorick)) compare_yorick=0;

  tmp = yoga_obj(random(n, n)-0.5 + diag(random(n)));
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);
  tmp = [];

  if(compare_yorick){
    write, format="%s", "doing i_mat= LUsolve(mat)... ";
    tic; i_mat= LUsolve(mat); tps1=tac();
    write, format="in %0.3fs\n", tps1;
  }

  write, "\ntest with double";
  write, format="%s", "doing yoga_getri, d_mat... ";
  tic; yoga_getri, d_mat; tps2=tac();
  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
    write, "max(abs(d_mat() - i_mat))";
    max(abs(d_mat() - i_mat));
  } else {
    write, format="in %0.3fs\n", tps2;
  }
  "max(abs( mat(,+)*d_mat(+,) - unit(n) ));";
  //max(abs(mat(,+)*d_mat()(+,)-unit(n)));
  max(abs( yoga_mm(yoga_obj(mat), d_mat, 'n', 'n') - unit(n) )); 
}

func check_potri(n, compare_yorick=, niter=)
{
  if(is_void(compare_yorick)) compare_yorick=0;
  if(is_void(niter)) niter=1;

  write, format="\ngeneration for the matrix %dx%d... ", n, n;
  tic; 
  tmp = yoga_obj(random(n, n)-0.5 + diag(random(n)));
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);
  tmp = [];
  tps1=tac()
  write, format="in %0.3fs\n", tps1;

  if(compare_yorick){
    tps1=0.;
    write, format="%s", "doing i_mat= LUsolve(mat)... ";
    for(iter=0; iter<niter; iter++){
      tic; i_mat= LUsolve(mat); tps1+=tac()/niter;
    }
    write, format="in %0.3fs\n", tps1;
  }

  write, format="%s", "doing yoga_potri, d_mat... ";
  tps2=0.;
  for(iter=0; iter<niter; iter++){
    d_mat = yoga_obj(mat);
    tic; yoga_potri, d_mat; tps2+=tac()/niter;
  }
  

  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
    write, "max(abs(d_mat() - i_mat))";
    max(abs(d_mat() - i_mat));
    window, 0;
    pli, d_mat();
    window, 1;
    pli, i_mat;
  } else {
    write, format="in %0.3fs\n", tps2;
  }

  "max(abs( mat(,+)*inv_mat(+,) - unit(n) ));";
  max(abs( yoga_mm(yoga_obj(mat), d_mat, 'n', 'n')() - unit(n) ));
  /*
  h_tmp=mat;
  write, "\ntest yoga_potri with double (multiGPU)";
  write, format="%s", "doing yoga_potri_mgpu, 1, mat, d_mat... ";
  tic; yoga_potri_mgpu, 1, h_tmp, d_mat; tps3=tac();

  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps3, tps1/tps3;
    write, "max(abs(d_mat() - i_mat))";
    max(abs(d_mat() - i_mat));
    window, 0;
    pli, d_mat();
    window, 1;
    pli, i_mat;
  } else {
    write, format="in %0.3fs\n", tps3;
  }

  tmp2=  yoga_mm(tmp, tmp, 'n', 't');
  tmp3=  yoga_mm(tmp2, d_mat, 'n', 't');
  "max(abs( mat(,+)*d_mat()(+,) - unit(n) ));";
  max(abs(tmp3()-unit(n)));
  */
  write, format="%s", "doing yoga_potri, h_mat... with CPU ";
  tps2=0.;
  for(iter=0; iter<niter; iter++){
    h_mat = mat;
    tic; yoga_potri, h_mat; tps2+=tac()/niter;
  }
  

  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
    write, "max(abs(h_mat - i_mat))";
    max(abs(h_mat - i_mat));
    window, 0;
    pli, h_mat;
    window, 1;
    pli, i_mat;
  } else {
    write, format="in %0.3fs\n", tps2;
  }

  "max(abs( mat(,+)*inv_mat(+,) - unit(n) ));";
  max(abs( yoga_mm(yoga_obj(mat), yoga_obj(h_mat), 'n', 'n')() - unit(n) ));
}

func check_syevd(n, compare_yorick=)
{
  if(is_void(compare_yorick)) compare_yorick=0;

  write, format="\ngeneration for the matrix %dx%d... ", n, n;
  tic; 
  tmp = yoga_obj(float(random(n, n)-0.5));
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);
  tmp = [];
  tps1=tac()
  write, format="in %0.3fs\n", tps1;

  if(compare_yorick){
    write, format="%s", "doing y_EV = SVdec(mat)... ";
    tic; y_EV = SVdec(mat); tps1=tac();
    write, format="in %0.3fs\n", tps1;
  }

  write, "\ntest with double";
  d_U = yoga_obj(float(mat*0.));
  h_EV = array(0.f, n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, d_U... ";
  tic; yoga_syevd, d_mat, h_EV, d_U; tps2=tac();
  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
    write, "max(abs(h_EV(::-1) - y_EV))";
    max(abs(h_EV(::-1) - y_EV));
  } else {
    write, format="in %0.3fs\n", tps2;
  }
  write, "Verif: max(abs( d_mat - d_U*unit(n)*h_EV(,-)*d_Ut";
  max(abs( mat - yoga_mm(yoga_mm(d_U, yoga_obj(float(unit(n)*h_EV(,-))), 'n', 'n'), d_U, 'n', 't')() ));
  
  write, "\ntest with double (inplace)";
  h_EV2 = array(0.f, n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV2... ";
  tic; yoga_syevd, d_mat, h_EV2; tps2=tac();

  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
    write, "max(abs(h_EV2(::-1) - y_EV))";
    max(abs(h_EV2(::-1) - y_EV));
  } else {
    write, format="in %0.3fs\n", tps2;
  }
  write, "Verif: max(abs( d_mat - d_U*unit(n)*h_EV(,-)*d_Ut";
  max(abs( mat - yoga_mm(yoga_mm(d_mat, yoga_obj(float(unit(n)*h_EV(,-))), 'n', 'n'), d_U, 'n', 't')() ));


  write, "\ntest with double (inplace, noComputeU)";
  d_mat = yoga_obj(float(mat));
  h_EV3 = array(0.f, n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV2... ";
  tic; yoga_syevd, d_mat, h_EV3, noComputeU=1; tps2=tac();

  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
    write, "max(abs(h_EV3(::-1) - y_EV))";
    max(abs(h_EV3(::-1) - y_EV));
  } else {
    write, format="in %0.3fs\n", tps2;
  }
  write, "Verif: max(abs(h_EV - h_EV3))";
  max(abs(h_EV - h_EV3));

}

func check_syevd_m(n, ngpu, compare_yorick=)
{
  if(is_void(ngpu)) ngpu=1;
  if(is_void(compare_yorick)) compare_yorick=0;
  
  write, format="\ngeneration for the matrix %dx%d... ", n, n;
  tic; 
  tmp = yoga_obj(random(n, n)-0.5 + diag(random(n)));
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);
  tmp = [];
  tps1=tac()
  write, format="in %0.3fs\n", tps1;

  if(compare_yorick){
    write, format="%s", "doing y_EV = SVdec(mat)... ";
    tic; y_EV = SVdec(mat); tps1=tac();
    write, format="in %0.3fs\n", tps1;
  }

  write, "\ntest with double";
  h_mat  = mat;
  h_U  = mat*0.;
  h_EV = array(0., n);
  write, format="doing yoga_syevd_m, ngpu=%d, mat, h_EV, h_U... ", ngpu;
  tic; yoga_syevd_m, ngpu, h_mat, h_EV, h_U; tps2=tac();
  if(compare_yorick){
    write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
    write, "max(abs(h_EV(::-1) - y_EV))";
    max(abs(h_EV(::-1) - y_EV));
    fma;
    plg,h_EV(::-1),color="red",  marks=0,width=4;
    plg,y_EV      ,color="green",marks=0,width=4;
    pltitle,"eigenvalues";
    logxy,0,1;
  } else {
    write, format="in %0.3fs\n", tps2;
  }

  write, "max(abs( mat - h_U*unit(n)*h_EV(,-)*h_Ut";
  d_U=yoga_obj(h_U);
  max(abs( mat - yoga_mm(yoga_mm(d_U, yoga_obj(unit(n)*h_EV(,-)), 'n', 'n'), d_U, 'n', 't')() ));
}

func compare_syevd(n, ngpu)
{
  if(is_void(ngpu)) ngpu=1;
  
  write, format="\ngeneration for the matrix %dx%d... ", n, n;
  tic; 
  tmp = yoga_obj(random(n, n)-0.5 + diag(random(n)));
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);
  tmp = [];
  tps1=tac()
  write, format="in %0.3fs\n", tps1;

  write, "\ntest singleGPU with double";
  d_U = yoga_obj(mat*0.);
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, d_U... ";
  tic; yoga_syevd, d_mat, h_EV, d_U; tps1=tac();
  write, format="in %0.3fs\n", tps1;

  write, "\ntest singleGPU with double inplace";
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV... ";
  tic; yoga_syevd, d_mat, h_EV; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  d_mat = yoga_obj(mat);
  write, "\ntest singleGPU with double (without U computation)";
  d_U = yoga_obj(mat*0.);
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, d_U, noComputeU=1... ";
  tic; yoga_syevd, d_mat, h_EV, d_U, noComputeU=1; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "\ntest singleGPU with double inplace (without U computation)";
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, noComputeU=1... ";
  tic; yoga_syevd, d_mat, h_EV, noComputeU=1; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "\ntest multiGPU with double";
  h_mat  = mat;
  h_U  = mat*0.;
  h_EV2 = array(0., n);
  write, format="doing yoga_syevd_m, ngpu=%d, mat, h_EV, h_U... ", ngpu;
  tic; yoga_syevd_m, ngpu, h_mat, h_EV2, h_U; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "\ntest multiGPU with double (without U computation)";
  h_mat  = mat;
  h_U  = mat*0.;
  h_EV2 = array(0., n);
  write, format="doing yoga_syevd_m, ngpu=%d, mat, h_EV, h_U, noComputeU=1... ", ngpu;
  tic; yoga_syevd_m, ngpu, h_mat, h_EV2, h_U, noComputeU=1; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "max(abs(h_EV - h_EV2))";
  max(abs(h_EV - h_EV2));

  write, "\ntest MKL with double";
  h_mat  = mat;
  h_U  = mat;
  h_EV2 = array(0., n);
  write, format="%s", "doing yoga_syevd, mat, h_EV, h_U... ";
  tic; yoga_syevd, h_U, h_EV2; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "\ntest MKL with double (without U computation)";
  h_mat  = mat;
  h_U  = mat*0.;
  h_EV2 = array(0., n);
  write, format="%s", "doing yoga_syevd, mat, h_EV, h_U, noComputeU=1... ";
  tic; yoga_syevd, h_mat, h_EV2, h_U, noComputeU=1; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "max(abs(h_EV - h_EV2))";
  max(abs(h_EV - h_EV2));

}

func check_svd(n,m)
{
  //n=4096; m=128;
  min_mn = min([m,n]);

  y_mat=random_n(m,n);

  write,"computing svd on CPU: y_S = SVdec(y_mat, y_U, y_Vt)";
  tic; y_S = SVdec(y_mat, y_U, y_Vt); tps1 = tac();
  write,format= "done in %0.3fs...\n", tps1;

//  d_mat=yoga_obj(y_mat);
//  d_U=yoga_obj(array(0.0,m,m));
//  d_Vt=yoga_obj(array(0.0,n,n));
//  d_S=yoga_obj(array(0.0,min_mn));
//  write,"computing svd on GPU: yoga_svd, d_mat, d_S, d_U, d_Vt";
//  tic; yoga_svd,d_mat,d_S,d_U,d_Vt; tps2 = tac();
//  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

//  write, "max(abs(d_S()-y_S))";
//  max(abs(d_S()-y_S));

  h_mat=yoga_host_obj(y_mat,pagelock=1);
  h_U=yoga_host_obj(array(0.0,m,m),pagelock=1);
  h_Vt=yoga_host_obj(array(0.0,n,n),pagelock=1);
  h_S=yoga_host_obj(array(0.0,min_mn),pagelock=1);
  write,"computing svd on GPU: yoga_svd_host,h_mat,h_S,h_U, h_Vt";
  tic; yoga_svd_host,h_mat, h_S, h_U, h_Vt; tps3 = tac();
  write, format="in %0.3fs (x%0.3f)\n", tps3, tps1/tps3;

  write, "max(abs(h_S()-y_S))";
  max(abs(h_S()-y_S));

  fma;
  //plg,d_S(),              marks=0,width=4;
  plg,h_S(),color="red",  marks=0,width=4;
  plg,y_S  ,color="green",marks=0,width=4;
  pltitle,"eigenvalues";
  logxy,0,1;
  write,format = "GPU full time : %f\n",tps3;
  write,format = "GPU svd only  : %f\n",tps2;
  write,format = "CPU svd time  : %f\n",tps1;
    
}

func check_mm_cpu(n,m,k){
	write, format="generation of A %dx%d\n", n, m;
	A=random_n(n, m);
	write, format="generation of B %dx%d\n", m, k;
	B=random_n(m, k);
	write, format="%s", "doing C=A(, +)*B(+,)...";
	tic; C=A(, +)*B(+,); tps1=tac();
	write, format=" in %0.3fs\n", tps1;
	write, format="%s", "doing C=yoga_mm_cpu(A,B)...";
	tic; C1=yoga_mm_cpu(A,B); tps2=tac();
	write, format=" in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
	write, "Verif : max(abs(C-C1))";
	max(abs(C-C1));
}

func check_axpy_cpu(n,m){
	write, format="generation of X %dx%d\n", n, m;
	X=random_n(n, m);
	write, format="generation of Y %dx%d\n", n, m;
	Y=random_n(n, m);
	write, format="%s", "doing C=X+Y...";
	tic; C=X+Y; tps1=tac();
	write, format=" in %0.3fs\n", tps1;
	write, format="%s", "doing C=yoga_axpy_cpu(1,X,Y)...";
	tic; C1=yoga_axpy_cpu(1,X,Y); tps2=tac();
	write, format=" in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;
	write, "Verif : max(abs(C-C1))";
	max(abs(C-C1));
}
/*
mat = fits_read("cmaa24240.fits");
n=24240
h_U  = mat*0.;
h_EV = array(0., n);
ngpu=2
tic; yoga_syevd_m, ngpu, mat, h_EV, h_U; tps2=tac();
tps2
584.435



mat = fits_read("cmaa37920.fits");
n=37920
ngpu=2
h_U  = mat*0.;
h_EV = array(0., n);
tic; yoga_syevd_m, ngpu, mat, h_EV, h_U; tps2=tac();
tps2


hippo2:~$ nvidia-smi 
Thu Feb 27 17:48:03 2014       
+------------------------------------------------------+                       
| NVIDIA-SMI 5.319.37   Driver Version: 319.37         |                       
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla M2090         Off  | 0000:11:00.0     Off |                  Off |
| N/A   N/A    P0   145W /  N/A |     6016MB /  6143MB |     98%      Default |
+-------------------------------+----------------------+----------------------+
|   1  Tesla M2090         Off  | 0000:14:00.0     Off |                  Off |
| N/A   N/A    P0   142W /  N/A |     6016MB /  6143MB |     98%      Default |
+-------------------------------+----------------------+----------------------+
                                                                               
+-----------------------------------------------------------------------------+
| Compute processes:                                               GPU Memory |
|  GPU       PID  Process name                                     Usage      |
|=============================================================================|
|    0     23380  yorick                                             11876MB  |
|    1     23380  yorick                                             11876MB  |
+-----------------------------------------------------------------------------+


 */
