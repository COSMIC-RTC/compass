require,"yoga.i";

require,"util_fr.i";

func bench_evd(n, niter)
{
  tmp = yoga_obj(random(n, n)-0.5);
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

func check_getri(n)
{
  tmp = yoga_obj(random(n, n)-0.5);
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);

  write, format="%s", "doing i_mat= LUsolve(mat)... ";
  tic; i_mat= LUsolve(mat); tps1=tac();
  write, format="in %0.3fs\n", tps1;

  write, "\ntest with double";
  write, format="%s", "doing yoga_getri, d_mat... ";
  tic; yoga_getri, d_mat; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "max(abs(d_mat() - i_mat))";
  max(abs(d_mat() - i_mat));

  "max(abs( mat(,+)*d_mat()(+,) - unit(n) ));"
  max(abs(mat(,+)*d_mat()(+,)-unit(n)));
}

func check_potri(n, compare_yorick=, niter=)
{
  if(is_void(compare_yorick)) compare_yorick=0
  if(is_void(niter)) niter=1

  tmp = yoga_obj(random(n, n)-0.5);
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);

  if(compare_yorick){
    tps1=0.;
    write, format="%s", "doing i_mat= LUsolve(mat)... ";
    for(iter=0; iter<niter; iter++){
      tic; i_mat= LUsolve(mat); tps1+=tac()/niter;
    }
    write, format="in %0.3fs\n", tps1;
  }

  write, "\ntest yoga_potri with double";
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

  tmp2=  yoga_mm(tmp, tmp, 'n', 't');
  tmp3=  yoga_mm(tmp2, d_mat, 'n', 't');
  "max(abs( mat(,+)*d_mat()(+,) - unit(n) ));";
  max(abs(tmp3()-unit(n)));
/*
  h_tmp=mat;
  write, "\ntest yoga_potri with double (multiGPU)";
  write, format="%s", "doing yoga_potri_mgpu, 2, mat, d_mat... ";
  tic; yoga_potri_mgpu, 2, h_tmp, d_mat; tps3=tac();

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
}

func check_syevd(n)
{
  tmp = yoga_obj(random(n, n)-0.5);
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);

  imat=fits_read("imat.fits");
  mat=imat(,+)*imat(,+);
  d_mat=yoga_obj(mat);
  n=dimsof(mat)(2);
  
  write, format="%s", "doing y_EV = SVdec(mat)... ";
  tic; y_EV = SVdec(mat); tps1=tac();
  write, format="in %0.3fs\n", tps1;

  write, "\ntest with double";
  d_U = yoga_obj(mat*0.f);
  h_EV = array(0.f, n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, d_U... ";
  tic; yoga_syevd, d_mat, h_EV, d_U; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "max(abs(h_EV(::-1) - y_EV))";
  max(abs(h_EV(::-1) - y_EV));

  write, "\ntest with double inplace";
  //d_U = yoga_obj(mat*0.);
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV... ";
  tic; yoga_syevd, d_mat, h_EV; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "max(abs(h_EV(::-1) - y_EV))";
  max(abs(h_EV(::-1) - y_EV));
/*
  fma;
  plg,h_EV(::-1),color="red",  marks=0,width=4;
  plg,y_EV      ,color="green",marks=0,width=4;
  pltitle,"eigenvalues";
  logxy,0,1;
*/
}

func check_syevd_m(n, ngpu)
{
  if(is_void(ngpu)) ngpu=1;
  
  tmp = yoga_obj(random(n, n)-0.5);
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);

  write, format="%s", "doing y_EV = SVdec(mat)... ";
  tic; y_EV = SVdec(mat); tps1=tac();
  write, format="in %0.3fs\n", tps2;

  write, "\ntest with double";
  h_U  = mat*0.;
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd_m, ngpu, mat, h_EV, h_U... ";
  tic; yoga_syevd_m, ngpu, mat, h_EV, h_U; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "max(abs(h_EV(::-1) - y_EV))";
  max(abs(h_EV(::-1) - y_EV));

  fma;
  plg,h_EV(::-1),color="red",  marks=0,width=4;
  plg,y_EV      ,color="green",marks=0,width=4;
  pltitle,"eigenvalues";
  logxy,0,1;
}

func compare_syevd(n, ngpu)
{
  if(is_void(ngpu)) ngpu=1;
  
  tmp = yoga_obj(random(n, n)-0.5);
  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  mat = d_mat(); //tmp(,+)*tmp(,+);

  write, "\ntest with double";
  d_U = yoga_obj(mat*0.);
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, d_U... ";
  tic; yoga_syevd, d_mat, h_EV, d_U; tps1=tac();
  write, format="in %0.3fs\n", tps1;

  write, "\ntest with double inplace";
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV... ";
  tic; yoga_syevd, d_mat, h_EV; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  d_mat = yoga_mm(tmp, tmp, 'n', 't')
  write, "\ntest with double (without U computation)";
  d_U = yoga_obj(mat*0.);
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, d_U, noComputeU=1... ";
  tic; yoga_syevd, d_mat, h_EV, d_U, noComputeU=1; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "\ntest with double inplace (without U computation)";
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, noComputeU=1... ";
  tic; yoga_syevd, d_mat, h_EV, noComputeU=1; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "\ntest with double";
  h_U  = mat*0.;
  h_EV2 = array(0., n);
  write, format="%s", "doing yoga_syevd_m, ngpu, mat, h_EV, h_U... ";
  tic; yoga_syevd_m, ngpu, mat, h_EV2, h_U; tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "\ntest with double (without U computation)";
  h_U  = mat*0.;
  h_EV2 = array(0., n);
  write, format="%s", "doing yoga_syevd_m, ngpu, mat, h_EV, h_U, noComputeU=1... ";
  tic; yoga_syevd_m, ngpu, mat, h_EV2, h_U, noComputeU=1; tps2=tac();
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

  d_mat=yoga_obj(y_mat);
  d_U=yoga_obj(array(0.0,m,m));
  d_Vt=yoga_obj(array(0.0,n,n));
  d_S=yoga_obj(array(0.0,min_mn));
  //write,"computing svd on GPU: yoga_svd, d_mat, d_S, d_U, d_Vt";
  //tic; yoga_svd,d_mat,d_S,d_U,d_Vt; tps2 = tac();
  //write, format="in %0.3fs (x%0.3f)\n", tps2, tps1/tps2;

  write, "max(abs(d_S()-y_S))";
  max(abs(d_S()-y_S));

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
