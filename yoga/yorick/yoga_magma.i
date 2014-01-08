require,"yoga.i";

require,"util_fr.i";

func check_evd(n)
{
  tmp = random(n, 128);
  tmp2= tmp(,+)*tmp(,+);
  d_mat = yoga_obj(tmp2);
  d_U = yoga_obj(tmp2*0.);
  h_EV = array(0., n);
  write, format="%s", "doing yoga_syevd, d_mat, h_EV, d_U... ";
  tic; yoga_syevd, d_mat, h_EV, d_U; tps1=tac();
  write, format="in %0.3fs\n", tps1;

  write, format="%s", "doing y_EV = SVdec(tmp2)... ";
  tic; y_EV = SVdec(tmp2); tps2=tac();
  write, format="in %0.3fs (x%0.3f)\n", tps2, tps2/tps1;

  write, "max(abs(h_EV(::-1) - y_EV))";
  max(abs(h_EV(::-1) - y_EV));
}

func check_svd(n,m)
{
  //n=4096; m=128;
  min_mn = min([m,n]);

  tmp=float(random(m,n));
  mat=yoga_obj(tmp);
  U=yoga_obj(array(0.0f,m,m));
  Vt=yoga_obj(array(0.0f,n,n));
  S=yoga_obj(array(0.0f,min_mn));
  write,"computing svd on GPU: yoga_svd,mat,S,U, Vt";
  tic; yoga_svd,mat,S,U,Vt; tps1 = tac();
  write,"done ...";
  //tmp=float(random(m,n));
  mat=yoga_host_obj(tmp,pagelock=1);
  U=yoga_host_obj(array(0.0f,m,m),pagelock=1);
  Vt=yoga_host_obj(array(0.0f,n,n),pagelock=1);
  S1=yoga_host_obj(array(0.0f,min_mn),pagelock=1);
  write,"computing svd on GPU: yoga_svd,mat,S1,U, Vt";
  tic; yoga_svd_host,mat,S1,U,Vt; tps2 = tac();
  write,"done ...";
  write,"computing svd on CPU: S2 = SVdec(mat(), U2, V2t)";
  tic; S2 = SVdec(tmp, U2, V2t); tps3 = tac();
  write,"done ...";

//write, "S";
//info, S();
//info, S1();
//info, S2;
//write, "U";
//info, U();
//info, U2;
//write, "Vt";
//info, Vt();
//info, V2t;

  write, "max(abs(S()-S2))";
  max(abs(S()-S2));
  write, "max(abs(S1()-S2))";
  max(abs(S1()-S2));
//write, "max(abs(U()-U2))";
//max(abs(U()-U2));
//write, "max(abs(Vt()-V2t))";
//max(abs(Vt()-V2t));
  fma;
  plg,S(),marks=0,width=4;
  plg,S1(),color="red",marks=0,width=4;
  plg,S2,color="green",marks=0,width=4;
  pltitle,"eigenvalues";
  logxy,0,1;
  write,format = "GPU full time : %f\n",tps1;
  write,format = "GPU svd only  : %f\n",tps2;
  write,format = "CPU svd time  : %f\n",tps3;
    
}
