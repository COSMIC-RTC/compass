require,"yoga.i";

require,"util_fr.i"

  func check_cula_svd(n,m,hostonly=,verif=)
{
//	imat = fits_read("imat.fits");
//	n=dimsof(imat)(2);
//	m=dimsof(imat)(3);

  if (hostonly == []) hostonly = 0;
  if (verif == []) verif = 1;
  //n=4096; m=128;
  min_mn = min([m,n]);
	
//	tmp=imat
  tmp=float(random(m,n));

  if (!hostonly) {
    mat=yoga_obj(tmp);
    U=yoga_obj(array(0.0f,m,m));
    Vt=yoga_obj(array(0.0f,n,n));
    S=yoga_obj(array(0.0f,min_mn));
    write,"computing svd on GPU: yoga_svd,mat,S,U, Vt";
    tic; yoga_cula_svd,mat,S,U,Vt; tps1 = tac();
    write,"done ...";
  }
  if (0) {
  //tmp=float(random(m,n));
  mat=yoga_host_obj(tmp,pagelock=1);
  if (!verif) tmp=[];
  U=yoga_host_obj(array(0.0f,m,m),pagelock=1);
  Vt=yoga_host_obj(array(0.0f,n,n),pagelock=1);
  S1=yoga_host_obj(array(0.0f,min_mn),pagelock=1);
  write,"computing svd on GPU: yoga_svd,mat,S1,U, Vt";
  tic; yoga_cula_svd_host,mat,S1,U,Vt; tps2 = tac();
  write,"done ...";
  }
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
//write, "max(abs(U()-U2))";
//max(abs(U()-U2));
//write, "max(abs(Vt()-V2t))";
//max(abs(Vt()-V2t));
  if (verif) {
    write,"computing svd on CPU: S2 = SVdec(mat(), U2, V2t)";
    tic; S2 = SVdec(tmp, U2, V2t); tps3 = tac();
    write,"done ...";
    fma;
    plg,S(),marks=0,width=4;
    plg,S1(),color="red",marks=0,width=4;
    plg,S2,color="green",marks=0,width=4;
    pltitle,"eigenvalues";
    logxy,0,1;
    write,format = "GPU culaDevice time : %f\n",tps1;
    write,format = "GPU culaHost time  : %f\n",tps2;
    write,format = "CPU svd time  : %f\n",tps3;
  }
  if (!hostonly) {
    write, "max(abs(S()-S2))";
    max(abs(S()-S2));
    write, "max(abs(S1()-S2))";
    //max(abs(S1()-S2));
  }
}

//check_cula_svd, 512,512
