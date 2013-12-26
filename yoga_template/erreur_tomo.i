require,"yoga_template.i";

func check_eigen(l,&nmf)
// checking eigen values
{
  if ( nmf==-1 ) {
    nn = where(l<0.8);
    if(is_array(nn)) {
      nmf = numberof(l)-min(nn)+1;
    } else {
      nmf = 0;
    }
    write,format="Number of filtered modes = %d \n",nmf;
  }
  if( disp==1) {
    window,0;
    fma; limits; logxy,0,1;
    plg, l;
    plg,[l(1),l(0)],numberof(l)-[nmf,nmf],marks=0,color="red";
  }
}

func invert_eigen( L, nfilt )
// invert eigen values vector
{
    if( nfilt>0 ) {
        L1 = L;
        L1(1-nfilt:0) = 1.0;
        L1 = 1./L1;
        L1(1-nfilt:0) = 0.0;
    } else {
        L1 = 1.0 / L;
    }
    return L1;
}


func erreur_tomo(lib=,nmf=,filename=)
{
  // choice of the library to perform computations
  if ((lib == "") || (lib == [])) lib = "yorick";
  // choice of the number of filtered modes
  if (nmf == []) nmf = -1;
  // choice of default file name
  if ((filename == "") || (filename == [])) filename = "-np20-d42";

  // reading covariance matrices
  cmm = fits_read("data/cmm"+filename+".fits");
  cpm = fits_read("data/cpm"+filename+".fits");
  cpp = fits_read("data/cpp"+filename+".fits");
  
  n = dimsof(cmm)(2);
  write,format="starting EVD %dx%d with %s librairy\n",n,n, lib;

  if (lib == "yorick"){
    l = SVdec(cmm, U, VT);
  }  
  if (lib == "magma") {
    N = dimsof(cmm)(2);
    ldda = ((N + 31)/32)*32;
    d_U= yoga_obj(array(0.0,N,ldda));
    l2 = yoga_magma_syevd('V', 'U', cmm, d_U);
    l = SVdec(cmm, U, VT); 
    error;
  } 
  if (lib == "cula") {
    // transfering cov matrix to GPU
    g_cmm = yoga_obj(double(cmm));
    // allocating space on GPU for the result
    g_U   = yoga_obj(double(cmm*0));
    g_VT  = yoga_obj(double(cmm*0));
    g_L   = yoga_obj(double(0*cmm(1,)));

    // compute svd using CULA
    yoga_cula_svd, g_cmm, g_L, g_VT, g_U;

    // switching to double precision
    g_U = yoga_obj(double(g_U()));
    g_VT = yoga_obj(double(g_VT()));
    g_L = yoga_obj(double(g_L()));
    // retrieving eigen values for filtering on CPU
    l = g_L();
  }

  
  // check eigen values
  check_eigen,l,nmf;
  
  write,"Pseudo-inverse computation";
  // inverting eigen values
  l1 = invert_eigen( l, nmf );

  if (lib == "yorick") {
    // computing the pseudo-inverse
    cmm1 =  (U*l1(-,))(,+) * VT(+,);
    
    write,"Reconstructor computation";
    // computing the reconstructor
    R = cpm(,+) * cmm1(+,);
  }

  if (lib == "cula") {
    // computing the pseudo-inverse
    n = numberof(l1);
    tmp = g_U()*l1(-,);
    g_U = [];
    g_m1 = yoga_obj(array(0.0,n,n));
    yoga_mm, g_m1, yoga_obj(tmp), g_VT;
    g_VT = [];
    tmp = [];

    // computing the reconstructor 
    g_cpm = yoga_obj(cpm);
    g_R = yoga_mm( g_cpm, g_m1, 'n', 't', 1.0, 0.0 );
    R = g_R();
  }
  
  return R;
}

