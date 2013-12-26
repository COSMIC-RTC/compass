require,"yoga_template.i";

func check_inv(n)
{
  if (n==[]) n = 1024;
  n2 = n * n;

  a = random_n(n,n);

  b = LUsolve(a);
  ldda = ((n + 31)/32)*32;
  d_b= yoga_obj(array(0.0,n,ldda));
  yoga_magma_getri,a,d_b;

  write,format="error : %f\n",max(abs(b-d_b()));
}
