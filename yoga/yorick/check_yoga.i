require,"yoga.i";
require,"util_fr.i";
require,"yoga_obj.i";
require,"yoga_fft.i";
require,"yoga_sort.i";
require,"yoga_scan.i";
require,"yoga_rng.i";
require,"yoga_magma.i";

func checkAll(size)
{  
  if(is_void(size)) size=512;
  
  check_yoga_fft,size;
  
  check_yoga_obj,size;

  checkall_yoga_prng,size;

  //check_yoga_sort,size;

  //check_yoga_scan,size;
}



/*
*/
