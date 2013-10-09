#include "yoga_template.i"

func test_rico(void)
{
  im1 = yoga_obj(float(random(100)));
  im2 = yoga_obj(float(random(100)));
  
  res = im1()*im2();
  
  yoga_immult,im1,im2;
  
  (res - im1())(*)(ptp);
}

