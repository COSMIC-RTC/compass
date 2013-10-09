require,"yoga_aolib.i";

func check_aotemplate(dim)
{

  if (dim == []) dim = 4096;

  // create aotemplate
  tmp = yoga_aotemplate(dim);
  // fill with random
  yoga_templatefill,tmp;
  // retreive aotemplate
  res=yoga_gettemplate(tmp);
  // do computation
  yoga_templatecomp,tmp;
  // retreive result
  res2=yoga_gettemplate(tmp,"res");
  // check result
  sin((res-roll(res,-1))*2*3.14159)(1:100)-res2(1:100);
  // delete template
  tmp = [];
}
