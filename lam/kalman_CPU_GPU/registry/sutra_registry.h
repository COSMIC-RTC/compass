//sutra_registry.h


#ifndef __SEGER__SUTRA_REGISTRY_H__
#define __SEGER__SUTRA_REGISTRY_H__

#include <carma_host_obj.h>
#include <string>
#include <map>
using namespace std;

class sutra_registry
{
 public:
   //functions for get parameters
   double                        get_d(string name) const;
   int                           get_i(string name) const;   
   const carma_host_obj<double>& get_chd(string name) const;
   const carma_host_obj<float>&  get_chf(string name) const;
     
   //function for set parameters  
   void set_d(string name, double d) {map_d[name] = d;};
   void set_i(string name, int    i) {map_i[name] = i;};
   
   //we copy matrices!
   void set_chd(string name, const carma_host_obj<double>&);
   void set_chf(string name, const carma_host_obj<float>&);
   
   //destructor (in future we remove it by not using pointers)
   ~sutra_registry();
 private:
   map<string, double> map_d;
   map<string, int>    map_i;
   map<string,  carma_host_obj<double>* > map_chd;
   map<string,  carma_host_obj<float>* >  map_chf;
};



#endif
