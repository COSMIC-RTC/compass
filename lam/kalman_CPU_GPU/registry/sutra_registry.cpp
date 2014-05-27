//sutra_registry.cpp

#include "sutra_registry.h"
#include <stdexcept>

double sutra_registry::get_d(string name) const
{
   map<string, double>::const_iterator it = map_d.find(name);
   if (it == map_d.end())
     throw runtime_error("Error | sutra_registry::get_d | cannot find: " + name);
   return it->second;
}
//                                                                                                
int sutra_registry::get_i(string name) const
{
   map<string, int>::const_iterator it = map_i.find(name);
   if (it == map_i.end())
     throw runtime_error("Error | sutra_registry::get_i | cannot find: " + name);
   return it->second;
}
//                                                                                                
const carma_host_obj<double>& sutra_registry::get_chd(string name) const
{
   map<string, carma_host_obj<double>*>::const_iterator it = map_chd.find(name);
   if (it == map_chd.end())
     throw runtime_error("Error | sutra_registry::get_chd | cannot find: " + name);
   return *(it->second);
}
//                                                                                                
const carma_host_obj<float>&  sutra_registry::get_chf(string name) const
{
   map<string, carma_host_obj<float>*>::const_iterator it = map_chf.find(name);
   if (it == map_chf.end())
     throw runtime_error("Error | sutra_registry::get_chf | cannot find: " + name);
   return *(it->second);
}
//                                                                                                
void sutra_registry::set_chd(string name, const carma_host_obj<double>& chd)
{
   map_chd[name] = new carma_host_obj<double>(&chd);
}
//                                                                                                
void sutra_registry::set_chf(string name, const carma_host_obj<float>& chf)
{
   map_chf[name] = new carma_host_obj<float>(&chf);
}
//                                                                                                
