//sp_parmap.h

/*  Copyright (C) 2003-2008 Sergey Rodionov
 *
 *  Author: Sergey Rodionov (seger@astro.spbu.ru)
 *
 *  This software is released under the GNU General Public License (GNU GPL).
 *  Please read the included file COPYING for more information.
*/

//sky predictor parameters map
//we work parameters like NEMO work with parameteres
//(defv value)

#ifndef __SEGER__SP_PARMAP_H__
#define __SEGER__SP_PARMAP_H__

#include <vector>
#include <string>
#include <map>
using namespace std;


const char SP_PARMAP_PARFILE[] = "parfile";

class sp_parmap
  //you should use this class in init() function
{
 public:
   //TODO: 
   //gcc 4.3.2 can't convert from "char** to const char**"
   //(from not const to const!!!!)
   sp_parmap(int argc, const char** argv, const char*defv[]);
   sp_parmap(int argc,       char** argv, const char*defv[]);

   string get_s(const string&);
   double get_d(const string&);
   int    get_i(const string&);
   bool   get_b(const string&); //boolean value (0|1 t|f true|false )
   vector<string>  get_vs(const string&, int quan = -1, char sep = ',');
   vector<double>  get_vd(const string&, int quan = -1, char sep = ',');
   vector<int>     get_vi(const string&, int quan = -1, char sep = ',');
   // if quan<0 --> quantity not fixed

   vector<string>  get_vs_fix(const string& key, int quan, char sep = ':')
     { return get_vs(key,quan,sep); }; 
   vector<double>  get_vd_fix(const string& key, int quan, char sep = ':')
     { return get_vd(key,quan,sep); }; 
   vector<int>     get_vi_fix(const string& key, int quan, char sep = ':')
     { return get_vi(key,quan,sep); }; 
   
   string get_full_arg()
     {return full_arg;};
   
 private:
   map<string,string> m;
   string full_arg;
   string pname; //name of program
 private:
   void _init(int argc, const char** argv, const char*defv[]);
   
   void load_defv(const char* defv[] , 
		  vector<string>& defv_param,
		  vector<string>& defv_defv,
		  vector<string>& defv_desc);
     //last string in defv MUST be NULL
     //pre last string in defv MUST be VERSION
   void add_to_map(string key,string param);  //if duplicate --> error
   bool separate(string in, string& out1, string& out2, char sep);
   void print_help(const vector<string>& defv_param);
   void print_full_help(const vector<string> & defv_param,
			const vector<string>& defv_desc);
   
   void error_exit();
};

#endif
