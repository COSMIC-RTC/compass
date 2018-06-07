// sp_parmap.cpp

/*  Copyright (C) 2003-2008 Sergey Rodionov
 *
 *  Author: Sergey Rodionov (seger@astro.spbu.ru)
 *
 *  This software is released under the GNU General Public License (GNU GPL).
 *  Please read the included file COPYING for more information.
 */

#include "sp_parmap.h"
#include <algorithm>
#include <cctype>
#include <iostream>
using namespace std;

sp_parmap::sp_parmap(int argc, const char** argv, const char* defv[]) {
  _init(argc, argv, defv);
}
//
sp_parmap::sp_parmap(int argc, char** argv, const char* defv[]) {
  _init(argc, (const char**)argv, defv);
}
//
string sp_parmap::get_s(const string& key) {
  map<string, string>::const_iterator it = m.find(key);
  if (it == m.end()) {
    cerr << "Error. Bad parameters: Can't find parameter \"" << key << "\""
         << endl;
    error_exit();
  }
  return it->second;
}
//
double sp_parmap::get_d(const string& key) {
  string s = get_s(key);
  char* uk;
  double rez = strtod(s.c_str(), &uk);
  if (*uk != 0 || s.empty()) {
    cerr << "Error. Bad float-point parameters \"" << key << "=" << s << "\""
         << endl;
    error_exit();
  }
  return rez;
}
//
int sp_parmap::get_i(const string& key) {
  string s = get_s(key);
  char* uk;
  int rez = strtol(s.c_str(), &uk, 10);
  if (*uk != 0 || s.empty()) {
    cerr << "Error. Bad integer parameters \"" << key << "=" << s << "\""
         << endl;
    error_exit();
  }
  return rez;
}
//
bool sp_parmap::get_b(const string& key) {
  string s = get_s(key);
  for (size_t i = 0; i < s.size(); i++) s[i] = tolower(s[i]);
  if (s == "t" || s == "true" || s == "1") return true;
  if (s == "f" || s == "false" || s == "0") return false;
  cerr << "Error. Bad boolean parameter \"" << key << "=" << s << "\"" << endl;
  cerr << "(must be \"t\"|\"f\"|\"true\"|\"false\" \"0\"|\"1\")" << endl;
  exit(EXIT_FAILURE);
  return false;
}
//
vector<string> sp_parmap::get_vs(const string& key, int quan, char sep) {
  vector<string> vec;
  vec.resize(0);
  string tmp = get_s(key);
  vec.push_back("");
  for (size_t i = 0; i < tmp.size(); i++) {
    if (tmp[i] == sep)
      vec.push_back("");
    else
      vec[vec.size() - 1] += tmp[i];
  }
  if (quan >= 0 && quan != (int)vec.size()) {
    cerr << "Bad number of vector of parameters \"" << key << "=" << tmp << "\""
         << endl;
    exit(EXIT_FAILURE);
  }
  return vec;
}
//
vector<double> sp_parmap::get_vd(const string& key, int quan, char sep) {
  vector<string> s_vec = get_vs(key, quan, sep);
  vector<double> d_vec;
  d_vec.resize(s_vec.size());
  for (size_t i = 0; i < d_vec.size(); i++) {
    char* uk;
    d_vec[i] = strtod(s_vec[i].c_str(), &uk);
    if (*uk != 0 || s_vec[i].empty()) {
      cerr << "Error. Bad float-point parameters in vector: \"" << key << "="
           << get_s(key) << "\"" << endl;
      exit(EXIT_FAILURE);
    }
  }
  return d_vec;
}
//
vector<int> sp_parmap::get_vi(const string& key, int quan, char sep) {
  vector<string> s_vec = get_vs(key, quan, sep);
  vector<int> i_vec;
  i_vec.resize(s_vec.size());
  for (size_t i = 0; i < i_vec.size(); i++) {
    char* uk;
    i_vec[i] = strtol(s_vec[i].c_str(), &uk, 10);
    if (*uk != 0 || s_vec[i].empty()) {
      cerr << "Error. Bad inteter parameters in vector: \"" << key << "="
           << get_s(key) << "\"" << endl;
      exit(EXIT_FAILURE);
    }
  }
  return i_vec;
}
//
void sp_parmap::_init(int argc, const char** argv, const char* defv[]) {
  if (argc < 1) {
    cerr << "Error | sp_parmap::sp_parmap | argc < 1" << endl;
    exit(EXIT_FAILURE);
  }
  pname = argv[0];

  vector<string> defv_param;
  vector<string> defv_defv;
  vector<string> defv_desc;
  load_defv(defv, defv_param, defv_defv, defv_desc);

  bool det_help = false;
  bool det_help_full = false;
  // load parameters from argumentrs
  m.clear();
  bool det_unnamed_zone = true;
  for (int i = 1; i < argc; i++) {
    string arg = argv[i];
    string arg1, arg2;
    if (separate(arg, arg1, arg2, '='))  // arg -- contaion =
    {
      det_unnamed_zone = false;
    } else {
      if (!det_unnamed_zone)  // error
      {
        cerr << "Error then load parameters: Ambiguous keyword =" << endl;
        error_exit();
      }
    }
    if (det_unnamed_zone) {
      if ((i - 1) >= (int)defv_param.size()) {
        cerr << "Error then load parameters: "
                "Too many un-named arguments"
             << endl;
        error_exit();
      }
      if (arg == "help")
        det_help = true;
      else
        add_to_map(defv_param[i - 1], arg);
    } else {
      if (arg1 == "help") {
        if (arg2.empty())
          det_help = true;
        else
          det_help_full = true;
      } else {
        if (find(defv_param.begin(), defv_param.end(), arg1) ==
            defv_param.end()) {
          cerr << "Error then load parameters: Parameter \"" << arg1
               << "\" unknown" << endl;
          error_exit();
        }
        add_to_map(arg1, arg2);
      }
    }
  }
  // get default parameters
  for (size_t i = 0; i < defv_param.size(); i++)
    if (!defv_defv[i].empty())  // if default value exist
    {
      if (m.find(defv_param[i]) == m.end())  // this parameters not load
      {
        add_to_map(defv_param[i], defv_defv[i]);
      }
    }
  // print help if we need it
  if (det_help_full) print_full_help(defv_param, defv_desc);
  if (det_help) print_help(defv_param);
  // create full_arg
  for (int i = 0; i < argc; i++) {
    full_arg += argv[i];
    if (i != argc - 1) full_arg += " ";
  }
}
//
void sp_parmap::load_defv(const char* defv[], vector<string>& defv_param_key,
                          vector<string>& defv_param_defv,
                          vector<string>& defv_desc) {
  defv_param_key.clear();
  defv_param_defv.clear();
  defv_desc.clear();
  int i = 0;
  while (defv[i] != NULL) {
    string param_key_defv, param_desc, param_key, param_defv;
    if (!separate(defv[i], param_key_defv, param_desc, '\n')) {
      cerr << "Error then load parameters: bad defv (can't find \\n)" << endl;
      exit(EXIT_FAILURE);
    }
    if (!separate(param_key_defv, param_key, param_defv, '=')) {
      cerr << "Error then load parameters: can't find = symbol" << endl;
      exit(EXIT_FAILURE);
    }
    param_desc.assign(param_desc, param_desc.find_first_not_of(" "),
                      string::npos);
    defv_param_key.push_back(param_key);
    defv_param_defv.push_back(param_defv);
    defv_desc.push_back(param_desc);
    i++;
  }
  //   if ( defv_param_key.back() != "VERSION" )
  //     {
  //	cerr<<"Error then load parameters: in defv last parameter "<<
  //	  "MUST be \"VERSION\""<<endl;
  //   }
}
//
void sp_parmap::add_to_map(string key, string param) {
  if (m.find(key) != m.end())  // find this key
  {
    cerr << "Error then load parameters: Parameter \"" << key << "\" duplicated"
         << endl;
    exit(EXIT_FAILURE);
  }
  m[key] = param;
}
//
bool sp_parmap::separate(string in, string& out1, string& out2, char sep) {
  size_t pos = in.find(sep);
  if (pos == string::npos) return false;
  out1.assign(in, 0, pos);
  out2.assign(in, pos + 1, string::npos);
  return true;
}
//
void sp_parmap::print_full_help(const vector<string>& defv_param,
                                const vector<string>& defv_desc) {
  size_t max_param_size = 0;
  size_t max_desc_size = 0;
  for (size_t i = 0; i < defv_param.size(); i++) {
    max_param_size = max(defv_param[i].size(), max_param_size);
    max_desc_size = max(defv_desc[i].size(), max_desc_size);
  }
  for (size_t i = 0; i < defv_param.size(); i++) {
    cout << defv_param[i]
         << string(max_param_size - defv_param[i].size() + 5, ' ') << ": "
         << defv_desc[i] << string(max_desc_size - defv_desc[i].size() + 3, ' ')
         << "[";
    if (m.find(defv_param[i]) != m.end()) cout << m[defv_param[i]];
    cout << "]" << endl;
  }
  exit(EXIT_SUCCESS);
}
//
void sp_parmap::print_help(const vector<string>& defv_param) {
  cout << pname;
  for (size_t i = 0; i < defv_param.size(); i++) {
    cout << "  " << defv_param[i] << "=";
    if (m.find(defv_param[i]) != m.end()) cout << m[defv_param[i]];
  }
  cout << endl;
  exit(EXIT_SUCCESS);
}
//
void sp_parmap::error_exit() {
  cout << "Type \"" << pname << " help=h\" for help" << endl;
  cout << "and  \"" << pname << " help=\" for get default parameters" << endl;
  exit(EXIT_FAILURE);
}
