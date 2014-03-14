/*! @file yoga_ao_api.h

 @brief Definition of functions callable by the Yorick interpreter

 Based on the standard API for interfacing Yorick packages to the interpreter
 Binding is realized through the definition of the corresponding "extern" wrappers
 in the libyoga_ao.i file (see yorick directory). Functions callable by the interpreter
 start with a Y_

 Several persistent objects are also defined :

 Copyright (c) 2011, Damien Gratadour & Arnaud Sevin

 This file is part of YoGA, the Yorick with GPU Acceleration plugin

 This program is free software; you can redistribute it and/or  modify it
 under the terms of the GNU General Public License  as  published  by the
 Free Software Foundation; either version 3 of the License,  or  (at your
 option) any later version.
 
 This program is distributed in the hope  that  it  will  be  useful, but
 WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _YOGA_AO_API_H_
#define _YOGA_AO_API_H_

#include <yapi.h>

/*
 _                        _       __ 
 | |_ _   _ _ __   ___  __| | ___ / _|
 | __| | | | '_ \ / _ \/ _` |/ _ \ |_ 
 | |_| |_| | |_) |  __/ (_| |  __/  _|
 \__|\__, | .__/ \___|\__,_|\___|_|  
 |___/|_|                        

 */
typedef struct tscreen_struct {
  void *sutra_tscreen;
  int device;
} tscreen_struct;

typedef struct atmos_struct {
  void *sutra_atmos;
  int device;
} atmos_struct;

typedef struct source_struct {
  void *sutra_source;
  int device;
} source_struct;

typedef struct target_struct {
  void *sutra_target;
  int device;
} target_struct;

typedef struct phase_struct {
  void *sutra_phase;
  int device;
} phase_struct;

typedef struct wfs_struct {
  void *sutra_wfs;
  int device;
} wfs_struct;

typedef struct acquisim_struct {
  void *sutra_acquisim;
  int device;
} acquisim_struct;

typedef struct sensors_struct {
  void *sutra_sensors;
  int device;
} sensors_struct;

typedef struct rtc_struct {
  void *sutra_rtc;
  int device;
} rtc_struct;

typedef struct dms_struct {
  void *sutra_dms;
  int device;
} dms_struct;

typedef struct telemetry_struct {
  void *sutra_telemetry;
  int device;
} telemetry_struct;

typedef struct aotemplate_struct {
  void *sutra_aotemplate;
  int device;
} aotemplate_struct;

extern "C" {

  /*
   _   _
   / \ | |_ _ __ ___   ___  ___
   / _ \| __| '_ ` _ \ / _ \/ __|
   / ___ \ |_| | | | | | (_) \__ \
  /_/   \_\__|_| |_| |_|\___/|___/

   */

  void
  atmos_free(void *obj);
  void
  atmos_print(void *obj);

  static y_userobj_t yAtmos = { const_cast<char*>("yAtmos Object"),
  /** 
   * @typedef Yorick API atmos userobj
   */
  &atmos_free, &atmos_print, 0, 0, 0 };

  void
  Y_yoga_atmos(int argc);
  void
  Y_get_spupil(int argc);

  /*
   _   ____                           
   | |_/ ___|  ___ _ __ ___  ___ _ __  
   | __\___ \ / __| '__/ _ \/ _ \ '_ \ 
   | |_ ___) | (__| | |  __/  __/ | | |
   \__|____/ \___|_|  \___|\___|_| |_|

   */

  void
  tscreen_free(void *obj);
  void
  tscreen_print(void *obj);

  static y_userobj_t yTscreen = { const_cast<char*>("yTscreen Object"),
  /** 
   * @typedef Yorick API tscreen userobj
   */
  &tscreen_free, 0, 0, 0, 0 };

  void
  Y_init_tscreen(int argc);
  void
  Y_get_tscreen(int argc);
  void
  Y_extrude_tscreen(int argc);
  void
  Y_get_tscreen_config(int argc);
  void
  Y_get_tscreen_update(int argc);
  void
  Y_tscreen_initvk(int argc);
  void
  Y_tscreen_genevk(int argc);

  /*
   ____                           
   / ___|  ___  _   _ _ __ ___ ___ 
   \___ \ / _ \| | | | '__/ __/ _ \
     ___) | (_) | |_| | | | (_|  __/
   |____/ \___/ \__,_|_|  \___\___|

   */
  void
  source_free(void *obj);
  void
  source_print(void *obj);

  static y_userobj_t ySource = { const_cast<char*>("ySource Object"),
      &source_free, 0, 0, 0, 0 };

  /*
   _____                    _   
   |_   _|_ _ _ __ __ _  ___| |_ 
   | |/ _` | '__/ _` |/ _ \ __|
   | | (_| | | | (_| |  __/ |_ 
   |_|\__,_|_|  \__, |\___|\__|
   |___/          
   */

  void
  target_free(void *obj);
  void
  target_print(void *obj);

  static y_userobj_t yTarget = { const_cast<char*>("yTarget Object"),
      &target_free, &target_print, 0, 0, 0 };

  void
  Y_yoga_target(int argc);
  void
  Y_target_addlayer(int argc);
  void
  Y_target_atmostrace(int argc);
  void
  Y_target_getimage(int argc);
  void
  Y_target_getphase(int argc);
  void
  Y_target_getphasetele(int argc);
  void
  Y_target_getamplipup(int argc);

  /*
   ____  _                    
   |  _ \| |__   __ _ ___  ___ 
   | |_) | '_ \ / _` / __|/ _ \
    |  __/| | | | (_| \__ \  __/
   |_|   |_| |_|\__,_|___/\___|

   */

  void
  phase_free(void *obj);
  void
  phase_print(void *obj);
  void
  phase_eval(void *obj, int n);

  static y_userobj_t yPhase = { const_cast<char*>("yPhase Object"), &phase_free,
      &phase_print, &phase_eval, 0, 0 };

  void
  Y_yoga_phase(int argc);
  void
  Y_phase_set(int argc);

  /*
   __        _______ ____  
   \ \      / /  ___/ ___| 
   \ \ /\ / /| |_  \___ \ 
   \ V  V / |  _|  ___) |
   \_/\_/  |_|   |____/ 

   */

  void
  wfs_free(void *obj);
  void
  wfs_print(void *obj);

  static y_userobj_t yWfs = { const_cast<char*>("yWfs Object"), &wfs_free,
      &wfs_print, 0, 0, 0 };

  void
  Y_yoga_wfs(int argc);
  void
  Y_wfs_initgs(int argc);

  /*
   *                        _     _
   *   __ _  ___ __ _ _   _(_)___(_)_ __ ___
   *  / _` |/ __/ _` | | | | / __| | '_ ` _ \
   * | (_| | (_| (_| | |_| | \__ \ | | | | | |
   *  \__,_|\___\__, |\__,_|_|___/_|_| |_| |_|
   *               |_|
   */

  void
  acquisim_free(void *obj);
  void
  acquisim_print(void *obj);

  static y_userobj_t yAcquisim = { const_cast<char*>("yAcquisim Object"),
      &acquisim_free, &acquisim_print, 0, 0, 0 };

  void
  Y_yoga_acquisim(int argc);
  void
  Y_acquisim_fillbcube(int argc);

  /*
   _       _                     _
   | |_ ___| | ___ _ __ ___   ___| |_ _ __ _   _
   | __/ _ \ |/ _ \ '_ ` _ \ / _ \ __| '__| | | |
   | ||  __/ |  __/ | | | | |  __/ |_| |  | |_| |
   \__\___|_|\___|_| |_| |_|\___|\__|_|   \__, |
   |___/
   */

  void
  telemetry_free(void *obj);
  void
  telemetry_print(void *obj);

  static y_userobj_t yTelemetry = { const_cast<char*>("yTelemetry Object"),
      &telemetry_free, &telemetry_print, 0, 0, 0 };

  void
  Y_yoga_telemetry(int argc);

  /*
   ____                                
   / ___|  ___ _ __  ___  ___  _ __ ___ 
   \___ \ / _ \ '_ \/ __|/ _ \| '__/ __|
   ___) |  __/ | | \__ \ (_) | |  \__ \
    |____/ \___|_| |_|___/\___/|_|  |___/

   */

  void
  sensors_free(void *obj);
  void
  sensors_print(void *obj);

  static y_userobj_t ySensors = { const_cast<char*>("ySensors Object"),
      &sensors_free, &sensors_print, 0, 0, 0 };

  void
  Y_yoga_sensors(int argc);
  void
  Y_sensors_initgs(int argc);
  void
  Y_sensors_addlayer(int argc);
  void
  Y_sensors_initarr(int argc);
  void
  Y_sensors_loadkernels(int argc);
  void
  Y_sensors_initlgs(int argc);
  void
  Y_sensors_updatelgsprof(int argc);
  void
  Y_sensors_updatelgs(int argc);
  void
  Y_sensors_compimg(int argc);
  void
  Y_sensors_compimg_tele(int argc);
  void
  Y_sensors_getimg(int argc);
  void
  Y_sensors_getdata(int argc);
  void
  Y_sensors_setphase(int argc);

  /*
   ____  __  __     
   |  _ \|  \/  |___ 
   | | | | |\/| / __|
   | |_| | |  | \__ \
    |____/|_|  |_|___/

   */

  void
  dms_free(void *obj);
  void
  dms_print(void *obj);

  static y_userobj_t yDMs = { const_cast<char*>("yDMs Object"), &dms_free,
      &dms_print, 0, 0, 0 };

  void
  Y_yoga_dms(int argc);
  void
  Y_yoga_addpzt(int argc);
  void
  Y_yoga_rmpzt(int argc);
  void
  Y_yoga_loadpzt(int argc);
  void
  Y_yoga_addkl(int argc);
  void
  Y_yoga_loadkl(int argc);
  void
  Y_yoga_rmkl(int argc);
  void
  Y_yoga_getkl(int argc);
  void
  Y_yoga_addtt(int argc);
  void
  Y_yoga_loadtt(int argc);
  void
  Y_yoga_setcomm(int argc);
  void
  Y_yoga_shapedm(int argc);
  void
  Y_yoga_resetdm(int argc);
  void
  Y_yoga_oneactu(int argc);
  void
  Y_yoga_getdm(int argc);
  void
  Y_yoga_setdm(int argc);
  void
  Y_target_dmtrace(int argc);
  void
  Y_dms_getdata(int argc);

  /*
   ____ _____ ____ 
   |  _ \_   _/ ___|
   | |_) || || |    
   |  _ < | || |___ 
   |_| \_\|_| \____|

   */

  void
  rtc_free(void *obj);
  void
  rtc_print(void *obj);

  static y_userobj_t yRTC = { const_cast<char*>("yRTC Object"), &rtc_free,
      &rtc_print, 0, 0, 0 };

  void
  Y_yoga_rtc(int argc);
  void
  Y_rtc_addcentro(int argc);
  void
  Y_rtc_addcontrol(int argc);
  void
  Y_rtc_rmcontrol(int argc);
  void
  Y_rtc_setthresh(int argc);
  void
  Y_rtc_setnmax(int argc);
  void
  Y_rtc_docentroids(int argc);
  void
  Y_rtc_doimat(int argc);
  void
  Y_rtc_setgain(int argc);
  void
  Y_rtc_loadmgain(int argc);
  void
  Y_rtc_setdelay(int argc);
  void
  Y_rtc_getimat(int argc);
  void
  Y_rtc_getcentroids(int argc);
  void
  Y_rtc_getcmat(int argc);
  void
  Y_rtc_buildcmat(int argc);
  void
  Y_rtc_imatsvd(int argc);
  void
  Y_rtc_framedelay(int argc);
  void
  Y_rtc_docontrol(int argc);
  void
  Y_controller_getdata(int argc);
  void
  Y_controller_initcured(int argc);

  /*
   ____  _                       
   / ___|| | ___  _ __   ___  ___ 
   \___ \| |/ _ \| '_ \ / _ \/ __|
   ___) | | (_) | |_) |  __/\__ \
    |____/|_|\___/| .__/ \___||___/
   |_|              

   */

  void
  Y_slopes_geom(int argc);
  void
  Y_sensors_getslopes(int argc);
  void
  Y_sensors_initnmax(int argc);
  void
  Y_sensors_initweights(int argc);
  void
  Y_sensors_initbcube(int argc);
  void
  Y_sensors_initcorr(int argc);
  void
  Y_sensors_loadcorrfnct(int argc);
  void
  Y_sensors_loadweights(int argc);
  void
  Y_sensors_compslopes(int argc);
  void
  Y_sensors_getnmax(int argc);
  void
  Y_centroider_getdata(int argc);

  /*
   ____ _       _           _ 
   / ___| | ___ | |__   __ _| |
   | |  _| |/ _ \| '_ \ / _` | |
   | |_| | | (_) | |_) | (_| | |
   \____|_|\___/|_.__/ \__,_|_|

   */

  void
  Y_move_atmos(int argc);
  void
  Y_move_sky(int argc);
  void
  Y_sensors_trace(int argc);
  void
  Y_yoga_add_telemetry_obj(int argc);
  void
  Y_yoga_add_telemetry_stream(int argc);

  /*
   _                       _       _       
   __ _  ___   | |_ ___ _ __ ___  _ __ | | __ _| |_ ___ 
   / _` |/ _ \  | __/ _ \ '_ ` _ \| '_ \| |/ _` | __/ _ \
| (_| | (_) | | ||  __/ | | | | | |_) | | (_| | ||  __/
   \__,_|\___/   \__\___|_| |_| |_| .__/|_|\__,_|\__\___|
   |_|                    

   */

  void
  aotemplate_free(void *obj);
  void
  aotemplate_print(void *obj);

  static y_userobj_t yAotemplate = { const_cast<char*>("yAotemplate Object"),
      &aotemplate_free, &aotemplate_print, 0, 0, 0 };

  void
  Y_yoga_aotemplate(int argc);
  void
  Y_yoga_templatefill(int argc);
  void
  Y_yoga_templatecomp(int argc);
  void
  Y_yoga_gettemplate(int argc);

}
#endif // _YOGA_AO_API_H_
