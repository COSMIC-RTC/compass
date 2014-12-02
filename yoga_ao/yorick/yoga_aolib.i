plug_in,"yoga_ao";

require,"yoga.i";

/*
    _    ___              _                 _                 
   / \  / _ \    _____  _| |_ ___ _ __  ___(_) ___  _ __    _ 
  / _ \| | | |  / _ \ \/ / __/ _ \ '_ \/ __| |/ _ \| '_ \  (_)
 / ___ \ |_| | |  __/>  <| ||  __/ | | \__ \ | (_) | | | |  _ 
/_/   \_\___/   \___/_/\_\\__\___|_| |_|___/_|\___/|_| |_| (_)
                                                              
                                        
 _   _  ___   __ _  __ _     __ _  ___  
| | | |/ _ \ / _` |/ _` |   / _` |/ _ \ 
| |_| | (_) | (_| | (_| |  | (_| | (_) |
 \__, |\___/ \__, |\__,_|___\__,_|\___/ 
 |___/       |___/     |_____|          

*/

// atmosphere model
extern yoga_atmos;
/* DOCUMENT yoga_atmos
   obj = yoga_atmos(nscreens,r0,size,size2,alt,wspeed,wdir,deltax,deltay,pupil[,ndevice])

   creates an yAtmos object on the gpu
   nscreens : # of layers
   r0       : r0 for each layer
   size     : linear size of the screen
   size2    : second dim of extrude matrix A : size x size2
   alt      : array of altitudes per layers
   wspeed   : array of wind speeds per layers
   wdir     : array of wind directions per layers
   deltax   : array of x displacement per iteration (one per layers)
   deltay   : array of y displacement per iteration (one per layers)
   pupil    : array containing the pupil
   
   SEE ALSO:
 */
extern init_tscreen;
/* DOCUMENT init_tscreen
   init_tscreen,yoga_atmos_obj,altitude,a,b,istencilx,istencily,seed
     
   loads on the gpu in an yAtmos object and for a given screen data needed for extrude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   a              : matrix A
   b              : matrix B
   istencilx      : stencil for x direction
   istencily      : stencil for y direction
   seed           : seed for random numbers
   
   SEE ALSO:
 */
extern get_tscreen;
/* DOCUMENT get_tscreen
   screen = get_tscreen(yoga_atmos_obj,altitude)
     
   returns the screen in an yAtmos object and for a given altitude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   
   SEE ALSO:
 */
extern get_spupil;
/* DOCUMENT get_spupil
   pup = get_spupil(yoga_atmos_obj)
     
   returns the pupil in an yAtmos object
   yoga_atmos_obj : the yAtmos object
   
   SEE ALSO:
 */
extern get_tscreen_config;
/* DOCUMENT get_tscreen_config
   arr = get_tscreen_config(yoga_atmos_obj,altitude,str)
     
   returns config data in an yAtmos object and for a given altitude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   str            : type of data "A" or "B" or "istx" or "isty" or "values"
   
   SEE ALSO:
 */
extern get_tscreen_update;
/* DOCUMENT get_tscreen_update
   vect = get_tscreen_update(yoga_atmos_obj,altitude)
     
   returns only the update vector in an yAtmos object and for a given altitude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   
   SEE ALSO:
 */
extern extrude_tscreen;
/* DOCUMENT extrude_tscreen
   extrude_tscreen,yoga_atmos_obj,altitude[,dir]
     
   executes one col / row screen extrusion for a given altitude in an yAtmos object 
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   dir            : 0 (default) = column / 1 = row
   
   SEE ALSO:
 */

extern tscreen_initvk;
/* DOCUMENT tscreen_initvk
   tscreen_initvk,yoga_atmos_obj,altitude,pupd[,seed]
     
   initialize fourier screen computation with von karman statistics for a given altitude in an yAtmos object 
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   pupd           : the pupil diameter in pixels
   
   SEE ALSO:
 */

extern tscreen_genevk;
/* DOCUMENT tscreen_genevk
   tscreen_genevk,yoga_atmos_obj,altitude[,l0,nalias]
     
   generate screen with von karman statistics using fourier for a given altitude in an yAtmos object 
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   l0             : external scale (default inf)
   nalias         : number of pass for alias computation (default 0)
   
   SEE ALSO:
 */

// targets
extern yoga_target;
/* DOCUMENT yoga_target
   obj = yoga_target(ntargets,xpos,ypos,lambda,mag,sizes,pupil[,ndevice])
     
   creates an yTarget object on the gpu
   ntargets : # of targets
   xpos     : array of targets x positions
   ypos     : array of targets y positions
   lambda   : array of wavelengths
   mag      : array of magnitudes
   sizes    : array of linear # of pixels in the pupil
   pupil    : the pupil for image computation
   
   SEE ALSO:
 */
extern target_addlayer;
/* DOCUMENT target_addlayer
   target_addlayer,yoga_target_obj,ntarget,type,alt,xoff,yoff
   with type = "wfs" or "img"
     
   adds a layer of disturbances for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target
   type            : type of layer "atmos" or "optics"
   alt             : altitude of layer
   xoff            : array of x reference for raytracing
   yoff            : array of y reference for raytracing
   
   SEE ALSO:
 */
extern target_atmostrace;
/* DOCUMENT target_atmostrace
   target_atmostrace,yoga_target_obj,ntarget
     
   does raytracing through the atmos layers for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target
   
   SEE ALSO:
 */
extern target_dmtrace;
/* DOCUMENT target_dmtrace
   target_dmtrace,yoga_target_obj,ntarget,yoga_dm_obj
     
   does raytracing through the DMs layers for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target
   yoga_dm_obj     : the yDMs object
   
   SEE ALSO:
 */
extern target_init_strehlmeter;
/* DOCUMENT target_init_strehlmeter
   target_init_strehlmeter,yoga_target_obj,ntarget
     
   initialize structures for strehl meter
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target

   SEE ALSO:
 */
extern target_getstrehl;
/* DOCUMENT target_getstrehl
   strehl = target_getstrehl(yoga_target_obj,ntarget)
     
   initialize structures for strehl meter
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target

   SEE ALSO:
 */
extern target_getimage;
/* DOCUMENT target_getimage
   img = target_getimage(yoga_target_obj,ntarget,type)
     
   get given target image for a yTarget object and a yAtmos object
   yoga_target_obj : the yTarget object
   ntarget         : index of given target
   type            : the type of image ("se" : short exp., "le" : long exp.)

   SEE ALSO:
 */
 
extern target_getphase;
/* DOCUMENT target_getphase
   screen = target_getphase(yoga_target_obj,ntarget)
     
   returns the phase for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target

   SEE ALSO:
 */
extern target_getphasetele;
/* DOCUMENT target_getphasetele
   screen = target_getphasetele(yoga_target_obj,ntarget)
     
   returns the phase telemetry for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target

   SEE ALSO:
 */
extern target_getamplipup;
/* DOCUMENT target_getamplipup
   ampli = target_getamplipup(yoga_target_obj,ntarget)
     
   returns the complex amplitude in the pupil plane for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target

   SEE ALSO:
 */

// phase
extern yoga_phase;
//extern phase_copy;
extern phase_set;

// wfs
extern yoga_wfs;
/* DOCUMENT yoga_wfs
   obj = yoga_wfs(nxsub,nvalid,npix,nphase,nrebin,nfft,ntot,npup,pdiam,nphot,lgs[,ndevice])
     
   creates an yWfs object on the gpu
   nxsub  : linear # of subaps
   nvalid : number of valid subaps
   npix   : linear number of cam pixels in a subap
   nphase : linear size of phase screens in a subap
   nrebin : rebin factor from hr img to cam pixels
   nfft   : linear size of fft arrays for psf computations for a subap
   ntot   : linear size of hr image (total FoV) for a subap
   npup   : linear size of total pupil image
   pdiam  : linear diameter of a subap (m)
   nphot  : number of photons for a subap
   lgs    : flag for lgs mode : 1 if lgs wfs
      
   SEE ALSO:
 */
extern wfs_initgs;
/* DOCUMENT wfs_initgs
   wfs_initgs,yoga_wfs_obj,xpos,ypos,lambda,mag,size
     
   inits the guide star for an yWfs object
   yoga_wfs_obj : the yWfs object
   xpos         : x position of the guide star
   ypos         : y position of the guide star
   lambda       : observing wavelength for the guide star
   mag          : magnitude of the guide star
   size         : linear size of total image
   
   SEE ALSO:
 */

// acquisim
extern yoga_acquisim;
/* DOCUMENT yoga_acquisim
   obj = yoga_acquisim(yoga_sensors_obj)
     
   creates an yAcquisim object on the gpu
   yoga_sensors_obj : the ySensors object
      
   SEE ALSO:
 */
extern acquisim_fillbcube;
/* DOCUMENT acquisim_fillbcube
   acquisim_fillbcube,yoga_sensor_obj, sensor_number,image
     
   inits the guide star for an yWfs object
   yoga_sensors_obj : the ySensors object
   sensor_number: which sensor fill
   image        : image used to fill WFS bincube
   
   SEE ALSO:
 */

extern yoga_sensors;
/* DOCUMENT yoga_sensors
   obj = yoga_sensors(nsensors,nxsub,nvalid,npix,nphase,nrebin,nfft,ntot,npup,pdiam,nphot,lgs[,ndevice])
     
   creates an ySensors object on the gpu
   nsensors : # of wfs
   nxsub    : array of linear # of subaps
   nvalid   : array of numbers of valid subaps
   npix     : array of linear numbers of cam pixels in a subap
   nphase   : array of linear sizes of phase screens in a subap
   nrebin   : array of rebin factors from hr img to cam pixels
   nfft     : array of linear sizes of fft arrays for psf computations for a subap
   ntot     : array of linear sizes of hr image (total FoV) for a subap
   npup     : array of linear sizes of total pupil image
   pdiam    : linear diameter of a subap (m)
   nphot    : number of photons for a subap
   lgs      : array of flags for lgs mode : 1 if lgs wfs
      
   SEE ALSO:
 */
extern sensors_initgs;
/* DOCUMENT sensors_initgs
   sensors_initgs,yoga_sensors_obj,xpos,ypos,lambda,mag,size
     
   inits the guide stars for an ySensors object
   yoga_sensors_obj : the ySensors object
   xpos             : array of x positions of the guide stars
   ypos             : array of y positions of the guide stars
   lambda           : array of observing wavelengths for the guide stars
   mag              : array of magnitudes of the guide stars
   size             : array of linear sizes of total images
   
   SEE ALSO:
 */
extern sensors_rmlayer;
extern sensors_addlayer;
/* DOCUMENT sensors_addlayer
   sensors_addlayer,yoga_sensors_obj,nsensor,type,alt,xoff,yoff
   type not used can be "atmos"
     
   add a disturbance layer for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   type             : type of layer "atmos" or "optics"
   alt              : altitude of layer
   xoff            : array of x reference for raytracing
   yoff            : array of y reference for raytracing
  
   SEE ALSO:
 */
extern sensors_initarr;
/* DOCUMENT sensors_initarr
   sensors_initarr,yoga_sensors_obj,nsensor,phasemap,hrmap,binmap,offsets,pupil,fluxpersub,isvalid,
   validsubsx,validsubsy,istart,jstart
     
   init arrays for image computation for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   phasemap         : array of pixels transform from phase screen into
                      subaps phase screens
   hrmap            : array of pixels transform from minimal FoV image to
                      full FoV image (array of 0 if the same)
   binmap           : array of pixels transform from full FoV hr images to
                      binned images
   offsets          : array of pixels offsets for subaps phase screens
   pupil            : the pupil array
   fluxpersub       : array of flux per subaperture
   isvalid          : array nxsub x nxsub of 0/1 for valid subaps
   validsubsx       : array nvalid x coordinates of valid subs in a nxsub x nxsub array
   validsubsy       : array nvalid y coordinates of valid subs in a nxsub x nxsub array
   istart           : i index of first phase elem of subap
   istart           : j index of first phase elem of subap
   
   SEE ALSO:
 */
extern sensors_compimg;
/* DOCUMENT sensors_compimg
   sensors_compimg,yoga_sensors_obj,nsensor
     
   image computation for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs

   SEE ALSO:
 */

extern sensors_compimg_tele;
/* DOCUMENT sensors_compimg_tele
   sensors_compimg_tele,yoga_sensors_obj,nsensor
     
   image computation for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs

   SEE ALSO:
 */

extern slopes_geom;
/* DOCUMENT slopes_geom
   slopes_geom,yoga_sensors_obj,nsensor,lambOverD
     
   geometric slopes computation for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   lambOverD        : lambda over d for one subap

   SEE ALSO:
 */
extern sensors_getimg;
/* DOCUMENT sensors_getimg
   img = sensors_getimg(yoga_sensors_obj,nsensor)
     
   returns the total image for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs

   SEE ALSO:
 */
extern sensors_getslopes;
/* DOCUMENT sensors_getslopes
   slopes = sensors_getslopes(yoga_sensors_obj,nsensor)
     
   returns the slopes for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs

   SEE ALSO:
 */
extern sensors_getdata;
/* DOCUMENT sensors_getdata
   arr = sensors_getdata(yoga_sensors_obj,nsensor,type)
   type = 
     
   returns the configuration data for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   type             : type of data "amplipup" or "amplifoc" or "hrimg" or
                      "bincube" or "phase" or "totimg" or "lgskern" or "ftlgskern"

   SEE ALSO:
 */

extern sensors_setphase;
/* DOCUMENT sensors_setphase
   arr = sensors_setphase(yoga_sensors_obj,nsensor,phase)
     
   sets phase of a given sensor in a ySensors object to the array phase
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   phase            : the new phase array

   SEE ALSO:
 */

extern sensors_loadkernels;
/* DOCUMENT sensors_loadkernels
   sensors_loadkernels,yoga_sensors_obj,nsensor,kernels
     
   load lgs kernels for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   kernels          : array of lgs kernels for each subap

   SEE ALSO:
 */

extern sensors_initlgs;
/* DOCUMENT sensors_initlgs
   sensors_initlgs,yoga_sensors_obj,nsensor,nprof,hg,h0,deltah,pixsize,doffaxis,prof1d,profcum,azimuth
     
   load all arrays needed for lgs spot computation for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   nprof            : # of elements in profile
   hg               : altitude in m of profile cog
   h0               : altitude in m of first element in profile
   deltah           : profile resolution in m
   pixsize          : high resolution image pixel size (in arcsec)
   doffaxis         : array of off-axis distance
   prof1d           : the sodium layer profile
   profcum          : the sodium layer profile integrated
   azimuth          : angles of rotation for each spot (rad)

   SEE ALSO:
 */

extern sensors_updatelgsprof;
/* DOCUMENT sensors_updatelgsprof
   sensors_updatelgsprof,yoga_sensors_obj,nsensor,prof1d,profcum,ideriv,interpx,interpw
     
   update arrays for lgs spot computation for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   prof1d           : the sodium layer profile
   profcum          : the sodium layer profile integrated
   ideriv           : array of 0/1 for subaps for which profile has to be integrated before interp
   interpx          : map of coordinates for profile interpolation
   interpw          : map of weights for profile interpolation

   SEE ALSO:
 */

extern sensors_updatelgs;
/* DOCUMENT sensors_updatelgs
   sensors_updatelgs,yoga_sensors_obj,nsensor
     
   update lgs spot for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs

   SEE ALSO:
 */


extern yoga_rtc;
/* DOCUMENT yoga_rtc
   obj = yoga_rtc(device)
   creates an yRTC object on the gpu
      
   SEE ALSO:
 */

extern rtc_addcentro;
/* DOCUMENT rtc_addcentro
   rtc_addcentro,yoga_rtc_obj,nvalid,type
   add centroider to yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_addcontrol;
/* DOCUMENT rtc_addcontrol
   rtc_addcontrol,yoga_rtc_obj,nactu,delay,type
   add controller to yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_rmcontrol;
/* DOCUMENT rtc_rmcontrol
   rtc_rmcontrol,yoga_rtc_obj
   rm last controller of yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_setthresh;
/* DOCUMENT rtc_setthresh
   rtc_setthresh,yoga_rtc_obj,ncentro,thresh
   set threshold for a centroider in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_setnmax;
/* DOCUMENT rtc_setnmax
   rtc_setnmax,yoga_rtc_obj,ncentro,nmax
   set nmax for a centroider in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_docentroids;
/* DOCUMENT rtc_docentroids
   rtc_docentroids,yoga_rtc_obj,ncontroller,yoga_wfs_obj
   do centroiding for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_docentroids_geom;
/* DOCUMENT rtc_docentroids_geom
   rtc_docentroids_geom,yoga_rtc_obj,ncontroller,yoga_wfs_obj
   do centroiding using non diffractive wfs model
   for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_doCmm;
/* DOCUMENT rtc_doCmm
   rtc_doCmm,yoga_rtc_obj,ncontroller,yoga_wfs_obj,yoga_atmos_obj,diamTel,cobs
   do Cmm matrix for a mv controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_doCphim;
/* DOCUMENT rtc_doCmm
   rtc_doCphim,yoga_rtc_obj,ncontroller,yoga_wfs_obj,yoga_atmos_obj,yoga_dm_obj,L0,cn2,alphaX
                alphaY,xactu,yactu,diamTel,k2;
   do Cphim matrix for a mv controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_doimatkl;
extern rtc_doimat;
/* DOCUMENT rtc_doimat
   rtc_doimat,yoga_rtc_obj,ncontroller,yoga_wfs_obj,yoga_dm_obj[,geomtype]
   do imat for a controller in a yoga_rtc_obj object
   geomtype : optional, type of geometric slopes computation
      
   SEE ALSO:
 */
extern rtc_getnoisemat;
extern rtc_getimat;
/* DOCUMENT rtc_getimat
   res = rtc_getimat(yoga_rtc_obj,ncontroller)
   get imat for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_setimat;
/* DOCUMENT rtc_setimat
   rtc_setimat,yoga_rtc_obj,ncontroller,imat
   set imat for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_setcovmat;
/* DOCUMENT rtc_covmat
   rtc_setcovmat,yoga_rtc_obj,ncontroller,covmat
   set covmat (Cphi) for a MV controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_setcentroids;
/* DOCUMENT rtc_setcentroids
   rtc_setcentroids,yoga_rtc_obj,ncontroller,imat
   set centroids for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_getcentroids;
/* DOCUMENT rtc_getcentroids
   res = rtc_getcentroids(yoga_rtc_obj,ncontroller)
   get centroids for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_getcom;
/* DOCUMENT rtc_getcom
   res = rtc_getcom(yoga_rtc_obj,ncontroller)
   get commands for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_getcmat;
/* DOCUMENT rtc_getcmat
   res = rtc_getcmat(yoga_rtc_obj,ncontroller)
   get cmat for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_setcmat;
/* DOCUMENT rtc_setcmat
   rtc_setcmat,yoga_rtc_obj,ncontroller,cmat
   set cmat for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_buildcmatmv;
extern rtc_buildcmat;
/* DOCUMENT rtc_buildcmat
   rtc_buildcmat,yoga_rtc_obj,ncontroller,nfilt,filt_tt
   do generalized inverse of imat for a controller in a yoga_rtc_obj object
   nfilt   : number of modes filtered (lowest eigenvalues)
   filt_tt : 0 (default) / 1 : flag to keep or not tt in the inversion
   SEE ALSO:
 */

extern rtc_setgain;
/* DOCUMENT rtc_setgain
   rtc_setgain,yoga_rtc_obj,ncontroller,gain
   set loop gain for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_docovmat;
/* DOCUMENT rtc_docovmat
   rtc_docovmat,yoga_rtc_obj,ncontroller,yoga_dm_obj,type,alt,indx_pup,dim,xpos,ypos,norm,ampli,method
   compute the covariance matrix on a yoga_rtc_obj with MV controller
*/
extern rtc_loadcovmat;
/* DOCUMENT rtc_loadcovmat
   rtc_loadcovmat,yoga_rtc_obj,ncontroller,covmat
   set the inverse of the covariance matrix on a yoga_rtc_obj with MV controller

   SEE ALSO:
*/
extern rtc_loadklbasis;
/* DOCUMENT rtc_loadklbasis
   rtc_loadklbasis,yoga_rtc_obj,ncontroller,klbasis
   set the KL basis on a yoga_rtc_obj with MV controller

   SEE ALSO:
*/
extern rtc_getklbasis;
/* DOCUMENT rtc_getklbasis
   rtc_getklbasis,yoga_rtc_obj,ncontroller
   retrieve the KL basis from a yoga_rtc_obj with MV controller

   SEE ALSO:
*/
extern rtc_getCmm;
/* DOCUMENT rtc_getcovmat
   rtc_getCmm,yoga_rtc_obj,ncontroller
   retrieve the inverse of the covariance matrix  from a yoga_rtc_obj with MV controller

   SEE ALSO:
*/
extern rtc_getCphim;
/* DOCUMENT rtc_getcovmat
   rtc_getCmm,yoga_rtc_obj,ncontroller
   retrieve the inverse of the covariance matrix  from a yoga_rtc_obj with MV controller

   SEE ALSO:
*/
extern rtc_getcovmat;
/* DOCUMENT rtc_getcovmat
   rtc_getcovmat,yoga_rtc_obj,ncontroller
   retrieve the inverse of the covariance matrix  from a yoga_rtc_obj with MV controller

   SEE ALSO:
*/
extern rtc_loadmgain;
/* DOCUMENT rtc_loadmgain
   rtc_loadmgain,yoga_rtc_obj,ncontroller,mgain
   set modal gain for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_loadnoisemat;
/* DOCUMENT rtc_loadnoisemat
   rtc_loadnoisemat,yoga_rtc_obj,ncontroller,noise
   set noise matrix for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */
extern rtc_setdelay;
/* DOCUMENT rtc_setdelay
   rtc_setdelay,yoga_rtc_obj,ncontroller,delay
   set loop delay for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_imatsvd;
/* DOCUMENT rtc_imatsvd
   rtc_imatsvd,yoga_rtc_obj,ncontroller
   compute svdec of imat for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_framedelay;
/* DOCUMENT rtc_framedelay
   rtc_framedelay,yoga_rtc_obj,ncontroller
   introduce frame delay a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_docontrol;
/* DOCUMENT rtc_docontrol
   rtc_docontrol,yoga_rtc_obj,ncontrolleryoga_dm_obj
   compute and apply command (including delay) for a controller in a yoga_rtc_obj object
      
   SEE ALSO:
 */

extern rtc_kalmancalculategain;
/* DOCUMENT rtc_kalmancalculategain
   rtc_kalmancalculategain,yoga_rtc_obj,ncontroller,bruit, SigmaV, atur, btur
      
   SEE ALSO:
 */

extern rtc_initkalman;
/* DOCUMENT rtc_initkalman
   rtc_initkalman,yoga_rtc_obj,ncontroller,bruit, D_Mo, N_Act, PROJ, SigmaV, atur, btur, is_zonal, is_sparse, is_GPU
      
   SEE ALSO:
 */

extern rtc_kalmangettime;
/* DOCUMENT rtc_kalmangettime
   rtc_initkalman,yoga_rtc_obj,ncontroller
      
   SEE ALSO:
 */

extern yoga_dms;
/* DOCUMENT yoga_dms
   obj = yoga_dms(ndm)
   creates an yDMs object on the gpu containing ndm DMs
      
   SEE ALSO:
 */

extern yoga_addpzt;
/* DOCUMENT yoga_addpzt
   yoga_addpzt,yoga_dms_obj,alt,dim,ninflu,influsize,ninflupos,n_npts,push4imat[,device]
   add a pzt dm to yoga_dms_obj object
      
   alt         : conjugated altitude
   dim         : linear size in pixels of the dm support
   ninflu      : number of actuators
   influsize   : linear size of the influence functions support
   ninflupos   : number of influ index  for shape computation
   n_npts      : number of influ points for shape computation
   push4imat  : nominal voltage for imat
   device      : optional : device on which to create the dm
   SEE ALSO:
 */

extern yoga_rmpzt;
/* DOCUMENT yoga_rmpzt
   yoga_rmpzt,yoga_dms_obj,alt
   remove a pzt dm from yoga_dms_obj object
      
   alt         : conjugated altitude
   SEE ALSO:
 */


extern yoga_loadpzt;
/* DOCUMENT yoga_loadpzt
   yoga_loadpzt,yoga_dms_obj,alt,influ,influpos,npoints,istart,xoff,yoff
   load data for a pzt dm from a yoga_dms_obj object
      
   alt        : conjugated altitude
   influ      : influence functions (if)
   influpos   : position of each point of the if contributing to the dm shape
   npoints    : number of if contributing to each point of the dm shape
   istart     : start index for dm shape computation
   xoff       : x offset for each if
   yoff       : y offset for each if
   SEE ALSO:
 */


extern yoga_addkl;
/* DOCUMENT yoga_addkl
   yoga_addkl,yoga_dms_obj,alt,dim,ninflu,influsize,nr,np,push4imat[,device]
   add a pzt dm to yoga_dms_obj object
      
   alt         : conjugated altitude
   dim         : linear size in pixels of the dm support
   ninflu      : number of actuators
   influsize   : linear size of the influence functions support
   ninflupos   : number of influ index  for shape computation
   nr          : number of radial points
   np          : number of azimuth points
   device      : optional : device on which to create the dm
   SEE ALSO:
 */

extern yoga_rmkl;
/* DOCUMENT yoga_rmkl
   yoga_rmkl,yoga_dms_obj,alt
   remove a kl dm from yoga_dms_obj object
      
   alt         : conjugated altitude
   SEE ALSO:
 */
extern yoga_getflokl;
extern yoga_getkl;
/* DOCUMENT yoga_getkl
   yoga_getkl,yoga_dms_obj,alt,nkl
   
   returns kl #nkl
      
   alt         : conjugated altitude
   nkl         : 
   SEE ALSO:
 */

extern yoga_floloadkl;
extern yoga_loadkl;
/* DOCUMENT yoga_loadkl
   yoga_loadkl,yoga_dms_obj,alt,rasbas,azbas,ord,cr,cp,xoff,yoff
   load data for a kl dm from a yoga_dms_obj object
      
   alt        : conjugated altitude
   rasbas     : radial basis of kl
   azbas      : azimuth basis of kl
   order      : radial order of kl
   cr         : cartesian for radial
   cp         : cartesian for azimuth
   xoff       : x offset for each if
   yoff       : y offset for each if
   SEE ALSO:
 */


extern yoga_addtt;
/* DOCUMENT yoga_addtt
   yoga_addtt,yoga_dms_obj,alt,dim,push4imat[,device]
   add a dm to yoga_dms_obj object
      
   alt         : conjugated altitude
   dim         : linear size in pixels of the dm support
   push4imat   : nominal voltage for imat
   device      : optional : device on which to create the dm
   SEE ALSO:
 */


extern yoga_loadtt;
/* DOCUMENT yoga_loadtt
   yoga_loadtt,yoga_dms_obj,alt,influ
   load data for a tt dm from a yoga_dms_obj object
      
   alt        : conjugated altitude
   influ      : influence functions (if)
   SEE ALSO:
 */


extern yoga_setcomm;
/* DOCUMENT yoga_setcomm
   yoga_setcomm,yoga_dms_obj,type,alt,comm
   set the command vector of a dm in a yoga_dms_obj object
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   comm        : command vector
   SEE ALSO:
 */

extern yoga_getcomm;
/* DOCUMENT yoga_getcomm
   comm = yoga_getcomm(yoga_dms_obj,type,alt)
   get the command vector of a dm in a yoga_dms_obj object
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   comm        : command vector
   SEE ALSO:
 */

extern yoga_getdm;
/* DOCUMENT yoga_getdm
   ima = yoga_getdm(yoga_dms_obj,type,alt)
   retreive the current shape of a dm in a yoga_dms_obj object
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   SEE ALSO:
 */

extern yoga_setdm;
/* DOCUMENT yoga_setdm
   yoga_setdm,yoga_dms_obj,type,alt,shape
   set the shape of a dm in a yoga_dms_obj object
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   shape       : dm shape
   SEE ALSO:
 */

extern yoga_shapedm;
/* DOCUMENT yoga_shapedm
   yoga_shapedm,yoga_dms_obj,type,alt
   computes the shape of a dm in a yoga_dms_obj object
   using the current command vector
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   SEE ALSO:
 */

extern yoga_resetdm;
/* DOCUMENT yoga_resetdm
   yoga_resetdm,yoga_dms_obj,type,alt
   reset the shape of a dm in a yoga_dms_obj object
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   SEE ALSO:
 */

extern yoga_oneactu;
/* DOCUMENT yoga_shapedm
   yoga_shapedm,yoga_dms_obj,type,alt
   computes the shape of a dm in a yoga_dms_obj object
   using the current command vector
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   nactu       : the actuator to be moved
   ampli       : the amplitude
   SEE ALSO:
 */

extern dms_getdata;
/* DOCUMENT dms_getdata
   dms_getdata,yoga_dms_obj,type,alt,typedata
   get data from a yDM
      
   type        : type of dm : pzt, tt
   alt         : conjugated altitude
   typedata    : type of data
   SEE ALSO:
 */

extern dms_comp_shape;
/* DOCUMENT dms_getdata
   dms_getdata,yoga_dms_obj,com
   
   com : command vector computed with sutra_controller
   SEE ALSO:
 */

extern sensors_compslopes;
/* DOCUMENT sensors_compslopes
   sensors_compslopes,yoga_sensors_obj,nsensor,yoga_rtc_obj,ncentro[,nmax/threshold]
     
   general slopes computation for a given sensor in a ySensors object
   
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   yoga_rtc_obj     : the yRTC object
   ncentro          : index of given centroider
   nmax             : nmax brightest pixels (if not provided set to default in centroider struct)
   threshold        : threshold value (if not provided set to default in centroider struct)

   SEE ALSO:
 */

extern sensors_initnmax;
/* DOCUMENT sensors_initnmax
   sensors_initnmax,yoga_sensors_obj,nsensor,nmax
     
   init structures for nmax brightest pixels centroiding for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   nmax             : nmax brightest pixels

   SEE ALSO:
 */

extern sensors_initweights;
/* DOCUMENT sensors_initweights
   sensors_initweights,yoga_sensors_obj,nsensor,yoga_rtc_obj,ncentro,weights
     
   init structures for weighted cog centroiding for a given sensor in a ySensors object
   and loads corresponding weighting functions
   
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   yoga_rtc_obj     : the yRTC object
   ncentro          : index of given centroider
   weights          : weighting function can be 2d (same for all subaps) or 3d array

   SEE ALSO:
 */

extern sensors_initbcube;
/* DOCUMENT sensors_initbcube
   sensors_initbcube,yoga_sensors_obj,nsensor,yoga_rtc_obj,ncentro
     
   init structures for bincube in centroider  for a given sensor in a ySensors object
   
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   yoga_rtc_obj     : the yRTC object
   ncentro          : index of given centroider

   SEE ALSO:
 */

extern sensors_loadweights;
/* DOCUMENT sensors_loadweights
   sensors_loadweights,yoga_sensors_obj,yoga_rtc_obj,ncentro,weights
     
   loads weighting functions for weighted cog centroiding for a given sensor in a ySensors object
   
   yoga_sensors_obj : the ySensors object
   yoga_rtc_obj     : the yRTC object
   ncentro          : index of given centroider
   weights          : weighting function can be 2d (same for all subaps) or 3d array

   SEE ALSO:
 */

extern sensors_initcorr;
/* DOCUMENT sensors_initcorr
   sensors_initcorr,yoga_sensors_obj,nsensor,yoga_rtc_obj,ncentro,corrfnct,corr_norm,sizex,sizey,interpmat
     
   init structures for correlation centroiding for a given sensor in a ySensors object
   and loads corresponding corr functions
   
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   yoga_rtc_obj     : the yRTC object
   ncentro          : index of given centroider
   corrfnct         : function to use for correlation can be 2d (same for all subaps) or 3d array
   corr_norm        : norm factor for correlation
   sizex            : window x size for interpolation
   sizey            : window y size for interpolation
   interpmat        : matrix used for interpolation

   SEE ALSO:
 */

extern sensors_loadcorrfnct;
/* DOCUMENT sensors_loadcorrfnct
   sensors_loadcorrfnct,yoga_sensors_obj,yoga_rtc_obj,ncentro,corrfnct,corr_norm
     
   loads functions to use for correlationcentroiding for a given sensor
   in a ySensors object
   
   yoga_sensors_obj : the ySensors object
   yoga_rtc_obj     : the yRTC object
   ncentro          : index of given centroider
   corrfnct         : function to use for correlation can be 2d
                      (same for all subaps) or 3d array
   corr_norm        : norm factor for correlation

   SEE ALSO:
 */

extern sensors_getnmax;
/* DOCUMENT sensors_getnmax
   sensors_getnmax,yoga_sensors_obj,nsensor
     
   search for nmax brightest pixels per subap for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs

   SEE ALSO:
 */

extern centroider_getdata;
/* DOCUMENT centroider_getdata
   arr = centroider_getdata(yoga_rtc_obj,ncentro,type)
     
   returns the configuration data for a given centroider in a yRTC object
   
   yoga_rtc_obj : the yRTC object
   ncentro      : index of given centroider
   type         : type of data "weights" 

   SEE ALSO:
 */

extern controller_getdata;
/* DOCUMENT controller_getdata
   arr = controller_getdata(yoga_rtc_obj,ncontrol,type)
     
   returns the configuration data for a given controller in a yRTC object
   
   yoga_rtc_obj : the yRTC object
   ncontrol     : index of given controller
   type         : type of data "eigenvals" 

   SEE ALSO:
 */

extern controller_setdata;
/* DOCUMENT controller_setdata
   ontroler_setdata,yoga_rtc_obj,ncontrol,type,data
     
   returns the configuration data for a given controller in a yRTC object
   
   yoga_rtc_obj : the yRTC object
   ncontrol     : index of given controller
   type         : type of data "eigenvals", "mes2mod", "mod2act" 

   SEE ALSO:
 */

extern controller_initcured;
/* DOCUMENT controller_initcured
   controller_initcured,yoga_rtc_obj,ncontrol,nxsubs,isvalid,ndivs,tt
     
   init a cured controller with data
   
   yoga_rtc_obj : the yRTC object
   ncontrol     : index of given controller
   nxsubs       : number of subaps across pupil (int)
   isvalid      : array of valid subaps in a nxsubs x nxsubs configuration (int)
   ndivs        : number of subdivision levels for cured
   tt           : 1 if separate tt, O is no tt mirror

   SEE ALSO:
 */

// global
extern move_atmos;
/* DOCUMENT move_atmos
   move_atmos,yoga_atmos_obj
     
   multi-layer extrude process for a yAtmos object

   SEE ALSO:
 */
extern move_sky;
/* DOCUMENT move_sky
   move_sky,yoga_atmos_obj,yoga_target_obj
     
   multi-layer & multi-target extrude process for a yAtmos object and a yTarget object

   SEE ALSO:
 */
extern sensors_trace;
/* DOCUMENT sensors_trace
   sensors_trace,yoga_sensors_obj,nsensor,yoga_atmos_obj
   or
   sensors_trace,yoga_sensors_obj,nsensor,yoga_dms_obj[,rst]
   
   do raytracing for a given sensor in a ySensors object and given turbulence in an yAtmos object
   or a given yDMs object
   
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   type             : type of structure through which to trace : atmos / dm / ...
   yoga_atmos_obj   : the yAtmos object
   yoga_dms_obj     : the yDMs object
   rst              : optional (0) reset phase before dm computation

   SEE ALSO:
 */

extern yoga_telemetry;
/* DOCUMENT yoga_telemetry
   yTele = yoga_telemetry(dimsObj, nbStreams)
   creates an yoga_telemetry object containing :
      - one yoga_obj with dimsObj as dimensions
      - nbStreams streams
   SEE ALSO:
 */

 extern yoga_add_telemetry_obj;
/* DOCUMENT yoga_add_telemetry_obj
   yoga_add_telemetry_obj, yTele, name, obj

   SEE ALSO:
 */

 extern yoga_add_telemetry_stream;
/* DOCUMENT yoga_add_telemetry_stream
   yoga_add_telemetry_stream, yTele

   SEE ALSO:
 */

 extern yoga_aotemplate;
/* DOCUMENT yoga_aotemplate
   yTemplate = yoga_aotemplate(nbelem[,device])
   creates an yoga_aotemplate object
   SEE ALSO:
 */

 extern yoga_templatefill;
/* DOCUMENT yoga_templatefill
   yoga_templatefill, yTemplate[, data]
   fill a yTemplate with optional data (else filled with random numbers)

   SEE ALSO:
 */

 extern yoga_templatecomp;
/* DOCUMENT yoga_templatecomp
   yoga_templatecomp, yTemplate
   do computation on a yTemplate

   SEE ALSO:
 */

 extern yoga_gettemplate;
/* DOCUMENT yoga_gettemplate
   res = yoga_gettemplate(yTemplate[,type])
   get content of yTemplate (type optional can be void, "data" or "res")

   SEE ALSO:
 */
extern rtc_doimatkl4pzt;
