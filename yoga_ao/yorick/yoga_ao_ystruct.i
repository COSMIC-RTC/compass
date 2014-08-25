struct geom_struct
{
  long  ssize;       // linear size of full image (in pixels)
  float zenithangle; // observations zenith angle (in deg)
  
  // internal keywords
  long  pupdiam;     // linear size of total pupil (in pixels)
  float cent;        // central point of the simulation
  pointer _ipupil;   // total pupil (include full guard band)
  pointer _mpupil;   // medium pupil (part of the guard band)
  pointer _spupil;   // small pupil (without guard band)
  long  _p1;         // min x,y for valid points in mpupil
  long  _p2;         // max x,y for valid points in mpupil
  long  _n;          // linear size of mpupil
  long  _n1;         // max x,y for valid points in ipupil
  long  _n2;         // min x,y for valid points in ipupil
};

struct tel_struct
{
  float diam;        // telescope diameter (in meters)
  float cobs;        // central obstruction ratio
  float t_spiders;   // secondary supports ratio
  string spiders_type; // secondary supports type: "four" or "six"
};

struct atmos_struct
{
  long    nscreens;    // number of turbulent layers
  float   r0;          // global r0 @ 0.5µm
  float   pupixsize;   // pupil piwel size (in meters)
  pointer   L0;        // L0 per layers in meters
  pointer dim_screens; // linear size of phase screens
  pointer alt;         // altitudes of each layer
  pointer winddir;     // wind directions of each layer
  pointer windspeed;   // wind speeds of each layer
  pointer frac;        // fraction of r0 for each layer
  pointer deltax;      // x translation speed (in pix / iteration) for each layer
  pointer deltay;      // y translation speed (in pix / iteration) for each layer
  pointer seeds;       // 
};

struct target_struct
{
  long    ntargets;  // number of targets
  pointer lambda;    // observation wavelength for each target
  pointer xpos;      // x positions on sky (in arcsec) for each target
  pointer ypos;      // y positions on sky (in arcsec) for each target
  pointer mag;       // magnitude for each target
};

struct wfs_struct
{
  string type;           // type of wfs : "sh" or "pyr"
  long   nxsub;          // linear number of subaps
  long   npix;           // number of pixels per subap
  float  pixsize;        // pixel size (in arcsec) for a subap
  float  lambda;         // observation wavelength (in µm) for a subap
  float  optthroughput;  // wfs global throughput
  float  fracsub;        // minimal illumination fraction for valid subaps
  long   openloop;       // 1 if in "open-loop" mode (i.e. does not see dm)
  float  fssize;       //size of field stop in arcsec
  string fstop;       //size of field stop in arcsec
  
  //target kwrd
  float xpos;      // guide star x position on sky (in arcsec) 
  float ypos;      // guide star x position on sky (in arcsec) 
  float gsalt;     // altitude of guide star (in m) 0 if ngs 
  float gsmag;     // magnitude of guide star
  float zerop;     // detector zero point
  float noise;     // desired noise : < 0 = no noise / 0 = photon only / > 0 photon + ron 
  
  float kernel;      // 
  pointer _ftkernel; // 

  // lgs only
  float lgsreturnperwatt;  // return per watt factor (high season : 10 ph/cm2/s/W)
  float laserpower;        // laser power in W
  float lltx;              // x position (in meters) of llt
  float llty;              // y position (in meters) of llt
  string proftype;         // type of sodium profile "gauss", "exp", etc ...
  float beamsize;          // laser beam fwhm on-sky (in arcsec)

  //internal kwrd
  long  _pdiam;          // pupil diam for a subap (in pixel)
  long  _Nfft;           // array size for fft for a subap (in pixel)
  long  _Ntot;           // total size of hr image for a subap (in pixel)
  long  _nrebin;         // rebin factor from hr to binned image for a subap 
  long  _nvalid;         // number of valid subaps

  float   _nphotons;     // number of photons per subap
  float   _subapd;       // subap diameter (m)
  pointer _fluxPerSub;   // fraction of nphotons per subap
  float   _qpixsize;     // quantum pixel size for the simulation
  
  pointer _istart;       // x start indexes for cutting phase screens 
  pointer _jstart;       // y start indexes for cutting phase screens 
  pointer _validsubs;    // (i,j) indices of valid subaps
  pointer _isvalid;      // array of 0/1 for valid subaps
  pointer _phasemap;     // array of pixels transform from phase screen into
                         // subaps phase screens
  pointer _hrmap;        // array of pixels transform from minimal FoV image to
                         // full FoV image (array of 0 if the same)
  pointer _binmap;       // array of pixels transform from full FoV hr images to
                         // binned images
  pointer _halfxy;       // phase offset for 1/2 pixel shift in (x,y)
  
  pointer _submask;       // fieldstop for each subap
  
  pointer _lgskern;      // lgs kernels for each subap
  pointer _profna;       // sodium profile
  pointer _altna;        // corresponding altitude
  pointer _prof1d;       // hr profile
  pointer _profcum;      // hr profile cumulated
  pointer _beam;         // 1d beam function
  pointer _ftbeam;       // 1d beam function fft
  pointer _azimuth;      // angles of rotation for each spot

  // pyramid-nly kwrds
  float   pyr_ampl;   // pyramid wfs modulation amplitude radius [arcsec]
  long    pyr_npts;   // total number of point along modulation circle [unitless]
  pointer pyr_pos;    // positions for modulation, overwrites ampl and npts [arcsec]
  string  pyr_loc;    // Location of modulation, before/after the field stop.
                          // valid value are "before" or "after" (default "after")
  string  pyrtype;    // Type of pyramid, either 0 for "Pyramid" or 1 for "RoofPrism".

  // pyramid internal kwrds
  pointer _pyr_offsets;   //
  pointer _pyr_cx;   //
  pointer _pyr_cy;   //
  
};

struct dm_struct
{
  string  type;            // type of dm
  long    nact;            // number of actuators in the diameter 
  float   alt;             // conjugaison altitude (im m)
  float   thresh;          // threshold on response for selection (<1)
  float   coupling;        // actuators coupling (<0.3)
  float   hyst;            // actuators hysteresis (<1.)
  float   margin;          //
  float   pupoffset(2);    // global offset in pupil of whole actuator pattern [m]
  float   unitpervolt;     // Influence function sensitivity in unit/volt. Optional [0.01]
                           // Stackarray: mic/volt, Tip-tilt: arcsec/volt.
  float   push4imat;       // nominal voltage for imat
  
  long    nkl;             // number of kl modes
  
  //internal kwrd
  long    _pitch;     // inter-actuator space in pixels
  long    _ntotact;   // total number of actuators
  long    _influsize; // total number of actuators
  long    _n1;        // position of leftmost pixel in largest support
  long    _n2;        // position of rightmost pixel in largest support
  long    _puppixoffset(2);
  pointer _influ;     // influence functions
  pointer _xpos;      // x positions of influ functions
  pointer _ypos;      // y positions of influ functions
  pointer _i1;            //
  pointer _j1;            //
  pointer _pupil;     // pupil mask for this dm
  pointer _com;       // current command
  pointer _influpos;  //
  pointer _ninflu;    // 
  pointer _influstart;//
  
  pointer _klbas;           // pointer to a kl struct
};

struct rtc_struct
{
  long    nwfs;      // number of wfs
  pointer centroiders; // an array of centroiders
  pointer controllers; // an array of controllers
};

struct centroider_struct
{
  long    nwfs;      // index of wfs in y_wfs structure on which we want to do centroiding
  string  type;      // type of centroiding "cog", "tcog", "bpcog", "wcog", "corr"
  string  type_fct;  // type of ref function "gauss", "file", "model"
  pointer weights;   // optional reference function(s) used for centroiding
  long    nmax;      //
  float   thresh;    //
  float   width;     // width of the Gaussian
  long    sizex;     //
  long    sizey;     //
  pointer interpmat; // optional reference function(s) used for centroiding
};

struct controller_struct
{
  string  type;     // type of controller
  pointer nwfs;     // index of wfss in controller
  pointer nvalid;   // number of valid subaps per wfs
  pointer ndm;      // index of dms in controller
  pointer nactu;    // number of controled actuator per dm
  pointer imat;     // full interaction matrix
  pointer cmat;     // full control matrix
  float   maxcond;  // max condition number
  long    delay;    // 
  float   gain;     //
  long    nkl;      // Florain features : number of KL modes used for computation of covmat in case of minimum variance controller
  long    cured_ndivs; // subdivision levels in cured
};

struct loop_struct
{
  long  niter;      // number of iterations
  float ittime;     // iteration time (in sec)
};

struct telemetry_struct
{
  pointer nwfs;     // wfs tab
};
