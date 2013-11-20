yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

require,yoga_ao_top+"/yorick/yoga_ao.i";
//require,yoga_ao_top+"/ywidgets/widget_wfs.i";

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data/";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

YOGA_AO_PARPATH = YOGA_AO_SAVEPATH+"par/";
mkdirp,YOGA_AO_PARPATH;

filename = YOGA_AO_PARPATH+"1wfs8x8_1layer_rtc_dm.par";

read_parfile,filename;

wfs_init;
sensors_compimg_tele, g_wfs, 0;

aaa = array(1.f, 64, 64);
aaa(*) = indgen(64*64);

g_as = yoga_acquisim(g_wfs,0);

window, 0;
acquisim_fillbcube, g_as, aaa;
sensors_compimg_tele, g_wfs, 0;
pli, sensors_getimg(g_wfs, 0);

window, 1;
acquisim_fillbcube, g_as, dist(64);
sensors_compimg_tele, g_wfs, 0;
pli, sensors_getimg(g_wfs, 0);

