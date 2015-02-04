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
//sensors_compimg, g_wfs, 0;

aaa = array(1.f, 64, 64);
aaa(*) = indgen(64*64);

g_as = yoga_acquisim(g_wfs,0);

num_ssp = (*(y_wfs(0)._isvalid))(*)(cum)(2:) * (*(y_wfs(0)._isvalid))(*);

//window, 0;
//acquisim_fillbcube, g_as, aaa;
//sensors_fillbinimage, g_wfs, 0;
//pli, sensors_getimg(g_wfs, 0);

//window, 1;
//acquisim_fillbcube, g_as, dist(64);
//sensors_fillbinimage, g_wfs, 0;
//pli, sensors_getimg(g_wfs, 0);

//window, 2;
//acquisim_fillbcube_2D, g_as, aaa, num_ssp;
//sensors_fillbinimage, g_wfs, 0;
//pli, sensors_getimg(g_wfs, 0);

//window, 3;
//acquisim_fillbcube_2D, g_as, dist(64), num_ssp;
//sensors_fillbinimage, g_wfs, 0;
//pli, sensors_getimg(g_wfs, 0);

acquisim_fillbcube, g_as, aaa;
tic; for(i=0; i<1000; i++) acquisim_fillbcube, g_as, aaa;
write, format="acquisim_fillbcube %f\n", tac();

acquisim_fillbcube_2D, g_as, aaa;
tic; for(i=0; i<1000; i++) acquisim_fillbcube_2D, g_as, aaa, num_ssp;
write, format="acquisim_fillbcube_2D %f\n", tac();

//b = array(int, 8,8)
//b(*) = num_ssp
//pli, b
//id = 1; 
//for(j=0; j<8; j++) { 
//  for(i=0; i<8; i++) { 
//    if(b(i+1,j+1)==id) {  
//      plt, swrite(format="%d",id), i+0.5,j+.5, tosys=1, color="red"; 
//      id++; 
//    } 
//  } 
//}
