func update_dm_prop(numdm)
{
  extern dm_disp,y_dm;
  
  if (y_dm == []) {
    pyk,swrite(format=dm_disp._cmd+"y_init_dm(%d)",1);
  } else {
    if (numdm > numberof(y_dm)) error,"not a valid dm";
    else {
      if (y_dm(numdm).type == "pzt") typedm = 0;
      else if (y_dm(numdm).type == "kl") typedm = 2;
      else typedm = 1;
      
      if (y_dm(numdm).type == "kl")
        pyk,swrite(format=dm_disp._cmd+"y_update_dm_gui(%d, %f, %d , %f, %f, %f, %f)",
                   long(typedm),(y_dm.alt)(numdm),(y_dm.nkl)(numdm),
                   (y_dm.coupling)(numdm),(y_dm.hyst)(numdm),(y_dm.thresh)(numdm),
                   (y_dm.unitpervolt)(numdm));
      else 
        pyk,swrite(format=dm_disp._cmd+"y_update_dm_gui(%d, %f, %d , %f, %f, %f, %f)",
                   long(typedm),(y_dm.alt)(numdm),(y_dm.nact)(numdm),
                   (y_dm.coupling)(numdm),(y_dm.hyst)(numdm),(y_dm.thresh)(numdm),
                   (y_dm.unitpervolt)(numdm));
    }
  }
}

func init_dm_prop(numdm,typedm,alt,nactu,coupling,hyst,thresh,unitv)
{
  extern y_dm;

  if (y_dm == []) return;
  else {
    if (typedm == 0) y_dm(numdm).type = "pzt";
    if (typedm == 1) y_dm(numdm).type = "tt";
    if (typedm == 2) y_dm(numdm).type = "kl";
    
    y_dm(numdm).alt         = alt;
    if (typedm == 2) y_dm(numdm).nkl = nactu;
    else y_dm(numdm).nact   = nactu;
    y_dm(numdm).coupling    = coupling;
    y_dm(numdm).hyst        = hyst;
    y_dm(numdm).thresh      = thresh;
    y_dm(numdm).unitpervolt = unitv;
    y_dm(numdm).push4imat   = 1.0f / unitv;
  }
}

func create_dm(numdm,nsub)
{
  extern y_dm;
  
  y_dm  = array(dm_struct(),numdm);
  
  if (numdm > 1) {
    typedm = array(string,numdm);
    typedm(*) = "pzt";
    typedm(0) = "tt";
  } else {
    typedm = "pzt";
  }
  for (i=1;i<=numdm;i++) {
    y_dm(i).type          = typedm(i);
    y_dm(i).alt           = 0.;
    y_dm(i).nact          = nsub+1;
    y_dm(i).coupling      = 0.2;
    y_dm(i).hyst          = 0.1;
    y_dm(i).thresh        = 0.3;
    y_dm(i).unitpervolt   = 0.01;
    y_dm(i).push4imat     = 1.0f / y_dm(i).unitpervolt;
  }
  
  pyk,swrite(format=dm_disp._cmd+"y_dm_clear(%d)",1);    
  //pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').clear()");
  for (cc=1;cc<=numdm;cc++)
    pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').insert_text(%d,'%s')",
               0,swrite(format="DM # %d",numdm-cc+1));
  
  update_dm_prop,0;
}

func update_ndm(numdm,nsub)
{
  extern y_dm;
  
  if (y_dm != []) {
    if (numdm > numberof(y_dm)) {
      tmp = y_dm;
      tmp2 = array(dm_struct,numberof(y_dm)+1);
      tmp2(1:numberof(y_dm)) = y_dm;
      y_dm = tmp2;
      
      y_dm(0).type         = "pzt";
      y_dm(0).alt          = 0.;
      y_dm(0).nact         = nsub+1;
      y_dm(0).coupling     = 0.2;
      y_dm(0).hyst         = 0.1;
      y_dm(0).thresh       = 0.3;
      y_dm(0).unitpervolt  = 0.01;
      y_dm(0).push4imat    = 1.0f / y_dm(0).unitpervolt;
      
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').insert_text(%d,'%s')",
                 numberof(y_dm)-1,swrite(format="DM # %d",numberof(y_dm)));
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').set_active(%d)",numberof(y_dm)-1);
    
    } else if (numdm <  numberof(y_dm)) {
      y_dm = y_dm(1:-1);
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').remove_text(%d)",numberof(y_dm));
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').set_active(%d)",numberof(y_dm)-1);
    } else {
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').set_active(%d)",0);
    }
    if (numdm != numberof(y_dm)) update_ndm,numdm,nsub;
  }
}

func remove_dm(numdm)
{
  extern y_dm;

  if (y_dm == []) return;
  else {
    if (numberof(y_dm)>1) {
      if (numdm != numberof(y_dm)) {
        if (numdm == 1) {
          tmp = y_dm(2:);
          y_dm = tmp;
        } else {
          tmp = y_dm(1:-1);
          tmp(1:numdm-1) = y_dm(1:numdm-1);
          tmp(numdm:-1) = y_dm(numdm+1:);
          y_dm = tmp;
        }
      }
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('ndm').set_value(%d)",numberof(y_dm)-1);
      pyk,swrite(format=dm_disp._cmd+"y_dm_clear(%d)",1);    
      //pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').clear()");
      for (cc=1;cc<=numdm;cc++)
        pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').insert_text(%d,'%s')",
                   0,swrite(format="DM # %d",numdm-cc+1));
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').set_active(%d)",numberof(y_dm)-1);
    }
  }
}


func load_default_dm(tconf,nsub)
{
  extern y_dm;
  
  if (tconf == 1) {
    ndm = 2;
    y_dm                  = array(dm_struct(),ndm); // clean start
    
    y_dm(1).type          = "pzt";
    y_dm(1).nact          = nsub+1;
    y_dm(1).alt           = 0.;
    y_dm(1).thresh        = 0.3;
    y_dm(1).hyst          = 0.1;
    y_dm(1).coupling      = 0.2;
    y_dm(1).unitpervolt   = 0.01;
    y_dm(1).push4imat     = 100.0f;
    
    y_dm(2).type          = "tt";
    y_dm(2).alt           = 0.;
    y_dm(2).unitpervolt   = 0.005;
    y_dm(2).push4imat     = 10.0f;
  }
  if (tconf == 2) { 
    ndm = 3;
    y_dm                  = array(dm_struct(),ndm); // clean start
    
    y_dm(1:2).type          = "pzt";
    y_dm(1).nact            = nsub+1;
    y_dm(2).nact            = long((nsub+1)*1.3);
    y_dm(1).alt             = 0.;
    y_dm(2).alt             = 10000.;
    y_dm(1:2).thresh        = 0.3;
    y_dm(1:2).hyst          = 0.1;
    y_dm(1:2).coupling      = 0.2;
    y_dm(1:2).unitpervolt   = 0.01;
    y_dm(1:2).push4imat     = 100.0f;
    
    y_dm(3).type          = "tt";
    y_dm(3).alt           = 0.;
    y_dm(3).unitpervolt   = 0.005;
    y_dm(3).push4imat     = 10.0f;
   }
  if (tconf == 4) { 
    ndm = 1;
    y_dm                  = array(dm_struct(),ndm); // clean start
    
    y_dm(1).type          = "kl";
    y_dm(1).nkl           = long(nsub*nsub/1.5);
    y_dm(1).alt           = 0.;
    y_dm(1).thresh        = 0.3;
    y_dm(1).unitpervolt   = 0.01;
    y_dm(1).push4imat     = 0.01f;
  }
  if (tconf < 8) {
    pyk,swrite(format=dm_disp._cmd+"y_dm_clear(%d)",1);    
    //pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').clear()");
    for (cc=1;cc<=numberof(y_dm);cc++)
      pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').insert_text(%d,'%s')",0,
                 swrite(format="DM # %d",numberof(y_dm)-cc+1));
    pyk,swrite(format=dm_disp._cmd+"glade.get_widget('dm_select').set_active(%d)",0);
    pyk,swrite(format=dm_disp._cmd+"glade.get_widget('ndm').set_value(%d)",numberof(y_dm));
    
  }
}
