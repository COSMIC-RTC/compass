#!/usr/bin/env python

import gobject
import gtk
import gtk.glade
import os
import sys
import time
import os, fcntl, errno

yoga_ao_top = os.environ['YOGA_AO_TOP']

class wao:

   def get_text(self,wdg):     
      self.pyk_resume(self.glade.get_widget(wdg).get_text())

   def get_value(self,wdg):      
      self.pyk_resume(str(self.glade.get_widget(wdg).get_value()))

   def destroy(self, wdg, data=None):
      self.py2yo('quit')
      raise SystemExit

   def __init__(self,path=os.path.join(yoga_ao_top, 'glade'), parent=None,py2yo=None):
      
      self.path=path
      
      self.gladefile = 'widget_ao.glade'
      
      self.gladet = gtk.glade.XML(os.path.join(path,self.gladefile))
      self.window = self.gladet.get_widget('window1')
      self.glade = gtk.glade.XML(os.path.join(path,self.gladefile), root='top')
      self.top = self.glade.get_widget('top')
      self.par_popup = self.gladet.get_widget('window2')
      if (self.par_popup):
         self.par_popup.connect('delete_event', self.on_popup_quit_activate)

      # handle destroy event
      if (self.top):
         self.top.connect('destroy', self.destroy)

      self.glade.signal_autoconnect(self)
      self.gladet.signal_autoconnect(self)

      if parent:
         parent.foreach(parent.remove)
         parent.add(self.top)

      # Do not open stdin event if not in standalone mode
      if not py2yo:
         # set stdin non blocking, this will prevent readline to block
         fd = sys.stdin.fileno()
         flags = fcntl.fcntl(fd, fcntl.F_GETFL)
         fcntl.fcntl(fd, fcntl.F_SETFL, flags | os.O_NONBLOCK)
         
         # add stdin to the event loop (yorick input pipe by spawn)
         gobject.io_add_watch(sys.stdin,gobject.IO_IN|gobject.IO_HUP,self.yo2py,None)

      #self.glade.get_widget('lgs_expander').set_expanded(0)
      #self.glade.get_widget('expander2').set_expanded(0)
      #self.glade.get_widget('layer_select').set_active(0)
      #self.glade.get_widget('target_select').set_active(0)
      self.pardir = "/"
      self.parfile = ""
      self.buffer  = gtk.TextBuffer()

   ######################################################
   # METHODS ASSOCIATED TO GLADE
   ######################################################
   def on_drawingarea_ao_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_ao')
      mwid1 = drawingarea.window.xid;
      self.py2yo('ao_win_init %d' % (mwid1))

   def on_drawingarea_turbu_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_turbu')
      mwid2 = drawingarea.window.xid;
      self.py2yo('turbu_win_init %d' % (mwid2))

   def on_drawingarea_lgs_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_lgs')
      mwid3 = drawingarea.window.xid;
      self.py2yo('lgs_win_init %d' % (mwid3))

   def on_drawingarea_centro_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_centro')
      mwid4 = drawingarea.window.xid;
      self.py2yo('centro_win_init %d' % (mwid4))

   def on_drawingarea_control_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_control')
      mwid5 = drawingarea.window.xid;
      self.py2yo('control_win_init %d' % (mwid5))

   def on_drawingarea_sky_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_sky')
      mwid6 = drawingarea.window.xid;
      self.py2yo('target_win_init %d' % (mwid6))

   def on_type_system_changed(self,wdg):
      nconf = wdg.get_active()
      teldiam = self.glade.get_widget('teldiam').get_value()
      self.py2yo('load_default_system %d %f' % (nconf+1,teldiam))
      
   def on_drawingarea_SR_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_SR')
      mwid7 = drawingarea.window.xid;
      self.py2yo('SR_win_init %d' % (mwid7))

   ######################################################
   # Atmosphere pane
   ######################################################
   def on_layer_select_changed(self,wdg):
      nlayer = wdg.get_active()
      self.y_update_layer(nlayer)

   def y_update_layer(self,nlayer):
      self.py2yo('update_layer_prop %d' % (nlayer))

   def y_layer_clear(self,clear):
      cmodel = self.glade.get_widget('layer_select').get_model()
      cmodel.clear()

   def y_init_layer(self,nlayer):
      alt       = self.glade.get_widget('alt').get_value()
      r0frac    = self.glade.get_widget('r0frac').get_value()
      l0        = self.glade.get_widget('l0').get_value()
      windspeed = self.glade.get_widget('windspeed').get_value()
      winddir   = self.glade.get_widget('winddir').get_value()
      pupdiam   = 128
      self.py2yo('init_layer_prop %d %f %f %f %f %f %d' % (nlayer,alt,r0frac,l0,windspeed,winddir,pupdiam))

   def y_update_layer_gui(self,alt,r0frac,l0,windspeed,winddir,dim_screen):
      self.glade.get_widget('alt').set_value(alt)
      self.glade.get_widget('r0frac').set_value(r0frac)
      self.glade.get_widget('l0').set_value(l0)
      self.glade.get_widget('windspeed').set_value(windspeed)
      self.glade.get_widget('winddir').set_value(winddir)
      self.glade.get_widget('screen_size').set_text(dim_screen)

   def y_init_atmos(self,dummy):
      teldiam = self.glade.get_widget('teldiam').get_value()
      pupdiam = 128
      zenith = self.glade.get_widget('zenith').get_value()
      nlayers = self.glade.get_widget('nlayers').get_value()
      r0 = self.glade.get_widget('r0').get_value()
      self.py2yo('create_atmos %f %d %d %d %f' % (teldiam,pupdiam,zenith,nlayers,r0))

   def on_nlayers_value_changed(self,wdg):
      nlayers = wdg.get_value()
      pupdiam = 128
      if (nlayers > 1):
         self.glade.get_widget('rm_screen').set_sensitive(1)
      else:
         self.glade.get_widget('rm_screen').set_sensitive(0)
      self.py2yo('update_nlayers %d %d' % (pupdiam,nlayers))
   
   def on_set_screen_clicked(self,wdg):
      nlayer = self.glade.get_widget('layer_select').get_active()
      self.y_init_layer(nlayer+1)

   def on_rm_screen_clicked(self,wdg):
      nlayer = self.glade.get_widget('layer_select').get_active()
      self.py2yo('remove_layer %d' % (nlayer+1))

   def on_default_atm_changed(self,wdg):
      type_conf = wdg.get_active()
      pupdiam = 128
      teldiam = self.glade.get_widget('teldiam').get_value()
      self.py2yo('load_default_atmos %d %d %f' % (type_conf+1,pupdiam,teldiam))

   def on_layers_plot_changed(self,wdg):
      type_plot = wdg.get_active()
      nlayer = self.glade.get_widget('layer_select').get_active()
      self.py2yo('select_plot_layers %d %d' % (type_plot,nlayer+1))

   def y_layers_plot_update(self,dummy):
      type_plot = self.glade.get_widget('layers_plot').get_active()
      nlayer = self.glade.get_widget('layer_select').get_active()
      self.py2yo('select_plot_layers %d %d' % (type_plot,nlayer+1))

   def on_init_atmos_clicked(self,wdg):
      teldiam = self.glade.get_widget('teldiam').get_value()
      pupdiam = 128
      zenith = self.glade.get_widget('zenith').get_value()
      cobs = self.glade.get_widget('cobs').get_value()
      r0 = self.glade.get_widget('r0').get_value()
      freq = self.glade.get_widget('freq').get_value()
      self.py2yo('init_gatmos %d %f %f %f %f %f %d' % (pupdiam,zenith,teldiam,cobs,r0,freq,1))

   ######################################################
   # WFS pane
   ######################################################
   def on_wfs_type_changed(self,clear):
      wtype = self.glade.get_widget('wfs_type').get_active()
      if wtype == 0:
         self.glade.get_widget('npix').set_sensitive(1)
         self.glade.get_widget('pixsize').set_sensitive(1)
      if wtype == 1:
         self.glade.get_widget('npix').set_sensitive(0)
         self.glade.get_widget('pixsize').set_sensitive(0)


   def y_wfs_clear(self,clear):
      cmodel = self.glade.get_widget('wfs_select').get_model()
      cmodel.clear()

   def on_wfs_select_changed(self,wdg):
      nwfs = wdg.get_active()
      self.y_update_wfs(nwfs+1)

   def y_update_wfs(self,nwfs):
      self.py2yo('update_wfs_prop %d' % (nwfs))

   def y_update_lgs(self,true):
      if true == 1:
         self.glade.get_widget('checklgs').set_active(True)
      if true == 0:
         self.glade.get_widget('checklgs').set_active(False)

   def y_init_wfs_prop(self,nwfs):
      typewfs    = self.glade.get_widget('wfs_type').get_active()
      lgsf       = self.glade.get_widget('checklgs').get_active()
      nsub       = self.glade.get_widget('nsub').get_value()
      npix       = self.glade.get_widget('npix').get_value()
      pixsize    = self.glade.get_widget('pixsize').get_value()
      frac       = self.glade.get_widget('frac').get_value()
      mag        = self.glade.get_widget('mag').get_value()
      xpos       = self.glade.get_widget('xpos').get_value()
      ypos       = self.glade.get_widget('ypos').get_value()
      zp         = self.glade.get_widget('zp').get_value()
      throughput = self.glade.get_widget('throughput').get_value()
      lambdai    = self.glade.get_widget('lambda').get_value()
      noise      = self.glade.get_widget('noise').get_value()
      self.py2yo('init_wfs_prop %d %d %d %d %f %f %f %f %f %f %f %f %f' % (typewfs,nwfs,nsub,npix,pixsize,mag,xpos,ypos,lambdai,frac,zp,throughput,noise))
      if (lgsf):
         gsalt    = self.glade.get_widget('lgsalt').get_value()
         lltx     = self.glade.get_widget('lltx').get_value()
         llty     = self.glade.get_widget('llty').get_value()
         power    = self.glade.get_widget('power').get_value()
         wreturn  = self.glade.get_widget('lgsreturn').get_value()
         proftype = self.glade.get_widget('lgs_prof').get_active_text()
         beam     = self.glade.get_widget('beam').get_value()
         self.py2yo('init_wfs_prop_lgs %f %f %f %f %f "%s" %f' % (gsalt,lltx,llty,power,wreturn,proftype,beam))

   def y_init_wfs(self,dummy):
      teldiam = self.glade.get_widget('teldiam').get_value()
      typewfs    = self.glade.get_widget('wfs_type').get_active()
      self.py2yo('create_wfs %d %d %d' % (1,typewfs,teldiam))

   def y_update_wfs_gui(self,typewfs,nsub,npix,pixsize,frac,xpos,ypos,lambdai,mag,zp,thhroughput,noise):
      self.glade.get_widget('wfs_type').set_active(typewfs)
      self.glade.get_widget('nsub').set_value(nsub)
      self.glade.get_widget('npix').set_value(npix)
      self.glade.get_widget('pixsize').set_value(pixsize)
      self.glade.get_widget('frac').set_value(frac)
      self.glade.get_widget('xpos').set_value(xpos)
      self.glade.get_widget('ypos').set_value(ypos)
      self.glade.get_widget('lambda').set_value(lambdai)
      self.glade.get_widget('mag').set_value(mag)
      self.glade.get_widget('zp').set_value(zp)
      self.glade.get_widget('throughput').set_value(thhroughput)
      self.glade.get_widget('noise').set_value(noise)

   def on_nwfs_value_changed(self,wdg):
      nwfs = wdg.get_value()
      teldiam = self.glade.get_widget('teldiam').get_value()
      if (nwfs > 1):
         self.glade.get_widget('rm_wfs').set_sensitive(1)
      else:
         self.glade.get_widget('rm_wfs').set_sensitive(0)
      self.py2yo('update_nwfs %d %d' % (nwfs,teldiam))

   def on_default_wfs_changed(self,wdg):
      type_conf = wdg.get_active()
      teldiam = self.glade.get_widget('teldiam').get_value()
      self.py2yo('load_default_wfs %d %f' % (type_conf+1,teldiam))

   def on_set_wfs_clicked(self,wdg):
      nwfs = self.glade.get_widget('wfs_select').get_active()
      self.y_init_wfs_prop(nwfs+1)

   def on_rm_wfs_clicked(self,wdg):
      nwfs = self.glade.get_widget('wfs_select').get_active()
      self.py2yo('remove_wfs %d' % (nwfs+1))

   def on_lgs_plot_changed(self,wdg):
      nwfs = self.glade.get_widget('wfs_select').get_active()
      type_conf = wdg.get_active()
      self.py2yo('lgs_plot_update %d %d' % (nwfs+1,type_conf+1))

   ######################################################
   # RTC pane
   ######################################################
   def y_centro_clear(self,clear):
      cmodel = self.glade.get_widget('centro_select').get_model()
      cmodel.clear()

   def on_centro_select_changed(self,wdg):
      ncentro = wdg.get_active()
      self.y_update_centro(ncentro+1)

   def on_type_centro_changed(self,wdg):
      typec = wdg.get_active()
      self.glade.get_widget('choose_centro').set_sensitive(0)
      self.glade.get_widget('fct_width').set_sensitive(0)
      self.glade.get_widget('label_width').set_sensitive(0)
      self.glade.get_widget('type_func').set_sensitive(0)
      self.glade.get_widget('label_fct').set_sensitive(0)
      self.glade.get_widget('nmax').set_sensitive(0)
      self.glade.get_widget('label_nmax').set_sensitive(0)
      self.glade.get_widget('thresh').set_sensitive(0)
      self.glade.get_widget('label_thresh').set_sensitive(0)
      if (typec == 1):
         self.glade.get_widget('thresh').set_sensitive(1)
         self.glade.get_widget('label_thresh').set_sensitive(1)
      if (typec == 2):
         self.glade.get_widget('nmax').set_sensitive(1)
         self.glade.get_widget('label_nmax').set_sensitive(1)
      if (typec == 3):
         self.glade.get_widget('type_func').set_sensitive(1)
         self.glade.get_widget('label_fct').set_sensitive(1)
      if (typec == 4):
         self.glade.get_widget('type_func').set_sensitive(1)
         self.glade.get_widget('label_fct').set_sensitive(1)

   def on_control_plot_changed(self,wdg):
      typec = wdg.get_active()
      self.py2yo('update_control_plot %d 0' % typec)
  
   def on_type_func_changed(self,wdg):   
      self.glade.get_widget('choose_centro').set_sensitive(0)
      self.glade.get_widget('fct_width').set_sensitive(0)
      self.glade.get_widget('label_width').set_sensitive(0)
      typef = wdg.get_active()
      if (typef == 0):
         self.glade.get_widget('fct_width').set_sensitive(1)
         self.glade.get_widget('label_width').set_sensitive(1)
      if (typef == 2):
         self.glade.get_widget('choose_centro').set_sensitive(1)

   def y_update_centro(self,ncentro):
      self.py2yo('update_centro_prop %d' % (ncentro))

   def y_init_centro(self,dummy):
      self.glade.get_widget('default_centro').set_active(0)

   def y_update_centro_gui(self,typec,nmax,thresh,type_fct,width):
      if (typec == "cog"):
         self.glade.get_widget('type_centro').set_active(0)
      if (typec == "tcog"):
         self.glade.get_widget('type_centro').set_active(1)
      if (typec == "bpcog"):
         self.glade.get_widget('type_centro').set_active(2)
      if (typec == "wcog"):
         self.glade.get_widget('type_centro').set_active(3)
      if (typec == "pyr"):
         self.glade.get_widget('type_centro').set_active(4)
      self.glade.get_widget('thresh').set_value(thresh)
      self.glade.get_widget('nmax').set_value(nmax)
      self.glade.get_widget('fct_width').set_value(width)
      if (type_fct == "gauss"):
         self.glade.get_widget('type_func').set_active(0)
      if (type_fct == "file"):
         self.glade.get_widget('type_func').set_active(1)
      if (type_fct == "model"):
         self.glade.get_widget('type_func').set_active(2)

   def on_default_centro_changed(self,wdg):
      type_centro = wdg.get_active_text()
      type_fct    = self.glade.get_widget('type_func').get_active_text()
      self.py2yo('load_default_centro "%s" "%s"' % (type_centro,type_fct))

   def on_set_centro_clicked(self,wdg):
      centro   = self.glade.get_widget('centro_select').get_active()
      typec    = self.glade.get_widget('type_centro').get_active_text()
      thresh   = self.glade.get_widget('thresh').get_value()
      nmax     = self.glade.get_widget('nmax').get_value()
      type_fct = self.glade.get_widget('type_func').get_active_text()
      width    = self.glade.get_widget('fct_width').get_value()
      self.py2yo('init_centro_prop %d "%s" %d %f "%s" %f' % (centro+1,typec,nmax,thresh,type_fct,width))

   def on_update_cmat_clicked(self,wdg):
      control = self.glade.get_widget('type_control').get_active()
      maxcond = self.glade.get_widget('maxcond').get_value()
      self.py2yo('cmat_update %d %f' % (control,maxcond))

   def on_spydrit_clicked(self,wdg):
      typec = self.glade.get_widget('control_plot').get_active()
      self.py2yo('spydr_control_plot %d 0' % typec)

   def y_update_control(self,dummy):
      self.py2yo('update_control_prop')

   def y_init_control(self,dummy):
      self.glade.get_widget('default_control').set_active(0)

   def y_update_control_gui(self,typec,maxcond,delay,gain):
      if (typec == "ls"):
         self.glade.get_widget('type_control').set_active(0)
      if (typec == "mv"):
         self.glade.get_widget('type_control').set_active(1)
      if (typec == "la"):
         self.glade.get_widget('type_control').set_active(2)
      if (typec == "lqg"):
         self.glade.get_widget('type_control').set_active(3)
      self.glade.get_widget('maxcond').set_value(maxcond)
      self.glade.get_widget('TTcond').set_value(maxcond)
      self.glade.get_widget('delay').set_value(delay)
      self.glade.get_widget('loop_gain').set_value(gain)

   def on_default_control_changed(self,wdg):
      type_control = wdg.get_active_text()
      
      if (type_control == "MV"):
         self.glade.get_widget('label68').set_visible(1)
         self.glade.get_widget('TTcond').set_visible(1)
         self.glade.get_widget('filterTT').set_visible(1)
      else:
         self.glade.get_widget('label68').set_visible(0)
         self.glade.get_widget('TTcond').set_visible(0)
         self.glade.get_widget('filterTT').set_visible(0)
         
      self.py2yo('load_default_control "%s"' % (type_control))

   def on_type_control_changed(self,wdg):
      type_control = wdg.get_active_text()
      if (type_control == "MV"):
         self.glade.get_widget('label68').set_visible(1)
         self.glade.get_widget('TTcond').set_visible(1)
         self.glade.get_widget('filterTT').set_visible(1)
      else:
         self.glade.get_widget('label68').set_visible(0)
         self.glade.get_widget('TTcond').set_visible(0)
         self.glade.get_widget('filterTT').set_visible(0)

   def on_set_control_clicked(self,wdg):
      control = self.glade.get_widget('type_control').get_active()
      maxcond = self.glade.get_widget('maxcond').get_value()
      delay   = self.glade.get_widget('delay').get_value()
      gain    = self.glade.get_widget('loop_gain').get_value()
      self.py2yo('init_control_prop %d %f %f %f' % (control+1,maxcond,delay,gain))

   def on_filterTT_clicked(self,wdg):
      maxcond = self.glade.get_widget('TTcond').get_value()
      self.py2yo('cmat_filterTT %f' % maxcond)

   def on_cmat_up_clicked(self,wdg):
      control = self.glade.get_widget('type_control').get_active()
      maxcond = self.glade.get_widget('maxcond').get_value()
      delay   = self.glade.get_widget('delay').get_value()
      gain    = self.glade.get_widget('loop_gain').get_value()
      self.py2yo('update_control %d %f %f %f' % (control+1,maxcond,delay,gain))


   ######################################################
   # Target pane
   ######################################################
   def y_target_clear(self,clear):
      cmodel = self.glade.get_widget('target_select').get_model()
      cmodel.clear()

   def on_target_select_changed(self,wdg):
      ntarget = wdg.get_active()
      self.y_update_target(ntarget+1)

   def y_update_target(self,ntarget):
      self.py2yo('update_target_prop %d' % (ntarget))

   def y_init_target_prop(self,ntarget):
      mag     = self.glade.get_widget('magt').get_value()
      xpos    = self.glade.get_widget('xpost').get_value()
      ypos    = self.glade.get_widget('ypost').get_value()
      lambdai = self.glade.get_widget('lambdat').get_value()
      self.py2yo('init_target_prop %d %f %f %f %f' % (ntarget,mag,xpos,ypos,lambdai))

   def y_init_target(self,dummy):
      self.py2yo('create_target %d' % (1))

   def y_update_target_gui(self,xpos,ypos,lambdai,mag):
      self.glade.get_widget('xpost').set_value(xpos)
      self.glade.get_widget('ypost').set_value(ypos)
      self.glade.get_widget('lambdat').set_value(lambdai)
      self.glade.get_widget('magt').set_value(mag)

   def y_target_plot_update(self,dummy):
      ntarget = self.glade.get_widget('target_select').get_active()
      self.py2yo('target_plot_update %d' % (ntarget+1))

   def on_ntargets_value_changed(self,wdg):
      ntargets = wdg.get_value()
      if (ntargets > 1):
         self.glade.get_widget('rm_source').set_sensitive(1)
      else:
         self.glade.get_widget('rm_source').set_sensitive(0)
      self.py2yo('update_ntargets %d' % (ntargets))

   def on_default_target_changed(self,wdg):
      type_conf = wdg.get_active()
      self.py2yo('load_default_target %d' % (type_conf+1))

   def on_set_source_clicked(self,wdg):
      ntarget = self.glade.get_widget('target_select').get_active()
      self.y_init_target_prop(ntarget+1)

   def on_rm_source_clicked(self,wdg):
      ntarget = self.glade.get_widget('target_select').get_active()
      self.py2yo('remove_target %d' % (ntarget+1))

   ######################################################
   # DM pane
   ######################################################
   def y_dm_clear(self,clear):
      cmodel = self.glade.get_widget('dm_select').get_model()
      cmodel.clear()

   def on_dm_select_changed(self,wdg):
      ndm = wdg.get_active()
      self.y_update_dm(ndm+1)

   def y_update_dm(self,ndm):
      self.py2yo('update_dm_prop %d' % (ndm))

   def y_init_dm_prop(self,ndm):
      typedm     = self.glade.get_widget('type_dm').get_active()
      alt        = self.glade.get_widget('dmalt').get_value()
      nactu      = self.glade.get_widget('nactu').get_value()
      coupling   = self.glade.get_widget('coupling').get_value()
      hyst       = self.glade.get_widget('hyst').get_value()
      thresh     = self.glade.get_widget('dm_thresh').get_value()
      unitv      = self.glade.get_widget('unitv').get_value()
      self.py2yo('init_dm_prop %d %d %f %d %f %f %f %f' % (ndm,typedm,alt,nactu,coupling,hyst,thresh,unitv))
 
   def y_init_dm(self,dummy):
      nsub = self.glade.get_widget('nsub').get_value()
      self.py2yo('create_dm %d %d' % (1,nsub))

   def y_update_dm_gui(self,type_dm,dmalt,nactu,coupling,hyst,dm_thresh,unitv):
      self.glade.get_widget('type_dm').set_active(type_dm)
      self.glade.get_widget('dmalt').set_value(dmalt)
      self.glade.get_widget('nactu').set_value(nactu)
      self.glade.get_widget('coupling').set_value(coupling)
      self.glade.get_widget('hyst').set_value(hyst)
      self.glade.get_widget('dm_thresh').set_value(dm_thresh)
      self.glade.get_widget('unitv').set_value(unitv)
 
   def on_ndm_value_changed(self,wdg):
      ndm = wdg.get_value()
      nsub = self.glade.get_widget('nsub').get_value()
      if (ndm > 1):
         self.glade.get_widget('rm_dm').set_sensitive(1)
      else:
         self.glade.get_widget('rm_dm').set_sensitive(0)
      self.py2yo('update_ndm %d %d' % (ndm,nsub))

   def on_default_dm_changed(self,wdg):
      type_conf = wdg.get_active()
      nsub = self.glade.get_widget('nsub').get_value()
      self.py2yo('load_default_dm %d %f' % (type_conf+1,nsub))

   def on_set_dm_clicked(self,wdg):
      ndm = self.glade.get_widget('dm_select').get_active()
      self.y_init_dm_prop(ndm+1)

   def on_rm_dm_clicked(self,wdg):
      ndm = self.glade.get_widget('dm_select').get_active()
      self.py2yo('remove_dm %d' % (ndm+1))

   def on_reset_dm_clicked(self,wdg):
      ndm = self.glade.get_widget('dm_select').get_active()
      self.py2yo('reset_dm %d' % (ndm+1))

   ######################################################
   # Main pane
   ######################################################
   def on_reset_strehl_clicked(self,wdg):
      target = self.glade.get_widget('reset_target').get_value_as_int()
      self.py2yo('target_reset_strehlmeter g_target %d'%target)

   def y_win_clear(self,clear):
      cmodel = self.glade.get_widget('winselect_type').get_model()
      cmodel.clear()

   def y_winnum_clear(self,clear):
      cmodel = self.glade.get_widget('winselect_number').get_model()
      cmodel.clear()

   def on_device_value_changed(self,wdg):
      ndev = wdg.get_value()
      self.py2yo('activeDevice %d' % ndev)

   def on_unzoom_button_clicked(self,wdg):
      self.py2yo('ao_unzoom')

   def on_winselect_type_changed(self,wdg):
      mtype = wdg.get_active_text()
      self.py2yo('pyk_set aodisp_type "%s"' % mtype)
      self.py2yo('update_main_display2 "%s"' % (mtype))

   def on_winselect_number_changed(self,wdg):
      mtype = self.glade.get_widget('winselect_type').get_active_text()
      nlayer = wdg.get_active()
      if (nlayer < 0):
         nlayer = 0
      self.py2yo('pyk_set aodisp_num %d' % nlayer)
      self.py2yo('update_main "%s" %d' % (mtype,nlayer))
      self.py2yo('pyk_set aodisp_num %d' % nlayer)

   def on_enable_display_toggled(self,wdg):
      val = wdg.get_active()
      if (val == 1):
         disp = 1;
      else:
         disp = 0;
      self.py2yo('pyk_set enable_display %d' % disp)

   def on_SRmonitor_toggled(self,wdg):
      val = wdg.get_active()
      if (val == 1):
         disp = 1;
      else:
         disp = 0;
      self.glade.get_widget('drawingarea_SR').set_visible(disp)
      self.py2yo('pyk_set SR_monitor %d' % disp)

   def on_start_ao_clicked(self,wdg):
      wdg.set_sensitive(0)
      self.py2yo('start_ao_loop')

   def on_ao_stop_clicked(self,wdg):
      self.py2yo('pyk_set aoloop %d' % 0)
      self.py2yo('animate %d' % 0)
      self.glade.get_widget('start_ao').set_sensitive(1)

   def on_initall_clicked(self,wdg):
      teldiam = self.glade.get_widget('teldiam').get_value()
      zenith = self.glade.get_widget('zenith').get_value()
      cobs = self.glade.get_widget('cobs').get_value()
      r0 = self.glade.get_widget('r0').get_value()
      freq = self.glade.get_widget('freq').get_value()
      self.py2yo('init_all %f %f %f %f %f' % (zenith,teldiam,cobs,r0,freq))

   def on_imagemenuitem1_activate(self,wdg):
      self.par_popup.show()

   def on_imagemenuitem2_activate(self,wdg):
      chooser = gtk.FileChooserDialog(title='parameter file selection',action=gtk.FILE_CHOOSER_ACTION_OPEN,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))
      filter = gtk.FileFilter()
      filter.add_pattern('*.par')
      filter.set_name('parfiles')
      chooser.add_filter(filter)
      chooser.set_current_folder(os.path.join(yoga_ao_top, 'data/par'))
      res = chooser.run()
      if res == gtk.RESPONSE_OK:
         filename  = chooser.get_filename()
         self.pardir = os.path.dirname(filename)
         self.parfile = os.path.basename(filename)
         self.window.set_title(self.parfile)
         self.py2yo('load_parfile "%s" "%s"' % (filename,self.parfile))
      chooser.destroy()
      #self.par_popup.show()

   def on_imagemenuitem10_activate(self,wdg):
      dialog = self.gladet.get_widget('aboutdialog1')
      dialog.run()
      dialog.hide()

   def on_next_clicked(self,wdg):
      self.py2yo('pyk_set aoloop %d' % 1)
      self.py2yo('ao_loop %d' % (1))
      self.py2yo('pyk_set aoloop %d' % 0)

   def on_imagemenuitem7_activate(self,wdg):
      self.par_popup.set_title(self.parfile)
      self.par_popup.show()
      f = open(self.pardir+'/'+self.parfile,'r')
      params = f.read()
      f.close()
      self.buffer.set_text(params)
      self.buffer.connect('changed',self.on_modified_editor)
      textarea = self.gladet.get_widget('textarea_edit')
      textarea.set_buffer(self.buffer)

   def on_modified_editor(self,wdg):
      self.par_popup.set_title(' * '+self.parfile)   

   def on_editor_save_activate(self,wdg):
      # get buffer content
      textarea = self.gladet.get_widget('textarea_edit')
      self.buffer = textarea.get_buffer()
      params=self.buffer.get_text(self.buffer.get_start_iter(),self.buffer.get_end_iter())
      # save file
      f = open(self.pardir+'/'+self.parfile,'w')
      self.pardir+'/'+self.parfile
      f.write(params)
      f.close()
      # refresh editor window title
      self.par_popup.set_title(self.parfile)
      self.py2yo('load_parfile "%s" "%s"' % (self.pardir+'/'+self.parfile,self.parfile))

   def on_popup_quit_activate(self,wdg,*args):
      self.par_popup.hide()
      return True

   def on_save_as_activate(self,wdg):
      chooser = gtk.FileChooserDialog(title='Save parfile as...',action=gtk.FILE_CHOOSER_ACTION_SAVE,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_SAVE,gtk.RESPONSE_OK))
      res = chooser.run()
      if res == gtk.RESPONSE_CANCEL:
         chooser.destroy()
         return
      name = chooser.get_filename()
      self.pardir = os.path.dirname(name)
      self.parfile = os.path.basename(name)
      chooser.destroy()
      # get buffer content
      textarea = self.gladet.get_widget('textarea_edit')
      self.buffer = textarea.get_buffer()
      params=self.buffer.get_text(self.buffer.get_start_iter(),self.buffer.get_end_iter())
      # save the content in tmp.par
      f = open(self.pardir+'/'+self.parfile,'w')
      f.write(params)
      f.close()
      # get new name
      self.py2yo('load_parfile "%s" "%s"' % (self.pardir+'/'+self.parfile,self.parfile))

   def on_item_quit_activate(self,wdg):
      self.py2yo('quit')

   ######################################################
   # THIS IS WHERE YOU PLACE METHODS ASSOCIATED TO GLADE
   ######################################################

   # Example of a button

   #def on_button_test_clicked(self,wdg):
   #   self.py2yo('hello_yorick')

   #def update_status_test(self,txt):
   #   self.glade.get_widget('statusbar_test').push(1,txt)

   #def update_txt_test(self,txt):
   #   self.glade.get_widget('entry_test').set_text(txt)


   ###################################################
   # END OF GLADE RELATED METHOD DEFINITION
   # minimal wrapper for yorick/python communication
   ###################################################
      
   def yo2py_flush(self):
      sys.stdin.flush()
   
   def py2yo(self,msg):
      # sends string command to yorick's eval
      sys.stdout.write(msg+'\n')
      sys.stdout.flush()
   
   def yo2py(self,cb_condition,*args):
      if cb_condition == gobject.IO_HUP:
         raise SystemExit, "lost pipe to yorick"
      # handles string command from yorick
      # note: individual message needs to end with /n for proper ungarbling
      while 1:
         try:
            msg = sys.stdin.readline()
            msg = "self."+msg
            #  self.py2yo('\"%s\"' % msg)
            try:
               exec(msg)
            except Exception, e:
               exc_type, exc_obj, exc_tb = sys.exc_info()
               sys.stderr.write('yo2py eval: '+str(e)+'\n')
               if exc_type == SyntaxError:
                  raise SystemExit, "lost pipe to yorick"
         except IOError, e:
            if e.errno == errno.EAGAIN:
               # the pipe's empty, good
               break
            # else bomb out
            raise SystemExit, "yo2py unexpected IOError:" + str(e)
         except Exception, ee:
            raise SystemExit, "yo2py unexpected Exception:" + str(ee)
      # carefull with the ident here
      return True

####################################
# CODE FOR STANDALONE MODE
####################################
if __name__ == "__main__":
   #print 'standalone demo'
   demo=gtk.Window(type=gtk.WINDOW_TOPLEVEL)
   demo.connect('destroy', gtk.main_quit)
   demo.show()
   w = wao(parent=demo)
   gtk.main()
