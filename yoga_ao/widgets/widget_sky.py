#!/usr/bin/env python

import gobject
import gtk
import gtk.glade
import os
import sys
import time
import os, fcntl, errno

yoga_ao_top = os.environ['YOGA_AO_TOP']

class sky:

   def get_text(self,wdg):     
      self.pyk_resume(self.glade.get_widget(wdg).get_text())

   def get_value(self,wdg):      
      self.pyk_resume(str(self.glade.get_widget(wdg).get_value()))

   def destroy(self, wdg, data=None):
      self.py2yo('quit')
      raise SystemExit

   def __init__(self,path=os.path.join(yoga_ao_top, 'glade'), parent=None,py2yo=None):
      
      self.path=path
      
      self.gladefile = 'widget_sky.glade'
      
      self.glade = gtk.glade.XML(os.path.join(path,self.gladefile), root='top')
      self.top = self.glade.get_widget('top')

      # handle destroy event
      if (self.top):
         self.top.connect('destroy', self.destroy)

      self.glade.signal_autoconnect(self)

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

      #self.glade.get_widget('expander1').set_expanded(0)
      #self.glade.get_widget('expander2').set_expanded(0)
      self.glade.get_widget('layer_select').set_active(0)
      self.glade.get_widget('target_select').set_active(0)
   ######################################################
   # METHODS ASSOCIATED TO GLADE
   ######################################################
   def on_drawingarea_sky2_map(self,wdg,*args):
      drawingarea = self.glade.get_widget('drawingarea_sky')
      mwid1 = drawingarea.window.xid;
      drawingarea = self.glade.get_widget('drawingarea_sky1')
      mwid2 = drawingarea.window.xid;
      drawingarea = self.glade.get_widget('drawingarea_sky2')
      mwid3 = drawingarea.window.xid;
      self.py2yo('sky_win_init %d %d %d' % (mwid1,mwid2,mwid3))

   ######################################################
   # Atmosphere pane
   ######################################################
   def on_layer_select_changed(self,wdg):
      nlayer = wdg.get_active()
      self.y_update_layer(nlayer)

   def y_update_layer(self,nlayer):
      self.py2yo('update_layer_prop %d' % (nlayer))

   def y_init_layer(self,nlayer):
      alt       = self.glade.get_widget('alt').get_value()
      r0frac    = self.glade.get_widget('r0frac').get_value()
      windspeed = self.glade.get_widget('windspeed').get_value()
      winddir   = self.glade.get_widget('winddir').get_value()
      pupdiam = self.glade.get_widget('pupdiam').get_value()
      self.py2yo('init_layer_prop %d %f %f %f %f %d' % (nlayer,alt,r0frac,windspeed,winddir,pupdiam))

   def y_update_layer_gui(self,alt,r0frac,windspeed,winddir,dim_screen):
      self.glade.get_widget('alt').set_value(alt)
      self.glade.get_widget('r0frac').set_value(r0frac)
      self.glade.get_widget('windspeed').set_value(windspeed)
      self.glade.get_widget('winddir').set_value(winddir)
      self.glade.get_widget('screen_size').set_text(dim_screen)

   def y_init_atmos(self,dummy):
      teldiam = self.glade.get_widget('teldiam').get_value()
      pupdiam = self.glade.get_widget('pupdiam').get_value()
      zenith = self.glade.get_widget('zenith').get_value()
      nlayers = self.glade.get_widget('nlayers').get_value()
      r0 = self.glade.get_widget('r0').get_value()
      self.py2yo('create_atmos %f %d %d %d %f' % (teldiam,pupdiam,zenith,nlayers,r0))

   def on_nlayers_value_changed(self,wdg):
      nlayers = wdg.get_value()
      pupdiam = self.glade.get_widget('pupdiam').get_value()
      if (nlayers > 1):
         self.glade.get_widget('rm_screen').set_sensitive(1)
      else:
         self.glade.get_widget('rm_screen').set_sensitive(0)
      self.py2yo('update_nlayers %d %d' % (pupdiam,nlayers))
   
   def on_set_screen_clicked(self,wdg):
      nlayer = self.glade.get_widget('layer_select').get_active()
      self.y_init_layer(nlayer+1)

   def on_rm_screen_clicked(self,wdg):
      pupdiam = self.glade.get_widget('pupdiam').get_value()
      nlayer = self.glade.get_widget('layer_select').get_active()
      self.py2yo('remove_layer %d %d' % (nlayer+1,pupdiam))

   def on_default_atm_changed(self,wdg):
      type_conf = wdg.get_active()
      pupdiam = self.glade.get_widget('pupdiam').get_value()
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
      pupdiam = self.glade.get_widget('pupdiam').get_value()
      zenith = self.glade.get_widget('zenith').get_value()
      cobs = self.glade.get_widget('cobs').get_value()
      r0 = self.glade.get_widget('r0').get_value()
      freq = self.glade.get_widget('freq').get_value()
      self.py2yo('init_gatmos %d %f %f %f %f %f %d' % (pupdiam,zenith,teldiam,cobs,r0,freq,1))

   ######################################################
   # Target pane
   ######################################################
   def on_target_select_changed(self,wdg):
      ntarget = wdg.get_active()
      self.y_update_target(ntarget+1)

   def y_update_target(self,ntarget):
      self.py2yo('update_target_prop %d' % (ntarget))

   def y_init_target_prop(self,ntarget):
      mag     = self.glade.get_widget('mag').get_value()
      xpos    = self.glade.get_widget('xpos').get_value()
      ypos    = self.glade.get_widget('ypos').get_value()
      lambdai = self.glade.get_widget('lambda').get_value()
      self.py2yo('init_target_prop %d %f %f %f %f' % (ntarget,mag,xpos,ypos,lambdai))

   def y_init_target(self,dummy):
      self.py2yo('create_target %d' % (1))

   def y_update_target_gui(self,xpos,ypos,lambdai,mag):
      self.glade.get_widget('xpos').set_value(xpos)
      self.glade.get_widget('ypos').set_value(ypos)
      self.glade.get_widget('lambda').set_value(lambdai)
      self.glade.get_widget('mag').set_value(mag)

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
   # Main pane
   ######################################################
   def on_winselect_type_changed(self,wdg):
      mtype = wdg.get_active_text();
      self.py2yo('pyk_set skydisp_type "%s"' % mtype)
      self.py2yo('update_main_display2 "%s"' % (mtype))

   def on_winselect_number_changed(self,wdg):
      mtype = self.glade.get_widget('winselect_type').get_active_text();
      nlayer = wdg.get_active();
      self.py2yo('pyk_set skydisp_num %d' % nlayer)
      self.py2yo('update_main "%s" %d' % (mtype,nlayer))

   def on_enable_display_toggled(self,wdg):
      val = wdg.get_active()
      if (val == 1):
         mtype = self.glade.get_widget('winselect_type').get_active_text();
      else:
         mtype = "";
      self.py2yo('pyk_set skydisp_type "%s"' % mtype)

   def on_start_sky_clicked(self,wdg):
      self.py2yo('start_sky_loop')

   def on_sky_stop_clicked(self,wdg):
      self.py2yo('pyk_set skyloop %d' % 0)
      self.py2yo('animate %d' % 0)

   def on_initall_clicked(self,wdg):
      teldiam = self.glade.get_widget('teldiam').get_value()
      pupdiam = self.glade.get_widget('pupdiam').get_value()
      zenith = self.glade.get_widget('zenith').get_value()
      cobs = self.glade.get_widget('cobs').get_value()
      r0 = self.glade.get_widget('r0').get_value()
      freq = self.glade.get_widget('freq').get_value()
      self.py2yo('init_all %d %f %f %f %f %f' % (pupdiam,zenith,teldiam,cobs,r0,freq))

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
               sys.stderr.write('yo2py eval: '+str(e)+'\n')
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
   w = sky(parent=demo)
   gtk.main()
