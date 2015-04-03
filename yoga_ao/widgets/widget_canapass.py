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
      
      self.gladefile = 'widget_canapass.glade'
      
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

   ######################################################
   # Main pane
   ######################################################
   def on_toggle_atmos_toggled(self,wdg):
      self.py2yo('pyk_set y_see_atmos  %d' % wdg.get_active())

   def on_reset_strehl_clicked(self,wdg):
      target = self.glade.get_widget('reset_target').get_value_as_int()
      self.py2yo('target_reset_strehlmeter g_target 0')

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
      self.py2yo('pyk_set enable_display %d' % wdg.get_active())

   def on_start_ao_clicked(self,wdg):
      wdg.set_sensitive(0)
      self.py2yo('start_ao_loop')

   def on_ao_stop_clicked(self,wdg):
      self.py2yo('pyk_set aoloop %d' % 0)
      self.glade.get_widget('start_ao').set_sensitive(1)

   def on_initall_clicked(self,wdg):
      self.py2yo('init_all')

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
