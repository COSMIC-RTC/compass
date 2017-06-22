"""
Created on Tue Feb  2 09:39:35 2016

@author: fferreira

To launch it :

    - locally : 
        bokeh serve --show bokeh_display.py
    - as a server : 
        bokeh serve --port 8081 --host hippo6.obspm.fr:8081 bokeh_display.py
        then, open a web browser and connect to http://hippo6.obspm.fr:8081/bokeh_display.py
"""

import numpy as np
import glob
import os

import h5py
import datetime

from bokeh.plotting import Figure, figure
from bokeh.models import Range1d,ColumnDataSource
from bokeh.models.widgets import Select, Slider, CheckboxButtonGroup, Panel, Tabs, Button, Dialog, Paragraph
from bokeh.io import curdoc
from bokeh.models.layouts import HBox, VBox

import matplotlib.pyplot as plt
import matplotlib as mpl
######################################################################################
#  _      _ _      
# (_)_ _ (_) |_ ___
# | | ' \| |  _(_-<
# |_|_||_|_|\__/__/
######################################################################################     
class html_display:  
    def __del__(self):
        files = glob.glob("/home/fferreira/public_html/bokeh_display*")
        for f in files:
            os.remove(f)
            
    def __init__(self):


        self.datapath = "/home/fferreira/Data/"
        self.covmat = None
        self.f = h5py.File(self.datapath+"breakdown_classical.h5")
        self.f0 = h5py.File(self.datapath+"breakdown_imat_geom_meth0.h5")
        self.f1 = h5py.File(self.datapath+"breakdown_imat_geom_meth1.h5")
        self.fh = h5py.File(self.datapath+"breakdown_imat_hack.h5")
        
        self.files = {"Classical":self.f,"Imat geom 0":self.f0,"Imat geom 1":self.f1,"Imat hacked":self.fh}
        self.niter = self.f["com"][:].shape[1]
        self.nactus = self.f["com"][:].shape[0]
        self.meths = ["Classical","Imat geom 0","Imat geom 1","Imat hacked"]
        self.plot_type = ["Commands","Variance"]
        self.coms_list = list(self.f.keys())
        self.url = "http://hippo6.obspm.fr/~fferreira/bokeh_display"
        self.old = None
        
        ######################################################################################
        #         _    _          _      
        # __ __ _(_)__| |__ _ ___| |_ ___
        # \ V  V / / _` / _` / -_)  _(_-<
        #  \_/\_/|_\__,_\__, \___|\__/__/
        #               |___/            
        ######################################################################################
        self.dialog = Dialog(closable=False, visible= False, title="Dialog Box", content= "")        
        # Tab 1
        self.comsTags = Paragraph(text="Commands type", height=25)
        self.coms = CheckboxButtonGroup(labels=self.coms_list,active=[0])
        self.methTags = Paragraph(text="Matrix used", height=25)
        self.meth = CheckboxButtonGroup(labels=self.meths,active=[0])
        self.plot_select = Select(title="Plot type",value=self.plot_type[1],options=self.plot_type)
        self.iter_select = Slider(title="Iteration number",start=1,end=self.niter,step=1)
        # Tab 2
        self.A = Select(title="Commands A",value=self.coms_list[0],options=self.coms_list)
        self.B = Select(title="Commands B",value=self.coms_list[0],options=self.coms_list)
        self.imat = Select(title="Imat used",value=self.meths[0],options=self.meths)
        self.power = Slider(title="Abs(covmat)**X",start=1,end=10,step=1)
        self.cmin = Slider(title="vmin",start=1,end=10,step=1)
        self.cmax = Slider(title="vmax",start=1,end=10,step=1)
        self.cut = Button(label="Cut !",type="success")
        self.draw = Button(label="Draw !",type="success")
        
        self.colors = {"H_com":"green","alias_com":"blue","bp_com":"orange","com":"black","noise_com":"red","tomo_com":"purple","trunc_com":"cyan","wf_com":"magenta"}
        
        self.source1 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.source2 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.source3 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.source4 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.sourceLegend = ColumnDataSource(data=dict(legend=[]))
        
        #self.sourceC = ColumnDataSource(data=dict(C=[],x=[],y=[],dw=[],dh=[]))

        
        self.sources = {"Classical":self.source1,"Imat geom 0":self.source2,"Imat geom 1":self.source3,"Imat hacked":self.source4}

        self.p=Figure(plot_height=600, plot_width=800)#,tools=[hover])
        for c in self.colors:            
            self.p.line(legend=c,line_color=self.colors[c])
              
        self.p.multi_line("x","y",color="color",source=self.source1)
        self.p.multi_line("x","y",color="color",line_dash=[4,4],source=self.source2)
        self.p.multi_line("x","y",color="color",line_dash=[2,2],source=self.source3)
        self.p.multi_line("x","y",color="color",line_dash=[6,6],source=self.source4)
        self.p.xaxis.axis_label = "Actuators"
        
        self.xdr = Range1d(start=0,end=self.nactus)
        self.ydr = Range1d(start=0,end=self.nactus)
        self.p2=figure(x_range=self.xdr, y_range=self.ydr)
        self.p2.image_url(url=[], x=0, y=self.nactus,w=self.nactus,h=self.nactus)
        
        self.control_plot = [self.plot_select,self.iter_select]
        self.control_matrix = [self.A,self.B,self.imat,self.power]
        self.buttons = [self.coms,self.meth]
        for control in self.control_plot:
            control.on_change('value', self.update)
        for button in self.buttons:
            button.on_change('active',self.update)

        self.draw.on_click(self.update_matrix2)
        self.cut.on_click(self.cut_matrix)

        self.controlsTab1 = [self.coms,self.meth,self.plot_select,self.iter_select]
        self.controlsTab2 = [self.imat,self.A,self.B,self.power]
        self.inputs = HBox(VBox(self.comsTags,self.coms,self.methTags,self.meth,self.plot_select,self.iter_select), width=650)
        self.inputs2 = HBox(VBox(self.imat,self.A,self.B,self.power, self.draw, self.cmax,self.cmin,self.cut))#, width=350)
        self.tab1 = Panel(child=HBox(self.inputs,self.p), title="Commands")
        self.tab2 = Panel(child=HBox(self.inputs2,self.p2), title="Covmats")
        self.tabs = Tabs(tabs=[self.tab1,self.tab2])
        curdoc().clear()
        self.update(None,None,None)
        #self.update_matrix(None,None,None)

        curdoc().add_root(self.tabs)#hplot(inputs,p))#, p, p2)
        curdoc().add_root(self.dialog)
        
    
    ######################################################################################
    #   ___      _ _ _             _       
    #  / __|__ _| | | |__  __ _ __| |__ ___
    # | (__/ _` | | | '_ \/ _` / _| / /(_-<
    #  \___\__,_|_|_|_.__/\__,_\__|_\_\/__/
    #                                      
    ######################################################################################
    def update(self,attrname,old,new):
       # plot_val = plot_type.value
        self.source1.data = dict(x=[],y=[],color=[])
        self.source2.data = dict(x=[],y=[],color=[])
        self.source3.data = dict(x=[],y=[],color=[])
        self.source4.data = dict(x=[],y=[],color=[])
    
        coms_active = self.coms.active
        meth_active = self.meth.active
        plot_val = self.plot_select.value
        iteration = int(self.iter_select.value)
        
        for ii in meth_active:
            i = self.meths[ii]
            xi = []
            yi = []
            coloris = []
            for jj in coms_active:
                j = self.coms_list[jj]
                data=self.files[i][j][:]
                if(plot_val == "Commands"):
                    yi.append(data[:,iteration].tolist())
                    xi.append(list(range(len(data[:,iteration]))))
                    coloris.append(self.colors[j])
                    self.p.yaxis.axis_label = "Volts"
                    
                elif(plot_val == "Variance"):
                    yi.append(np.var(data,axis=1).tolist())
                    xi.append(list(range(len(np.var(data,axis=1)))))
                    coloris.append(self.colors[j])
                    self.p.yaxis.axis_label = "Variance"

            source = self.sources[i]
            source.data = dict(x=xi,y=yi,color=coloris)
            
        print("Plots updated")
    
        
    def update_matrix(self,attrname,old,new):
        A_val = self.A.value
        B_val = self.B.value
        ff = self.imat.value
        A_cov = self.files[ff][A_val][:]
        B_cov = self.files[ff][B_val][:]
        print("Values ok")
        self.covmat = (np.dot(A_cov,B_cov.T)/B_cov.shape[1])
        mpl.image.imsave("/home/fferreira/public_html/tmp.png",self.covmat)
        print("dot product ok")
        #self.sourceC.data = dict(url=[self.url],x=0,y=covmat.shape[0],dw=covmat.shape[0],dh=covmat.shape[0])
        #self.sourceC.data = dict(C=[covmat],x=[0],y=[0],dw=[covmat.shape[0]],dh=[covmat.shape[0]])
        
        print("Matrix updated")  
        
    def cut_matrix(self):
        self.dialog.visible = False
        vmin = self.cmin.value
        vmax = self.cmax.value
        self.dialog.content="Updating matrix..."
        self.dialog.visible = True
        if(self.old):
            os.remove(self.old)
        time = str(datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%f'))
        self.old = "/home/fferreira/public_html/bokeh_display"+time+".png"
        mpl.image.imsave(self.old,self.covmat,vmin=vmin,vmax=vmax)
        self.p2.image_url(url=dict(value=self.url+time+".png"), x=0, y=self.covmat.shape[0],w=self.covmat.shape[0],h=self.covmat.shape[0],retry_attempts=2,retry_timeout=1000)
        self.dialog.visible = False
        
    def update_matrix2(self): 
        self.dialog.visible = False
        if(self.old):
            os.remove(self.old)
        #self.draw.disabled = True
        A_val = self.A.value
        B_val = self.B.value
        ff = self.imat.value
        powa = np.ceil(self.power.value)
        self.dialog.content="Computing and loading matrix..."
        self.dialog.visible = True
        A_cov = self.files[ff][A_val][:]
        B_cov = self.files[ff][B_val][:]
        print("Values ok")
        self.covmat = (np.dot(A_cov,B_cov.T)/B_cov.shape[1])                
        print("dot product ok")
        if(powa != 1 and powa != 0):
            self.covmat = np.abs(self.covmat)**(1./powa)
            print("scale adjusted")
        self.cmin.start = self.covmat.min()
        self.cmin.end = self.covmat.max()
        self.cmin.value = self.cmin.start
        self.cmin.step = (self.cmin.end - self.cmin.start)/100.
        self.cmax.start = self.covmat.min()
        self.cmax.end = self.covmat.max()
        self.cmax.value = self.cmax.end     
        self.cmax.step = self.cmin.step
        time = str(datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%f'))
        self.old = "/home/fferreira/public_html/bokeh_display"+time+".png"
        mpl.image.imsave(self.old,self.covmat)
        self.p2.image_url(url=dict(value=self.url+time+".png"), x=0, y=self.covmat.shape[0],w=self.covmat.shape[0],h=self.covmat.shape[0],retry_attempts=2,retry_timeout=1000)

        #self.sourceC.data = dict(url=[self.url],x=0,y=covmat.shape[0],dw=covmat.shape[0],dh=covmat.shape[0])
        #self.draw.disabled = False
        print("Matrix updated2")
        self.dialog.visible = False
        

disp = html_display()
        

    
     # initial load of the data
    
    
