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

import h5py

from bokeh.plotting import Figure
from bokeh.models import Range1d,ColumnDataSource, HoverTool, CustomJS
from bokeh.models.widgets import Select, Slider, CheckboxButtonGroup, Panel, Tabs, Button
from bokeh.io import curdoc
from bokeh.mpl import to_bokeh
from bokeh.models.layouts import HBox, VBox, VBoxForm

import matplotlib.pyplot as pl
import matplotlib as mpl
######################################################################################
#  _      _ _      
# (_)_ _ (_) |_ ___
# | | ' \| |  _(_-<
# |_|_||_|_|\__/__/
######################################################################################     
class html_display:  
    def __init__(self):
        self.datapath = "/home/fferreira/Data/"
        
        self.f = h5py.File(self.datapath+"breakdown_classical.h5")
        self.f0 = h5py.File(self.datapath+"breakdown_imat_geom_meth0.h5")
        self.f1 = h5py.File(self.datapath+"breakdown_imat_geom_meth1.h5")
        self.fh = h5py.File(self.datapath+"breakdown_imat_hack.h5")
        
        self.files = {"Classical":self.f,"Imat geom 0":self.f0,"Imat geom 1":self.f1,"Imat hacked":self.fh}
        self.niter = self.f["com"][:].shape[1]
        self.nactus = self.f["com"][:].shape[0]
        self.meths = ["Classical","Imat geom 0","Imat geom 1","Imat hacked"]
        self.plot_type = ["Commands","Variance"]
        self.coms_list = self.f.keys()
        
        ######################################################################################
        #         _    _          _      
        # __ __ _(_)__| |__ _ ___| |_ ___
        # \ V  V / / _` / _` / -_)  _(_-<
        #  \_/\_/|_\__,_\__, \___|\__/__/
        #               |___/            
        ######################################################################################
        # Tab 1
        self.coms = CheckboxButtonGroup(labels=self.coms_list,active=[0])
        self.meth = CheckboxButtonGroup(labels=self.meths,active=[0])
        self.plot_select = Select(title="Plot type",value=self.plot_type[1],options=self.plot_type)
        self.iter_select = Slider(title="Iteration number",start=1,end=self.niter,step=1)
        # Tab 2
        self.A = Select(title="Commands A",value=self.coms_list[0],options=self.coms_list)
        self.B = Select(title="Commands B",value=self.coms_list[0],options=self.coms_list)
        self.imat = Select(title="Imat used",value=self.meths[0],options=self.meths)
        self.power = Slider(title="Abs(covmat)**X",start=1,end=10,step=1)
        self.draw = Button(label="Draw !",type="success")
        
        self.colors = {"H_com":"green","alias_com":"blue","bp_com":"orange","com":"black","noise_com":"red","tomo_com":"purple","trunc_com":"cyan","wf_com":"magenta"}
        
        self.source1 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.source2 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.source3 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.source4 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        
        self.sourceC = ColumnDataSource(data=dict(C=[],x=[],y=[],dw=[],dh=[]))

        
        self.sources = {"Classical":self.source1,"Imat geom 0":self.source2,"Imat geom 1":self.source3,"Imat hacked":self.source4}

        #
        self.p=Figure(plot_height=600, plot_width=800)#,tools=[hover])
        self.p.multi_line("x","y",color="color",source=self.source1)
        self.p.multi_line("x","y",color="color",line_dash=[4,4],source=self.source2)
        self.p.multi_line("x","y",color="color",line_dash=[2,2],source=self.source3)
        self.p.multi_line("x","y",color="color",line_dash=[6,6],source=self.source4)
        self.p.xaxis.axis_label = "Actuators"
        
        self.xdr = Range1d(start=0,end=self.nactus)
        self.ydr = Range1d(start=0,end=self.nactus)
        self.p2=Figure(x_range=self.xdr, y_range=self.ydr)
        self.p2.image(image="C",x="x",y="y",dw="dw",dh="dh",source=self.sourceC,palette="RdYlBu11")#"C",source=sourceC)#,x=0,y=0,dw=nactus,dh=nactus)#,source=sourceC)
        
        self.control_plot = [self.plot_select,self.iter_select]
        self.control_matrix = [self.A,self.B,self.imat,self.power]
        self.buttons = [self.coms,self.meth]
        for control in self.control_plot:
            control.on_change('value', self.update)
        for button in self.buttons:
            button.on_change('active',self.update)

        self.draw.on_click(self.update_matrix2)

        self.controlsTab1 = [self.coms,self.meth,self.plot_select,self.iter_select]
        self.controlsTab2 = [self.imat,self.A,self.B,self.power]
        self.inputs = HBox(VBox(self.coms,self.meth,self.plot_select,self.iter_select), width=650)
        self.inputs2 = HBox(VBox(self.imat,self.A,self.B,self.power, self.draw))#, width=350)
        self.tab1 = Panel(child=HBox(self.inputs,self.p), title="Commands")
        self.tab2 = Panel(child=HBox(self.inputs2,self.p2), title="Covmats")
        self.tabs = Tabs(tabs=[self.tab1,self.tab2])
        curdoc().clear()
        self.update(None,None,None)
        self.update_matrix(None,None,None)

        curdoc().add_root(self.tabs)#hplot(inputs,p))#, p, p2)
        
    
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
                    xi.append(range(len(data[:,iteration])))
                    coloris.append(self.colors[j])
                    self.p.yaxis.axis_label = "Volts"
                elif(plot_val == "Variance"):
                    yi.append(np.var(data,axis=1).tolist())
                    xi.append(range(len(np.var(data,axis=1))))
                    coloris.append(self.colors[j])
                    self.p.yaxis.axis_label = "Variance"
            source = self.sources[i]
            source.data = dict(x=xi,y=yi,color=coloris)
            
        print "Plots updated"
    
        
    def update_matrix(self,attrname,old,new): 
        A_val = self.A.value
        B_val = self.B.value
        ff = self.imat.value
        A_cov = self.files[ff][A_val][:]
        B_cov = self.files[ff][B_val][:]
        print "Values ok"
        covmat = (np.dot(A_cov,B_cov.T)/B_cov.shape[1])
        
        print "dot product ok"

        self.sourceC.data = dict(C=[covmat],x=[0],y=[0],dw=[covmat.shape[0]],dh=[covmat.shape[0]])
        
        print "Matrix updated"  
        
        
    def update_matrix2(self): 
        #self.draw.disabled = True
        A_val = self.A.value
        B_val = self.B.value
        ff = self.imat.value
        powa = np.ceil(self.power.value)
        A_cov = self.files[ff][A_val][:]
        B_cov = self.files[ff][B_val][:]
        print "Values ok"
        covmat = (np.dot(A_cov,B_cov.T)/B_cov.shape[1])
        mpl.image.imsave("tmp.png",covmat)
        
        print "dot product ok"
        if(powa != 1 and powa != 0):
            covmat = np.abs(covmat)**(1./powa)
            print "scale adjusted"
        self.sourceC.data = dict(C=[covmat],x=[0],y=[0],dw=[covmat.shape[0]],dh=[covmat.shape[0]])
        #self.draw.disabled = False
        print "Matrix updated"
        

disp = html_display()
        

    
     # initial load of the data
    
    
