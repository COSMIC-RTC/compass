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
import pandas
import datetime

from bokeh.plotting import Figure, figure
from bokeh.models import Range1d,ColumnDataSource, HoverTool
from bokeh.models.widgets import Select, Slider, CheckboxButtonGroup, Panel, Tabs, Button, Dialog, Paragraph, RadioButtonGroup,TextInput
from bokeh.io import curdoc
from bokeh.models.layouts import HBox, VBox
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn

import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.sparse import csr_matrix
######################################################################################
#  _      _ _
# (_)_ _ (_) |_ ___
# | | ' \| |  _(_-<
# |_|_||_|_|\__/__/
######################################################################################
class html_display:
    def __del__(self):
        files = glob.glob("/home/fferreira/public_html/breakdown_display*")
        for f in files:
            os.remove(f)

    def __init__(self):

        self.datapath = "/home/fferreira/Data/"
        self.covmat = None
        self.files = glob.glob("/home/fferreira/Data/breakdown_*.h5")
        self.files.sort()
        self.f_list = []
        for f in self.files:
            self.f_list.append(f.split('/')[-1])

        self.f = h5py.File(self.files[0])
        if(self.f.attrs.keys().count("target.Lambda")):
            self.Lambda_tar = self.f.attrs["target.Lambda"][0]
        else:
            self.Lambda_tar = 1.65

        self.Btt = self.f["Btt"][:]
        if(self.f.keys().count("IF")): #Dense case
            self.IF = self.f["IF"][:]
        else: #Sparse case
            self.IF = csr_matrix((self.f["IF.data"][:],self.f["IF.indices"][:],self.f["IF.indptr"][:]))
            self.IF = self.IF.T
        self.P = self.f["P"][:]#/np.sqrt(self.IF.shape[0])
        self.modes = self.IF.dot(self.Btt)#np.dot(self.f["IF"][:],self.Btt)
        self.swap = np.arange(self.modes.shape[1])-2
        self.swap[0:2] = [self.modes.shape[1]-2,self.modes.shape[1]-1]
#        self.modes = self.modes[:,self.swap]

        self.indx_pup = self.f["indx_pup"][:]
        self.pup = np.zeros((self.f["dm_dim"].value,self.f["dm_dim"].value))

        self.niter = self.f["com"][:].shape[1]
        self.nactus = self.f["com"][:].shape[0]
        self.nmodes = self.P.shape[0]

        self.plot_type = ["Commands","Variance"]
        self.coms_list = self.f.keys()
        self.coms_list.remove("Btt")
        self.coms_list.remove("P")
        self.coms_list.remove("dm_dim")
        self.coms_list.remove("indx_pup")
        if(self.f.keys().count("IF")): #Dense case
            self.coms_list.remove("IF")
        else:
            self.coms_list.remove("IF.data")
            self.coms_list.remove("IF.indices")
            self.coms_list.remove("IF.indptr")
        if(self.f.keys().count("SR")):
            self.coms_list.remove("SR")
        if(self.f.keys().count("SR2")):
            self.coms_list.remove("SR2")
        if(self.f.keys().count("fit_error")):
            self.coms_list.remove("fit_error")
        if(self.f.keys().count("cov")):
            self.cov = self.f["cov"][:]
            self.cor = self.f["cor"][:]
            self.coms_list.remove("cov")
            self.coms_list.remove("cor")
        else:
            self.cov, self.cor = self.cov_cor()

        self.basis = ["Actuators","Btt"]
        self.url = "http://hippo6.obspm.fr/~fferreira/breakdown_display"
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
        self.DB_select = Select(title="Database",value=self.f_list[0],options=self.f_list)
        self.DB_button = Button(label="Load DB",type="success")
        self.plot_select = Select(title="Plot type",value=self.plot_type[1],options=self.plot_type)
        self.basis_select1 = Select(title="Basis",value=self.basis[0],options=self.basis)
        self.iter_select = Slider(title="Iteration number",start=1,end=self.niter,step=1)
        self.plusTag = Paragraph(text="Add :", height=25)
        self.plus_select = CheckboxButtonGroup(labels=self.coms_list+["fitting"],active=[0,1,2,4,5,6,8])
        self.moinsTag = Paragraph(text="Substract :", height=25)
        self.moins_select = CheckboxButtonGroup(labels=self.coms_list+["fitting"],active=[])
        self.diff_button = Button(label="Sum !",type="success")
        # Tab 2
        self.A = Select(title="Commands A",value=self.coms_list[0],options=self.coms_list)
        self.B = Select(title="Commands B",value=self.coms_list[0],options=self.coms_list)
        self.basis_select2 = Select(title="Basis",value=self.basis[0],options=self.basis)
        self.power = Slider(title="Abs(covmat)**X",start=0.1,end=1.,step=0.1,value=1.)
        self.cmin = Slider(title="vmin",start=1,end=10,step=1)
        self.cmax = Slider(title="vmax",start=1,end=10,step=1)
        self.rescale = Button(label="Rescale !",type="primary")
        self.draw = Button(label="Draw !",type="success")
        self.diag = Button(label="Plot diag !",type="primary")
        self.cut = Button(label="Cut !",type="primary")
        self.axiscut = Slider(title="X/Y cut",start=0,end=1,step=1)
        self.XY = RadioButtonGroup(labels=["X","Y"],active=0)
        self.DataTableItems = ["Type","Noise","Truncature","Aliasing","FilteredModes","Bandwidth","Tomography"]
        self.table_cov_source = ColumnDataSource(data=dict(Type=[],
                                                       Noise=[],
                                                        Truncature=[],
                                                        Aliasing=[],
                                                        FilteredModes=[],
                                                        Bandwidth=[],
                                                        Tomography=[]))
        self.table_cor_source = ColumnDataSource(data=dict(Type=[],
                                                      Noise=[],
                                                        Truncature=[],
                                                        Aliasing=[],
                                                        FilteredModes=[],
                                                        Bandwidth=[],
                                                        Tomography=[]))
        self.cov_table, self.cor_table = self.createDataTables()
        self.updateDataTables()
        # Tab 3
        self.basis_select3 = Select(title="Basis",value=self.basis[0],options=self.basis)
        self.modes_select = Slider(title="Mode #",value=0,start=0,end=self.modes.shape[1],step=1)
        #self.modes_select = TextInput(value="0:"+str(self.modes.shape[1]-1),title="Enter a mode to display")
        self.draw_mode = Button(label="Draw !",type="success")
        self.inc_mode = Button(label="+",type="primary")
        self.desinc_mode = Button(label="-",type="primary")

        self.colors = {"H_com":"green","bp_com":"orange",
                       "com":"black","noise_com":"red","tomo_com":"purple",
                       "trunc_com":"cyan","wf_com":"magenta",
                       "alias_wfs_com":"blue"}

        self.source1 = ColumnDataSource(data=dict(x=[], y=[],color=[],typec=[]))
        self.source2 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
        self.source3 = ColumnDataSource(data=dict(x=[], y=[],color=[]))

        self.hover = HoverTool(tooltips=[("x","@x"),("y","@y"),("type","@typec")])
        self.hoverlog = HoverTool(tooltips=[("x","@x"),("y","@y"),("type","@typec")])
        TOOLS = "resize,save,pan,box_zoom,tap, box_select, wheel_zoom, lasso_select,reset"

        self.p=Figure(plot_height=600, plot_width=800,tools=[TOOLS,self.hover])
        self.plog=Figure(plot_height=600, plot_width=800,y_range=[1e-6,10],y_axis_type="log",tools=[TOOLS,self.hoverlog])
        self.psum=Figure(plot_height=600, plot_width=800)
        for c in self.colors:
            self.p.line(legend=c,line_color=self.colors[c])
            self.plog.line(legend=c,line_color=self.colors[c])

        self.p.multi_line("x","y",color="color",source=self.source1)
        self.plog.multi_line("x","y",color="color",source=self.source1)
        self.psum.line(legend="Image SR", line_color="red")
        self.psum.line(legend="Phase SR ", line_color="purple")
        self.psum.line(legend="Var(X+Y)", line_color="blue")
        self.psum.line(legend="Var(X)+var(Y)", line_color="green")

        self.psum.multi_line("x","y",color="color",source=self.source3)
        self.psum.yaxis.axis_label = "Strehl Ratio"

        self.xdr = Range1d(start=0,end=self.nactus)
        self.ydr = Range1d(start=self.nactus,end=0)
        self.p2=figure(x_range=self.xdr, y_range=self.ydr, x_axis_location="above")
        self.p2.image_url(url=[], x=0, y=0,w=self.nactus,h=self.nactus)
        self.p3=Figure(plot_height=600, plot_width=800)
        self.p3.line(x="x",y="y",source=self.source2)

        self.xdr2 = Range1d(start=0,end=self.pup.shape[0])
        self.ydr2 = Range1d(start=self.pup.shape[1],end=0)
        self.pmodes=figure(x_range=self.xdr2, y_range=self.ydr2, x_axis_location="above")
        self.pmodes.image_url(url=[], x=0, y=0,w=self.pup.shape[0],h=self.pup.shape[1])

        self.control_plot = [self.plot_select,self.iter_select,self.basis_select1]

        self.buttons = [self.coms]
        for control in self.control_plot:
            control.on_change('value', self.update)
        for button in self.buttons:
            button.on_change('active',self.update)

        self.draw.on_click(self.update_matrix2)
        self.draw_mode.on_click(self.update_mode)
        self.rescale.on_click(self.rescale_matrix)
        self.diag.on_click(self.get_diag)
        self.cut.on_click(self.cut_matrix)
        self.inc_mode.on_click(self.mode_increment)
        self.desinc_mode.on_click(self.mode_desincrement)
        self.diff_button.on_click(self.plot_sum)
        self.DB_button.on_click(self.loadDB)


        self.inputs = HBox(VBox(self.DB_select,self.DB_button,self.comsTags,self.coms,self.plot_select,self.basis_select1,self.iter_select,self.plusTag,self.plus_select,self.moinsTag,self.moins_select,self.diff_button), width=650)
        self.inputs2 = HBox(VBox(self.DB_select,self.DB_button,self.basis_select2,self.A,self.B,self.power, self.draw,self.cmax,self.cmin,self.rescale,self.axiscut,self.XY,self.cut,self.diag))#, width=350)
        self.inputs3 = HBox(VBox(self.DB_select,self.DB_button,self.basis_select3,VBox(VBox(HBox(self.modes_select,HBox(self.desinc_mode,self.inc_mode,height=40))),self.draw_mode)))
        self.tab1 = Panel(child=HBox(self.inputs,VBox(self.p,self.psum),self.plog), title="Commands")
        self.tab2 = Panel(child=HBox(VBox(HBox(self.inputs2,self.p2,self.p3),self.cov_table,self.cor_table)), title="Covmats")
        self.tab3 = Panel(child=HBox(self.inputs3,self.pmodes), title="Basis")
        self.tabs = Tabs(tabs=[self.tab1,self.tab2,self.tab3])
        curdoc().clear()
        self.update(None,None,None)

        curdoc().add_root(self.tabs)#hplot(inputs,p))#, p, p2)
        curdoc().add_root(self.dialog)


    ######################################################################################
    #   ___      _ _ _             _
    #  / __|__ _| | | |__  __ _ __| |__ ___
    # | (__/ _` | | | '_ \/ _` / _| / /(_-<
    #  \___\__,_|_|_|_.__/\__,_\__|_\_\/__/
    #
    ######################################################################################
    def loadDB(self):
        self.dialog.visible=False
        self.dialog.content="Loading database..."
        self.dialog.visible = True

        self.f = h5py.File(self.datapath + str(self.DB_select.value))
        if(self.f.attrs.keys().count("target.Lambda")):
            self.Lambda_tar = self.f.attrs["target.Lambda"][0]
        else:
            self.Lambda_tar = 1.65
        if(self.f.keys().count("IF")): #Dense case
            self.IF = self.f["IF"][:]
        else: #Sparse case
            self.IF = csr_matrix((self.f["IF.data"][:],self.f["IF.indices"][:],self.f["IF.indptr"][:]))
            self.IF = self.IF.T
        self.P = self.f["P"][:]
        self.Btt = self.f["Btt"][:]
        self.modes = self.IF.dot(self.Btt)
        self.swap = np.arange(self.modes.shape[1])-2
        self.swap[0:2] = [self.modes.shape[1]-2,self.modes.shape[1]-1]
        #self.modes = self.modes[:,self.swap]
        self.indx_pup = self.f["indx_pup"][:]
        self.pup = np.zeros((self.f["dm_dim"].value,self.f["dm_dim"].value))
        self.update(None,None,None)
        self.niter = self.f["com"][:].shape[1]
        self.nactus = self.f["com"][:].shape[0]
        self.nmodes = self.P.shape[0]
        if(self.f.keys().count("cov")):
            self.cov = self.f["cov"][:]
            self.cor = self.f["cor"][:]
        else:
            self.cov, self.cor = self.cov_cor()

        #self.cov_table, self.cor_table = self.createDataTables()
        self.updateDataTables()


        print "DB loaded"
        self.dialog.visible=False

    def update(self,attrname,old,new):
       # plot_val = plot_type.value
        self.source1.data = dict(x=[],y=[],color=[],typec=[])

        coms_active = self.coms.active
        plot_val = self.plot_select.value
        basis_val = self.basis_select1.value
        iteration = int(self.iter_select.value)

        yi=[]
        xi=[]
        typec=[]
        coloris=[]
        for jj in coms_active:
            j = self.coms_list[jj]
            data=self.f[j][:]
            if(plot_val == "Commands"):
                if(basis_val == "Actuators"):
                    yi.append(data[:,iteration].tolist())
                    self.p.xaxis.axis_label = "Actuators"
                elif(basis_val == "Btt"):
                    yi.append(np.dot(self.P,data[:,iteration])[self.swap].tolist())
                    self.p.xaxis.axis_label = "Modes"
                xi.append(range(len(data[:,iteration])))
                typec.append([j]*len(data[:,iteration]))
                coloris.append(self.colors[j])
                self.p.yaxis.axis_label = "Volts"

            elif(plot_val == "Variance"):
                if(basis_val == "Actuators"):
                    yi.append(np.var(data,axis=1).tolist())
                    self.p.xaxis.axis_label = "Actuators"
                elif(basis_val == "Btt"):
                    yi.append(np.var(np.dot(self.P,data),axis=1)[self.swap].tolist())
                    self.p.xaxis.axis_label = "Modes"
                xi.append(range(len(np.var(data,axis=1))))
                typec.append([j]*len(np.var(data,axis=1)))
                coloris.append(self.colors[j])
                self.p.yaxis.axis_label = "Variance"

        self.source1.data = dict(x=xi,y=yi,color=coloris,typec=typec)

        print "Plots updated"


    def rescale_matrix(self):
        self.dialog.visible = False
        vmin = self.cmin.value
        vmax = self.cmax.value
        self.dialog.content="Updating matrix..."
        self.dialog.visible = True
        if(self.old):
            os.remove(self.old)
        time = str(datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%f'))
        self.old = "/home/fferreira/public_html/breakdown_display"+time+".png"
        mpl.image.imsave(self.old,self.covmat,vmin=vmin,vmax=vmax)
        self.p2.image_url(url=dict(value=self.url+time+".png"), x=0, y=0,w=self.covmat.shape[0],h=self.covmat.shape[0])
        self.dialog.visible = False

    def get_diag(self):
        x=np.arange(self.covmat.shape[0])
        y=np.diag(self.covmat)
        self.source2.data = dict(x=x,y=y)

    def cut_matrix(self):
        XorY=self.XY.labels[self.XY.active]
        ax = self.axiscut.value
        if(XorY == "X"):
            data = self.covmat[ax,:]
        else:
            data = self.covmat[:,ax]
        x=np.arange(data.size)
        self.source2.data = dict(x=x,y=data)

    def update_matrix2(self):
        self.dialog.visible = False
        if(self.old):
            os.remove(self.old)
        #self.draw.disabled = True
        A_val = self.A.value
        B_val = self.B.value
        basis = self.basis_select2.value
        powa = self.power.value
        self.dialog.content="Computing and loading matrix..."
        self.dialog.visible = True

        A_cov = self.f[A_val][:]
        B_cov = self.f[B_val][:]
        A_cov -= np.tile(np.mean(A_cov,axis=1),(A_cov.shape[1],1)).T
        B_cov -= np.tile(np.mean(B_cov,axis=1),(B_cov.shape[1],1)).T
        if(basis == "Btt"):
            A_cov = np.dot(self.P,A_cov)
            B_cov = np.dot(self.P,B_cov)
        print "Values ok"
        self.covmat = (np.dot(A_cov,B_cov.T)/B_cov.shape[1])
        print "dot product ok"
        if(powa != 1):
            self.covmat = np.abs(self.covmat)**powa * np.sign(self.covmat)
            print "scale adjusted"
        self.cmin.start = self.covmat.min()
        self.cmin.end = self.covmat.max()
        self.cmin.value = self.cmin.start
        self.cmin.step = (self.cmin.end - self.cmin.start)/100.
        self.cmax.start = self.covmat.min()
        self.cmax.end = self.covmat.max()
        self.cmax.value = self.cmax.end
        self.cmax.step = self.cmin.step
        self.axiscut.end = self.covmat.shape[0]
        time = str(datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%f'))
        self.old = "/home/fferreira/public_html/breakdown_display"+time+".png"
        mpl.image.imsave(self.old,self.covmat)
        self.p2.image_url(url=dict(value=self.url+time+".png"), x=0, y=0,w=self.covmat.shape[0],h=self.covmat.shape[0])

        #self.sourceC.data = dict(url=[self.url],x=0,y=covmat.shape[0],dw=covmat.shape[0],dh=covmat.shape[0])
        #self.draw.disabled = False
        print "Matrix updated2"
        self.dialog.visible = False

    def update_mode(self):
        self.dialog.visible = False
        if(self.old):
            os.remove(self.old)
        N = self.modes_select.value
        if(N>=self.modes.shape[1]):
            N=self.modes.shape[1]-1
            self.modes_select.value = N
        basis = self.basis_select3.value
        self.dialog.content="Loading..."
        self.dialog.visible = True

        if(basis == "Actuators"):
            pup = self.pup.flatten()
            pup[self.indx_pup] = self.IF[:,N].toarray()#self.f["IF"][:][:,N]
            self.pup = pup.reshape(self.pup.shape)
        elif(basis == "Btt"):
            pup = self.pup.flatten()
            pup[self.indx_pup] = self.modes[:,N-2]
            self.pup = pup.reshape(self.pup.shape)
        time = str(datetime.datetime.now().strftime('%Y-%m-%d_%H_%M_%f'))
        self.old = "/home/fferreira/public_html/breakdown_display"+time+".png"
        mpl.image.imsave(self.old,self.pup)
        self.pmodes.image_url(url=dict(value=self.url+time+".png"), x=0, y=0,w=self.pup.shape[0],h=self.pup.shape[0])

        #self.sourceC.data = dict(url=[self.url],x=0,y=covmat.shape[0],dw=covmat.shape[0],dh=covmat.shape[0])
        #self.draw.disabled = False
        print "Mode updated"
        self.dialog.visible = False

    def mode_increment(self):
        if(self.modes_select.value < self.modes.shape[1]-1):
            self.modes_select.value = self.modes_select.value+1
        else:
            self.modes_select.value = self.modes.shape[1]-1
    def mode_desincrement(self):
        if(self.modes_select.value > 0):
            self.modes_select.value = self.modes_select.value-1
        else:
            self.modes_select.value = 0

    def plot_sum(self):

        self.dialog.visible = False
        self.dialog.content="Computing..."
        self.dialog.visible = True

        plus = self.plus_select.active
        moins = self.moins_select.active
        basis_val = self.basis_select1.value
        plot_val = self.plot_select.value
        iteration = int(self.iter_select.value)

        if(plot_val == "Commands"):
            data = np.zeros(self.nactus)
            x = range(self.nactus)
        elif(plot_val == "Variance"):
            data = np.zeros((self.nmodes,self.niter))#self.nmodes)
            data2 = np.zeros(self.nmodes)
            x = range(self.nmodes)
        fitp = False
        fitm = False
        for i in plus:
            if(self.plus_select.labels[i] == "fitting"):
                fitp=True
            else:
                if(plot_val == "Commands"):
                    data += np.dot(self.P,self.f[self.coms_list[i]][:][:,iteration])
                elif(plot_val == "Variance"):
                    data += np.dot(self.P,self.f[self.coms_list[i]][:])
                    data2 += np.var(np.dot(self.P,self.f[self.coms_list[i]][:]),axis=1)
        for i in moins:
            if(self.plus_select.labels[i] == "fitting"):
                fitm=True
            else:
                if(plot_val == "Commands"):
                    data -= np.dot(self.P,self.f[self.coms_list[i]][:][:,iteration])
                elif(plot_val == "Variance"):
                    data -= np.dot(self.P,self.f[self.coms_list[i]][:])
                    data2 -= np.var(np.dot(self.P,self.f[self.coms_list[i]][:]),axis=1)
#        if(basis_val == "Btt"):
#            data = np.dot(self.P,data)
#            data2 = np.dot(self.P,data2)
        if(plot_val == "Variance"):
            data = np.var(data,axis=1)
            data = np.cumsum(data[self.swap])
            data2 = np.cumsum(data2[self.swap])
            data2 = np.exp(-data2*(2*np.pi/self.Lambda_tar)**2)
            data = np.exp(-data*(2*np.pi/self.Lambda_tar)**2)
            if(fitp and self.f.keys().count("fit_error")):
                data *= np.exp(-self.f["fit_error"].value)
                data2 *= np.exp(-self.f["fit_error"].value)
            if(fitm and self.f.keys().count("fit_error")):
                data /= np.exp(-self.f["fit_error"].value)
                data2 /= np.exp(-self.f["fit_error"].value)
        if(self.f.keys().count("SR2")):
            self.source3.data = dict(x=[x,x,x,x],y=[data,np.ones(len(x))*self.f["SR"].value,np.ones(len(x))*self.f["SR2"].value,data2],color=["blue","red","purple","green"])
        else:
            if(self.f.keys().count("SR")):
                self.source3.data = dict(x=[x,x,x],y=[data,np.ones(len(x))*self.f["SR"].value,data2],color=["blue","red","green"])
            else:
                self.source3.data = dict(x=x,y=data)
        print "Sum plotted"
        self.dialog.visible = False

    def cov_cor(self):
        cov = np.zeros((6,6))
        bufdict = {"0":self.f["noise_com"][:],
                   "1":self.f["trunc_com"][:],
                   "2":self.f["alias_wfs_com"][:],
                   "3":self.f["H_com"][:],
                   "4":self.f["bp_com"][:],
                   "5":self.f["tomo_com"][:]}
        for i in range(cov.shape[0]):
            for j in range(cov.shape[1]):
                if(j>=i):
                    tmpi = self.P.dot(bufdict[str(i)])
                    tmpj = self.P.dot(bufdict[str(j)])
                    cov[i,j] = np.sum(np.mean(tmpi*tmpj,axis=1)-np.mean(tmpi,axis=1)*np.mean(tmpj,axis=1))
                else:
                    cov[i,j] = cov[j,i]

        s = np.reshape(np.diag(cov),(cov.shape[0],1))
        sst = np.dot(s,s.T)
        cor = cov/np.sqrt(sst)

        return cov, cor

    def createDataTables(self):

        tmp = [TableColumn(field="Type", title="Covariance")]
        for item in self.DataTableItems[1:]:
            tmp.append(TableColumn(field=item, title=item))
        columns = tmp

        cov_table = DataTable(source=self.table_cov_source, columns=columns, width=1200, height=280)
        tmp[0] = TableColumn(field="Type", title="Correlation")
        cor_table = DataTable(source=self.table_cor_source, columns=columns, width=1200, height=280)

        return cov_table, cor_table

    def updateDataTables(self):
        self.table_cov_source.data = dict(Type=self.DataTableItems[1:],
                                           Noise=self.cov[:,0],
                                            Truncature=self.cov[:,1],
                                            Aliasing=self.cov[:,2],
                                            FilteredModes=self.cov[:,3],
                                            Bandwidth=self.cov[:,4],
                                            Tomography=self.cov[:,5])
        self.table_cor_source.data = dict(Type=self.DataTableItems[1:],
                                          Noise=self.cor[:,0],
                                            Truncature=self.cor[:,1],
                                            Aliasing=self.cor[:,2],
                                            FilteredModes=self.cor[:,3],
                                            Bandwidth=self.cor[:,4],
                                            Tomography=self.cor[:,5])




files = glob.glob("/home/fferreira/public_html/breakdown_display*")
for f in files:
    os.remove(f)

disp = html_display()



     # initial load of the data

    
