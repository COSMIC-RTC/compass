"""
Created on Tue Feb  2 09:39:35 2016

@author: fferreira
"""

import numpy as np

import h5py

from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.models.widgets import Select, Slider, CheckboxButtonGroup, Panel, Tabs
from bokeh.io import curdoc, vform, hplot
######################################################################################
#  _      _ _      
# (_)_ _ (_) |_ ___
# | | ' \| |  _(_-<
# |_|_||_|_|\__/__/
######################################################################################
                  
datapath = "/home/fferreira/Data/"

f = h5py.File(datapath+"breakdown_classical.h5")
f0 = h5py.File(datapath+"breakdown_imat_geom_meth0.h5")
f1 = h5py.File(datapath+"breakdown_imat_geom_meth1.h5")
fh = h5py.File(datapath+"breakdown_imat_hack.h5")

files = {"Classical":f,"Imat geom 0":f0,"Imat geom 1":f1,"Imat hacked":fh}
niter = f["com"][:].shape[1]
nactus = f["com"][:].shape[0]
meths = ["Classical","Imat geom 0","Imat geom 1","Imat hacked"]
plot_type = ["Commands","Variance"]
coms_list = f.keys()

######################################################################################
#         _    _          _      
# __ __ _(_)__| |__ _ ___| |_ ___
# \ V  V / / _` / _` / -_)  _(_-<
#  \_/\_/|_\__,_\__, \___|\__/__/
#               |___/            
######################################################################################
# Tab 1
coms = CheckboxButtonGroup(labels=coms_list,active=[0])
meth = CheckboxButtonGroup(labels=meths,active=[0])
plot_select = Select(title="Plot type",value=plot_type[1],options=plot_type)
iter_select = Slider(title="Iteration number",start=1,end=niter,step=1)
# Tab 2
A = Select(title="Commands A",value=coms_list[0],options=coms_list)
B = Select(title="Commands B",value=coms_list[0],options=coms_list)
imat = Select(title="Imat used",value=meths[0],options=meths)
power = Slider(title="Abs(covmat)**X",start=1,end=10,step=1)

colors = {"H_com":"green","alias_com":"blue","bp_com":"orange","com":"black","noise_com":"red","tomo_com":"purple","trunc_com":"cyan","wf_com":"magenta"}

source1 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
source2 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
source3 = ColumnDataSource(data=dict(x=[], y=[],color=[]))
source4 = ColumnDataSource(data=dict(x=[], y=[],color=[]))

sourceC = ColumnDataSource(data=dict(C=[]))

sources = {"Classical":source1,"Imat geom 0":source2,"Imat geom 1":source3,"Imat hacked":source4}
#hover = HoverTool(tooltips=[("Datas","@data_type")])
#
p=Figure(plot_height=600, plot_width=800)#,tools=[hover])
p.multi_line("x","y",color="color",source=source1)
p.multi_line("x","y",color="color",line_dash=[4,4],source=source2)
p.multi_line("x","y",color="color",line_dash=[2,2],source=source3)
p.multi_line("x","y",color="color",line_dash=[6,6],source=source4)
p.xaxis.axis_label = "Actuators"

p2=Figure(x_range=[0,nactus], y_range=[0,nactus])
#p2.image(image=[np.random.random((1000,1000))],x=[0],y=[0],dw=[1000],dh=[1000])#,source=sourceC)
#p.yaxis.axis_label = "Noise on commands [V]"

######################################################################################
#   ___      _ _ _             _       
#  / __|__ _| | | |__  __ _ __| |__ ___
# | (__/ _` | | | '_ \/ _` / _| / /(_-<
#  \___\__,_|_|_|_.__/\__,_\__|_\_\/__/
#                                      
######################################################################################
def update(attrname,old,new):
   # plot_val = plot_type.value
    source1.data = dict(x=[],y=[],color=[])
    source2.data = dict(x=[],y=[],color=[])
    source3.data = dict(x=[],y=[],color=[])
    source4.data = dict(x=[],y=[],color=[])

    coms_active = coms.active
    meth_active = meth.active
    plot_val = plot_select.value
    iteration = int(iter_select.value)
    A_val = A.value
    B_val = B.value
    ff = imat.value
    powa = np.ceil(power.value)

    for ii in meth_active:
        i = meths[ii]
        xi = []
        yi = []
        coloris = []
        for jj in coms_active:
            j = coms_list[jj]
            data=files[i][j][:]
            if(plot_val == "Commands"):
                yi.append(data[:,iteration].tolist())
                xi.append(range(len(data[:,iteration])))
                coloris.append(colors[j])
                p.yaxis.axis_label = "Volts"
            elif(plot_val == "Variance"):
                yi.append(np.var(data,axis=1).tolist())
                xi.append(range(len(np.var(data,axis=1))))
                coloris.append(colors[j])
                p.yaxis.axis_label = "Variance"
        source = sources[i]
        source.data = dict(x=xi,y=yi,color=coloris)
    
    A_cov = files[ff][A_val][:]
    B_cov = files[ff][B_val][:]
    covmat = np.dot(A_cov,B_cov.T)/B_cov.shape[1]
    if(powa != 1 and powa != 0):
        covmat = np.abs(covmat)**(1./powa)
        print "scale adjusted"
    sourceC.data = dict(C=covmat)
    #p2.image(image=[covmat],x=[0],y=[0],dw=[covmat.shape[0]],dh=[covmat.shape[0]])
    
    print "Update done"     

controls = [plot_select,iter_select,A,B,imat,power]
buttons = [coms,meth]
for control in controls:
    control.on_change('value', update)
for button in buttons:
    button.on_change('active',update)

inputs = hplot(vform(coms,meth,plot_select,iter_select))#, width=350)
inputs2 = hplot(vform(imat,A,B,power))#, width=350)
tab1 = Panel(child=hplot(inputs,p), title="Commands")
tab2 = Panel(child=hplot(inputs2,p2), title="Covmats")
tabs = Tabs(tabs=[tab1,tab2])

update(None,None,None) # initial load of the data

curdoc().add_root(tabs)#hplot(inputs,p))#, p, p2)

