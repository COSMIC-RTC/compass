
def write_loop(filename,loop,controller):
    """Write (append) AO loop parameters to file for YAO

    loop        : Param_loop        : compass loop parameters
    controller  : Param_controller  : compass controller parameters
    """
    f=open(filename,"a+")
    f.write("\n\n//------------------------------")
    f.write("\n//LOOP  parameters")
    f.write("\n//------------------------------")
    f.write("\nloop.method     = "+"\"none\""+";")  
    f.write("\nloop.leak       = "+str(0.001)+";")
    f.write("\nloop.gain       = "+str(controller.gain)+";")
    f.write("\nloop.framedelay = "+str(controller.delay)+";")
    f.write("\nloop.niter      = "+str(loop.niter)+";")
    f.write("\nloop.ittime     = "+str(loop.ittime)+";")
    f.write("\nloop.skipevery  = "+str(100000)+";")
    f.write("\nloop.startskip  = "+str(30)+";")
    f.write("\nloop.skipby     = "+str(5000)+";")

    f.close()
