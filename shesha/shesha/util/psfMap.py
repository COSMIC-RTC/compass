import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import os
import shutil


class PSF_map:
    def __init__( self, SR=np.zeros((0,0)) , radius=0, wl=0, r0=None,frac=None,NGS=None, LGS=None,TS=None, AOtype="AO", sup=None, filename=None):
        """ SR_map class

        SR      : np.ndarray        : array of strehl (organised by position)
        radius  : float             : radius of the map (in arcsec)
        wl      : float             : wavelength of observation (in nm)
        NGS     : np.ndarray        : position (arcsec) of the NGS (1st line:x axis, 2nd:y axis)
        LGS     : np.ndarray        : position (arcsec) of the LGS (1st line:x axis, 2nd:y axis)
        TS      : np.ndarray        : position (arcsec) of the TS (1st line:x axis, 2nd:y axis)
        sup     : shesha.supervisor : shesha supervisor to retrieve all data from
        filename: str               : filename to read the psf map from

        warning:
        using 'sup' will take only the supervisor into consideration
        using 'filename' will overide previous values of the psf map (unless called with sup)
        
        """
        self.map=SR
        
        self._Rtar=radius
        self._Rngs=0
        self._Rlgs=0
        self._Rlgs=0
        if(NGS is not None):
            self.NGSx=NGS[0]
            self.NGSy=NGS[1]
            self._Rngs= max((self.NGSx.max(),self.NGSy.max()))
        else:
            self.NGSx=np.zeros((0))
            self.NGSy=np.zeros((0))

        if(LGS is not None):
            self.LGSx=LGS[0]
            self.LGSy=LGS[1]
            self._Rlgs= max(self.LGSx.max(),self.LGSy.max())
        else:
            self.LGSx=np.zeros((0))
            self.LGSy=np.zeros((0))

        if(TS is not None):
            self.TSx=TS[0]
            self.TSy=TS[1]
            self._Rts= max(self.TSx.max(),self.TSy.max())
        else:
            self._Rts=0
            self.TSx=np.zeros((0))
            self.TSy=np.zeros((0))
        
        self.wavelength=wl
        self.type=AOtype

        self.r0=r0
        self.frac=frac
        
        if(sup is not None):
            self.NGSx=np.array([t.xpos for t in sup.config.p_wfs_ngs[:-1]])
            self.NGSy=np.array([t.ypos for t in sup.config.p_wfs_ngs[:-1]])
            if(self.NGSx.size>0):
                self._Rngs= max((self.NGSx.max(),self.NGSy.max()))
            self.LGSx=np.array([t.xpos for t in sup.config.p_wfs_lgs])
            self.LGSy=np.array([t.ypos for t in sup.config.p_wfs_lgs])
            if( self.LGSx.size>0):
                self._Rlgs= max(self.LGSx.max(),self.LGSy.max())
            self.TSx=np.array([t.xpos for t in sup.config.p_wfs_ts])
            self.TSy=np.array([t.ypos for t in sup.config.p_wfs_ts])
            if(TS is not None):
                self._Rts= max(self.TSx.max(),self.TSy.max())
            self.wavelength=sup.config.p_targets[0].Lambda
            NTAR=len(sup.config.p_targets)
            NTAR_side=int(np.sqrt(NTAR))
            if( NTAR!= NTAR_side**2):
                raise ValueError("not a square nb of targets")
            self.map=np.zeros((NTAR_side,NTAR_side))
            for i in range(NTAR):
                self.map.itemset(i,sup._sim.getStrehl(i)[1])
                #self.map.itemset(i,sup._sim.getStrehl(i)[0])
                tar=sup._sim.tar.d_targets[i]
                self._Rtar=max(self._Rtar,tar.posx,tar.posy)
            self.r0=sup.config.p_atmos.r0
            self.frac=sup.config.p_atmos.frac


        elif(filename is not None):
            self.read(filename)


    def setNGS(self,NGS,NGSy=None):
        if(NGSy is None):
            self.NGSx=NGS[0]
            self.NGSy=NGS[1]
        else:
            self.NGS=NGS
            self.NGSy=NGS
        self._Rngs= max((self.NGSx.max(),self.NGSy.max()))


    def setLGS(self,LGS,LGSy=None):
        if(LGSy is None):
            self.LGSx=LGS[0]
            self.LGSy=LGS[1]
        else:
            self.LGS=LGS
            self.LGSy=LGS
        self._Rlgs= max(self.LGSx.max(),self.LGSy.max())

    def setTS(self,TS,TSy=None):
        if(TSy is None):
            self.TSx=LGS[0]
            self.TSy=LGS[1]
        else:
            self.TS=LGS
            self.TSy=LGS
        self._Rts= max(self.TSx.max(),self.TSy.max())

    def setWaveLength(self,wl):
        self.wavelength=wl

    def plot(self, title=False,GS=False,WFS=False,LGS=False,NGS=False,TS=False,colorbar=True,vMin=None,vMax=None,xMax=0,yMax=0,axes=None):
        if(self.map.shape[0]>0 and self.map.shape[1]>0):
            xMin=-max(xMax,self._Rtar)
            xMax=max(xMax,self._Rtar)
            yMin=-max(yMax,self._Rtar)
            yMax=max(yMax,self._Rtar)
            if(vMax is not None):
                vMax=min(1,vMax)
            else:
                vMax=self.map.max()
            if(vMin is not None):
                vMin=max(0,vMin)
            else:
                vMin=self.map.min()
            if(axes is None):
                im = plt.imshow(self.map,extent=[-self._Rtar,self._Rtar,-self._Rtar,self._Rtar],vmin=vMin,vmax=vMax)
                axes=plt.gca()
            else:
                im = axes.imshow(self.map,extent=[-self._Rtar,self._Rtar,-self._Rtar,self._Rtar],vmin=vMin,vmax=vMax)
            if(colorbar):
                plt.colorbar()
            if(GS or WFS or LGS):
                axes.scatter(self.LGSy,self.LGSx,color="red")
                xMin=min(xMin,self.LGSx.min())
                xMax=max(xMax,self.LGSx.max())
                yMin=min(yMin,self.LGSy.min())
                yMax=max(yMax,self.LGSy.max())
            if(GS or WFS or NGS):
                axes.scatter(self.NGSy,self.NGSx,color="blue")
                xMin=min(xMin,self.NGSx.min())
                xMax=max(xMax,self.NGSx.max())
                yMin=min(yMin,self.NGSy.min())
                yMax=max(yMax,self.NGSy.max())
            if((GS or TS) and self.TSx.size>0):
                axes.scatter(self.TSy,self.TSx,color="black",s=4)
                xMin=min(xMin,self.TSx.min())
                xMax=max(xMax,self.TSx.max())
                yMin=min(yMin,self.TSy.min())
                yMax=max(yMax,self.TSy.max())

            t=self.type+" Strehl"
            if(self.wavelength>0):
                t+=" @ {:.3f} nm".format(self.wavelength)
            
            if(self._Rts>0):
                t+=", {:.2f}".format(self._Rts)+" arcsec ring optim"
            if(title):
                plt.title(t)
            if(axes is None):
                ax=plt.gca()
                ax.set_xlim(xMin,xMax)
                ax.set_ylim(yMin,yMax)
            else:
                axes.set_xlim(xMin,xMax)
                axes.set_ylim(yMin,yMax)

            return im


    def save(self,name=""):
        hdu_map=fits.PrimaryHDU(self.map)
        hdu_map.header["TYPE"]  =self.type
        hdu_map.header["LAMBDA"]=self.wavelength
        hdu_map.header["RTAR"]  =self._Rtar
        hdu_map.header["RLGS"]  =self._Rlgs
        hdu_map.header["RNGS"]  =self._Rngs
        hdu_map.header["RTS"]   =self._Rts
        hdu_map.header["R0"]    =self.r0

        hdu_LGSX   = fits.ImageHDU(self.LGSx, name="LGSX")
        hdu_LGSY   = fits.ImageHDU(self.LGSy, name="LGSY")
        hdu_NGSX   = fits.ImageHDU(self.NGSx, name="NGSX")
        hdu_NGSY   = fits.ImageHDU(self.NGSy, name="NGSY")
        hdu_TSX    = fits.ImageHDU(self.TSx,  name="TSX")
        hdu_TSY    = fits.ImageHDU(self.TSy,  name="TSY")
        hdu_FRAC   = fits.ImageHDU(self.frac, name="FRAC")

        hdul=fits.HDUList([hdu_map,
                            hdu_LGSX,
                            hdu_LGSY,
                            hdu_NGSX,
                            hdu_NGSY,
                            hdu_TSX ,
                            hdu_TSY,
                            hdu_FRAC
                            ])
        t=self.type+"_StrehlMap"
        if(self.wavelength>0):
            t+="_{:.3f}nm".format(self.wavelength)
        if(self._Rts>0):
            t+="_{:.2f}arcsec".format(self._Rts)
        t+=".fits"
        if(name==""):
            name=t
        hdul.writeto(name,overwrite=1)

    def read(self,name):
        hdu_map        = fits.open(name)
        self.type      = hdu_map[0].header["TYPE"]
        self.wavelength= hdu_map[0].header["LAMBDA"]
        self._Rtar     = hdu_map[0].header["RTAR"]
        self._Rlgs     = hdu_map[0].header["RLGS"]
        self._Rngs     = hdu_map[0].header["RNGS"]
        self._Rts      = hdu_map[0].header["RTS"]
        self.r0        = hdu_map[0].header["R0"]
        self.map       = hdu_map[0].data
        self.LGSx      = hdu_map["LGSx"].data 
        self.LGSy      = hdu_map["LGSy"].data 
        self.NGSx      = hdu_map["NGSx"].data 
        self.NGSy      = hdu_map["NGSy"].data 
        self.TSx       = hdu_map["TSx"].data 
        self.TSy       = hdu_map["TSy"].data  
        self.frac      = hdu_map["FRAC"].data  




def drawMapMOAO(MAPDIR,FRACS,R0S,TS,WFS,xMax=0,yMax=0,baseName="map_RTAR*",savePath=None, caption=None): 
     nRows=len(WFS)
     nCols=len(TS)
     maps=[]
     cc=0
     ccc=0
     print(nRows,nCols)
     for frac in FRACS:
         for R in R0S:
             fig,axes=plt.subplots(nrows=nRows,ncols=nCols)
             title_str="r0:"+str(R) +", fraction : "+str(frac)
             #fig.suptitle(title_str)   
             atmSuffix="r0-"+str(R)   
             for k,f in enumerate(frac):    
                 atmSuffix+="_l"+str(k)+"-"+str(f)    
             maps=[]  
             vMin=1  
             vMax=0  
             for i,W in enumerate(WFS):   
                 for j,T in enumerate(TS):   
                     cc+=1   
                     fileName=baseName+"_RWFS"+str(W)+"_RTS"+str(T[0])+"_NTS"+str(T[1])+"_"+atmSuffix+".fits"   
                     maps.append(PSF_map())   
                     if os.path.exists(MAPDIR+"/"+fileName):   
                         ccc+=1   
                         maps[-1].read(MAPDIR+"/"+fileName)   
                         vMin=min(vMin,maps[-1].map.min())  
                         vMax=max(vMax,maps[-1].map.max())  
                     else :   
                         print(fileName)   
             im=[] 
             for i in range(len(WFS)):  
                 for j in range(len(TS)):  
                    k=i*nCols+j 
                    im.append(maps[k].plot(WFS=True,TS=True,colorbar=False,vMin=vMin,vMax=vMax,xMax=xMax,yMax=yMax,axes=axes[i][j])) 
             fig.subplots_adjust(right=0.85)  
             cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])  
             fig.colorbar(im[-1], cax=cbar_ax) 
             fig.set_size_inches(12.,8.31)
             #fig.set_size_inches(8.25,7.5)
             if(caption is not None):
                cap=caption+"\nAtmosphere : r0:"+str(R)
                cap+="\ncomposed of "+str(len(frac))+ " layers"
                cap+="\n    altitude (m)   : [0 1000 4500 9000]"
                cap+="\n    r0 fraction (%): "+np.array2string(frac)
                fig.text(0,0,cap)
             if(savePath is not None):
                fig.savefig(savePath+atmSuffix+".png",bbox_inches='tight')
     print(cc,ccc)


def drawMapLTAO(MAPDIR,FRACS,R0S,WFS,xMax=0,yMax=0,baseName="map_RTAR15.0",savePath=None, caption=None): 
    nRows=len(WFS)
    nCols=1 # len(TS)
    maps=[]
    atm_conf=[]
    cc=0    
    ccc=0    
    vMin=1
    vMax=0 
    for frac in FRACS: 
        for R in R0S: 
            fig,axes=plt.subplots(nrows=nRows,ncols=nCols)
            title_str="r0:"+str(R) +", fraction : "+str(frac)
            #fig.suptitle(title_str)   
            atmSuffix="r0-"+str(R)   
            for k,f in enumerate(frac):    
                atmSuffix+="_l"+str(k)+"-"+str(f)    
            maps=[]  
            vMin=1  
            vMax=0  
            for i,W in enumerate(WFS):   
                cc+=1   
                fileName=baseName+"_RWFS"+str(W)+"_"+atmSuffix+".fits"   
                maps.append(PSF_map())   
                if os.path.exists(MAPDIR+"/"+fileName):   
                    ccc+=1   
                    maps[-1].read(MAPDIR+"/"+fileName)   
                    vMin=min(vMin,maps[-1].map.min())  
                    vMax=max(vMax,maps[-1].map.max())  
                else :   
                    print(fileName)   
            im=[] 
            for i in range(len(WFS)):  
                im.append(maps[i].plot(WFS=True,colorbar=False,vMin=vMin,vMax=vMax,xMax=xMax,yMax=yMax,axes=axes[i])) 
            fig.subplots_adjust(right=0.7) 
            cbar_ax = fig.add_axes([0.7, 0.15, 0.05, 0.7]) 
            fig.colorbar(im[-1], cax=cbar_ax) 
            fig.set_size_inches(5.,10.0) 
            #fig.set_size_inches(8.25,7.5)
            if(caption is not None):
                cap=caption+"\nAtmosphere : r0:"+str(R)
                cap+="\ncomposed of "+str(len(frac))+ " layers"
                cap+="\n    altitude (m)   : [0 1000 4500 9000]"
                cap+="\n    r0 fraction (%): "+np.array2string(frac)
                fig.text(0,0,cap)

            if(savePath is not None): 
                fig.savefig(savePath+atmSuffix+".png",bbox_inches='tight')
            print(cc,ccc) 
                 
