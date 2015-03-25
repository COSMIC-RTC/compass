import numpy as np

#import matplotlib.pyplot as pl

EELT_data="./EELT_data/"
def make_pupil(dim,pupd,tel,xc=-1,
                yc=-1,real=0,cobs=-1):
    
    #if type_ap=="Generic":
    return make_pupil_generic(dim,pupd,tel.t_spiders,tel.spiders_type,
                                xc,yc,real,tel.cobs)




def make_pupil_generic(dim, pupd, t_spiders=0.01, spiders_type="six",
                    xc=-1, yc=-1, real=0, cobs=0):
    """
    cdef np.ndarray pup = dist(dim,xc,yc)
    cdef np.ndarray[ndim=2,dtype=np.float32_t] spiders_map
    cdef float angle
    """


    pup=dist(dim,xc,yc)
    
    if(real==1):
        pup  = np.exp(-(pup /(pupd*0.5))**60.0)**0.69314
    else:
        pup = (pup<(pupd+1.)/2.).astype(np.float32)

    if (cobs>0):
        if(real==1):
            pup -= np.exp(-(dist(dim,xc,yc)/(pupd*cobs*0.5))**60.)**0.69314;
        else:
            pup -= (dist(dim,xc,yc) < (pupd*cobs+1.)*0.5).astype(np.float32)

            step=1./dim
            first=0.5*(step-1)

            X=np.tile(np.arange(dim,dtype=np.float32)*step+first,(dim,1))

            if(t_spiders<0):
                t_spiders=0.01
            t_spiders=t_spiders*pupd/dim
        
            if (spiders_type=="four"):

                s4_2=2*np.sin(np.pi/4)
                t4=np.tan(np.pi/4)

                spiders_map = ( (X.T > (X+t_spiders/s4_2)*t4) + (X.T < (X-t_spiders/s4_2)*t4) ).astype(np.float32)
                spiders_map*= ( (X.T > (-X+t_spiders/s4_2)*t4)+ (X.T < (-X-t_spiders/s4_2)*t4)).astype(np.float32)
        
                pup = pup*spiders_map


            elif (spiders_type=="six"):
      
                angle = 180/15.
                s2ma_2=2*np.sin(np.pi/2-np.pi/angle)
                s6pa_2=2*np.sin(np.pi/6+np.pi/angle)
                s6ma_2=2*np.sin(np.pi/6-np.pi/angle)
                t2ma=np.tan(np.pi/2-np.pi/angle)
                t6pa=np.tan(np.pi/6+np.pi/angle)
                t6ma=np.tan(np.pi/6-np.pi/angle)

                spiders_map = ((X.T > (-X+t_spiders/s2ma_2)*t2ma)+ (X.T < (-X-t_spiders/s2ma_2)*t2ma ))
                spiders_map*= ((X.T > ( X+t_spiders/s6pa_2)*t6pa)+ (X.T < ( X-t_spiders/s6pa_2)*t6pa ))
                spiders_map*= ((X.T > (-X+t_spiders/s6ma_2)*t6ma)+ (X.T < (-X-t_spiders/s6ma_2)*t6ma ))
                pup = pup*spiders_map

    print "generic pupil created"
    return pup


def make_VLT(dim,pupd,tel):
    tel.diam=8.0
    tel.cobs=0.14
    tel.t_spiders=0.005
    angle= 50.5*np.pi/180


    X=MESH(1.,dim);
    R=np.sqrt(X**2+(X.T)**2)

    pup=((R<0.5) & (R>(tel.cobs/2)) ).astype(np.float32)
    

    spiders_map= ((X.T>(X-tel.cobs/2+tel.t_spiders/np.sin(angle))*np.tan(angle))+ (X.T<(X-tel.cobs/2)*np.tan(tel.pupangle))) *(X>0)*(X.T>0)
    spiders_map+= np.fliplr(spiders_map)
    spiders_map+= np.flipud(spiders_map)

    pup = pup*spiders_map;
    
    print "VLT pupil created"
    return pup;



def make_EELT(dim,pupd,tel,N_seg):#dim,pupd,type_ap,cobs,N_seg,nbr_miss_seg,std_ref_err,t_spiders,angle)
    
    EELT_file=EELT_data+type_ap+"_N"+str(dim)+"_COBS"+str(100*cobs)+"_CLOCKED"+str(angle)+"_TSPIDERS"+str(100*t_spiders)+"_MS"+str(nbr_miss_seg)+"_REFERR"+str(100*std_ref_err)+".npy"


    if( os.path.isfile(EELT_file) ):
        pup=np.load(EELT_file)
    else:  
        file= EELT_data+"Coord_"+type_ap+".npz"
        data=np.load(file)
        x_seg=data['x_seg']
        y_seg=data['y_seg']

        file= EELT_data+"EELT_MISSING_"+type_ap+".npy"
        kseg=np.load(file)

        W=1.45*cos(pi/6)
        tel.diam=40.

        X=MESH(tel.diam*dim/pupd,dim)

        if(tel.nbrmissing>0):
                k_seg=np.sort(k_seg[:tel.nbrmissing])


        file= EELT_data+"EELT_REF_ERROR"+".npy"
        ref_err=np.load(file)
        mean_ref = np.sum(ref_err)/798

        std_ref = np.sqrt(1./798.*np.sum((ref_err-mean_ref)**2))
        ref_err = ref_err * tel.referr/ std_ref

        if(tel.nbrmissing > 0):
            ref_err[k_seg]=1.0

        pup=np.zeros((dim,dim))

        t_3=np.tan(np.pi/3)

        for i in range(N_seg):
            Xt=X+x_seg[i]
            Yt=X.T+y_seg[i]
            pup+=(1-ref_err[i])*(Yt<0.5*W)*(Yt>=-0.5*W)*(0.5*(Yt+t_3)*Xt<0.5*W) \
                               *(0.5*(Yt+t_3)*Xt>=-0.5*W)*(0.5*(Yt-t_3)*Xt<0.5*W) \
                               *(0.5*(Yt-t_3)*Xt>=-0.5*W)

        t_spiders= tel.t_spiders*tel.diam*dim/pupd

        s2_6=2*np.sin(np.pi/6)
        t_6 =np.tan(np.pi/6)

        spiders_map = np.abs(X) > t_spiders/2
        spiders_map*= ((X.T>(X+t_spiders/s2_6)*t_6) + (X.T<(X-t_spiders/s2_6)*t_6))
        spiders_map*= ((X.T>(-X+t_spiders/s2_6)*t_6)+(X.T<(-X-t_spiders/s2_6)*t_6))

        pup = pup*spiders_map

        #TODO rotate2
        """
        yoga_ao/yorick/yoga_ao_utils.i:658
        if (angle != 0) pup=rotate2(pup,angle);
        "EELT pupil has been created.";
        """


        np.save(EELT_file,pup)

    print "EELT pupil created"
    return pup




def MESH(Range,Dim):
    last=(0.5*Range-0.25/Dim)
    step= (2*last)/(Dim-1)
    
    return np.tile(np.arange(Dim)*step-last,(Dim,1))


def pad_array(A,N):
    S=A.shape
    D1=(N-S[0])/2
    D2=(N-S[1])/2
    padded=np.zeros((N,N))
    padded[D1:D1+S[0],D2:D2+S[1]]=A
    return padded




def dist(dim, xc=-1, yc=-1):
    if(xc <0):
        xc=(dim-1)/2.
    else:
        xc-=1.
    if(yc <0):
        yc=(dim-1)/2.
    else:
        yc-=1.

    dx=np.tile(np.arange(dim)-xc,(dim,1))
    dy=dx.T

    d=np.sqrt(dx**2+dy**2)
    return np.asfortranarray(d)

'''
def rotate2(image,angle, xc=-1,yc=-1, splin=0,outside=0):
    """rotate2(image,angle,xc,yc,splin,outside)

    Rotate the input image. Angle is in degrees, CCW.

    KEYWORDS:
    xc, yc: Center for coordinate transform. Note that this is
    compatible with the center defined by dist(), but is
    offset by 0.5 pixels w.r.t what you read on the yorick graphic
    window. I.e. the center of the bottom- left pixel is (1,1) in this
    function's conventions, not (0.5,0.5).

    splin: use spline2() instead of bilinear() for the interpolation

    outside: value for outliers.
    """

    angle *= np.pi/180.

    x,y = indices(image.shape[0],image.shape[1])

    if (xc<0) xc=np.ceil(image.shape[0]/2.+0.5)
    if (yc<0) yc=np.ceil(image.shape[1]/2.+0.5)

  x-=xc
  y-=yc

  x =  np.cos(angle)*x + np.sin(angle)*y
  y = -np.sin(angle)*x + np.cos(angle)*y

  x +=xc;
  y +=yc;

  if (splin!=0) return spline2(image,x,y,outside=outside)
  return bilinear(image,x,y,outside=outside)
'''



