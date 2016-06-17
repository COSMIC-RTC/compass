import numpy as np
import os
import scipy.ndimage.interpolation as interp
import hdf5_utils as h5u


#import matplotlib.pyplot as pl

EELT_data=os.environ.get('SHESHA_ROOT')+"/data/apertures/"

def make_pupil(dim,pupd,tel,xc=-1,
                yc=-1,real=0,cobs=-1):

    #TODO ohter types
    if(tel.type_ap=="EELT-Nominal"):
        tel.set_cobs(0.3)
        N_seg=798
        return make_EELT(dim,pupd,tel,N_seg)
    elif(tel.type_ap=="EELT-BP1"):
        tel.set_cobs(0.369)
        N_seg=768
        return make_EELT(dim,pupd,tel,N_seg)
    elif(tel.type_ap=="EELT-BP3"):
        tel.set_cobs(0.503)
        N_seg=672
        return make_EELT(dim,pupd,tel,N_seg)
    elif(tel.type_ap=="EELT-BP5"):
        tel.set_cobs(0.632)
        N_seg=558
        return make_EELT(dim,pupd,tel,N_seg)
    elif(tel.type_ap=="VLT"):
        raise NotImplementedError("make_VLT")
    else:
        tel.set_type_ap("Generic")
        return make_pupil_generic(dim,pupd,tel.t_spiders,tel.spiders_type,
                                xc,yc,real,tel.cobs)




def make_pupil_generic(dim, pupd, t_spiders=0.01, spiders_type="six",
                    xc=0, yc=0, real=0, cobs=0):
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

                #angle = np.pi/(180/15.)
                angle = 0
                s2ma_2=2*np.sin(np.pi/2-angle)
                s6pa_2=2*np.sin(np.pi/6+angle)
                s6ma_2=2*np.sin(np.pi/6-angle)
                t2ma=np.tan(np.pi/2-angle)
                t6pa=np.tan(np.pi/6+angle)
                t6ma=np.tan(np.pi/6-angle)

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
    """TODO"""

    EELT_file=EELT_data+tel.type_ap+"_N"+str(dim)+"_COBS"+str(100*tel.cobs)+"_CLOCKED"+str(tel.pupangle)+"_TSPIDERS"+str(100*tel.t_spiders)+"_MS"+str(tel.nbrmissing)+"_REFERR"+str(100*tel.referr)+".h5"

    if( os.path.isfile(EELT_file) ):
        print "reading EELT pupil from file ", EELT_file
        pup=h5u.readHdf5SingleDataset(EELT_file)
    else:
        print "creating EELT pupil..."
        file= EELT_data+"Coord_"+tel.type_ap+".dat"
        data=np.fromfile(file,sep="\n")
        data=np.reshape(data,(data.size/2,2))
        x_seg=data[:,0]
        y_seg=data[:,1]


        file= EELT_data+"EELT_MISSING_"+tel.type_ap+".dat"
        k_seg=np.fromfile(file,sep="\n").astype(np.int32)

        W=1.45*np.cos(np.pi/6)

        #tel.set_diam(39.146)
        tel.set_diam(37.)

        X=MESH(tel.diam*dim/pupd,dim)

        #t_spiders=0.014
        t_spiders=0.06
        tel.set_t_spiders(t_spiders)

        if(tel.nbrmissing>0):
            k_seg=np.sort(k_seg[:tel.nbrmissing])


        file= EELT_data+"EELT_REF_ERROR"+".dat"
        ref_err=np.fromfile(file,sep="\n")

        #mean_ref = np.sum(ref_err)/798.
        #std_ref = np.sqrt(1./798.*np.sum((ref_err-mean_ref)**2))
        #mean_ref=np.mean(ref_err)
        std_ref=np.std(ref_err)

        ref_err = ref_err * tel.referr/ std_ref

        if(tel.nbrmissing > 0):
            ref_err[k_seg]=1.0

        pup=np.zeros((dim,dim))

        t_3=np.tan(np.pi/3.)

        for i in range(N_seg):
            Xt=X+x_seg[i]
            Yt=X.T+y_seg[i]
            pup+=(1-ref_err[i])*(Yt<0.5*W)*(Yt>=-0.5*W)*(0.5*(Yt+t_3*Xt)<0.5*W) \
                               *(0.5*(Yt+t_3*Xt)>=-0.5*W)*(0.5*(Yt-t_3*Xt)<0.5*W) \
                               *(0.5*(Yt-t_3*Xt)>=-0.5*W)

        t_spiders= tel.t_spiders*tel.diam*dim/pupd

        s2_6=2*np.sin(np.pi/6)
        t_6 =np.tan(np.pi/6)

        spiders_map = np.abs(X) > t_spiders/2
        spiders_map*= ((X.T>(X+t_spiders/s2_6)*t_6) + (X.T<(X-t_spiders/s2_6)*t_6))
        spiders_map*= ((X.T>(-X+t_spiders/s2_6)*t_6)+(X.T<(-X-t_spiders/s2_6)*t_6))

        pup = pup*spiders_map

        if (tel.pupangle != 0):
            pup=interp.rotate(pup,tel.pupangle,reshape=False,order=2)


        print "writing EELT pupil to file ", EELT_file
        h5u.writeHdf5SingleDataset(EELT_file,pup)

    print "EELT pupil created"
    return pup

def make_phase_ab(dim,pupd,tel,pup):
    """Compute the EELT M1 phase aberration
    :parameters:
	TODO
    """

    if(tel.type_ap=="Generic"):
	return np.zeros((dim,dim)).astype(np.float32)

    ab_file=EELT_data+"aberration_"+tel.type_ap+"_N"+str(dim)+"_NPUP"+str(np.where(pup)[0].size)+"_CLOCKED"+str(tel.pupangle)+"_TSPIDERS"+str(100*tel.t_spiders)+"_MS"+str(tel.nbrmissing)+"_REFERR"+str(100*tel.referr)+"_PIS"+str(tel.std_piston)+"_TT"+str(tel.std_tt)+".h5"
    if( os.path.isfile(ab_file) ):
        print "reading aberration phase from file ", ab_file
        phase_error=h5u.readHdf5SingleDataset(ab_file)
    else:
        print "computing M1 phase aberration..."

        std_piston=tel.std_piston
        std_tt=tel.std_tt

        W=1.45*np.cos(np.pi/6)

        file= EELT_data+"EELT_Piston_"+tel.type_ap+".dat"
        p_seg=np.fromfile(file,sep="\n")
        #mean_pis=np.mean(p_seg)
        std_pis=np.std(p_seg)
        p_seg=p_seg*std_piston/std_pis
        N_seg=p_seg.size

        file= EELT_data+"EELT_TT_"+tel.type_ap+".dat"
        tt_seg=np.fromfile(file,sep="\n")

        file= EELT_data+"EELT_TT_DIRECTION_"+tel.type_ap+".dat"
        tt_phi_seg=np.fromfile(file,sep="\n")

        phase_error=np.zeros((dim,dim))
        phase_tt=np.zeros((dim,dim))
        phase_defoc=np.zeros((dim,dim))

        file= EELT_data+"Coord_"+tel.type_ap+".dat"
        data=np.fromfile(file,sep="\n")
        data=np.reshape(data,(data.size/2,2))
        x_seg=data[:,0]
        y_seg=data[:,1]

        X=MESH(tel.diam*dim/pupd,dim)

        t_3=np.tan(np.pi/3.)

        for i in xrange(N_seg):
            Xt=X+x_seg[i]
            Yt=X.T+y_seg[i]
            SEG=(Yt<0.5*W)*(Yt>=-0.5*W)*(0.5*(Yt+t_3*Xt)<0.5*W) \
                               *(0.5*(Yt+t_3*Xt)>=-0.5*W)*(0.5*(Yt-t_3*Xt)<0.5*W) \
                               *(0.5*(Yt-t_3*Xt)>=-0.5*W)

            if(i==0):
                N_in_seg=np.sum(SEG)
                Hex_diam=2*np.max(np.sqrt(Xt[np.where(SEG)]**2+Yt[np.where(SEG)]**2))

            if(tt_seg[i]!=0):
                TT=tt_seg[i] * (np.cos(tt_phi_seg[i])*Xt+np.sin(tt_phi_seg[i])*Yt)
                mean_tt=np.sum(TT[np.where(SEG==1)])/N_in_seg
                phase_tt+=SEG*(TT-mean_tt)

            #TODO defocus

            phase_error += SEG*p_seg[i]


        N_EELT=np.where(pup)[0].size
        if(np.sum(phase_tt)!=0):
            phase_tt*=std_tt/np.sqrt(1./N_EELT*np.sum(phase_tt[np.where(pup)]**2))

        #TODO defocus

        phase_error+=phase_tt+phase_defoc

        if (tel.pupangle != 0):
            phase_error=interp.rotate(phase_error,tel.pupangle,reshape=False,order=2)

        print "phase aberration created"
        print "writing aberration filel to file ",ab_file
        h5u.writeHdf5SingleDataset(ab_file,phase_error)

    return phase_error

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
        xc=int(dim/2.)
    else:
        xc-=1.
    if(yc <0):
        yc=int(dim/2.)
    else:
        yc-=1.

    dx=np.tile(np.arange(dim)-xc,(dim,1))
    dy=np.tile(np.arange(dim)-yc,(dim,1)).T

    d=np.sqrt(dx**2+dy**2)
    return d #np.asfortranarray(d)

