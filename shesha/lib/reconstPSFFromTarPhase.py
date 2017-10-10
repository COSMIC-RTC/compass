
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pfits

lambdaTargetInMicrons = 2.2 # lambda used to generate PSF...
Dtel = 38.542
pup = pfits.getdata("pupTarget_29-08-2017_14:33:42_.fits") # Telescope pupil (small size support)
phase = pfits.getdata("TarPhase__29-08-2017_14:33:42_.fits")# data cube (niter, smallsizex, smallsizey)
phsize = pup.shape[0]

pixsizeInMeters = Dtel/phsize # pixel scale in meters
# here beware of the sampling for the PSF (make sure the PSF is well sampled compared to diffraction limit)
bigsize = 2048
pixsizeArcSec = (lambdaTargetInMicrons*1e-6)/(pixsizeInMeters*bigsize)*206265. # pixel scale of the PSF in arcseconds


pupBig = np.zeros((bigsize,bigsize))
phaseBig = pupBig*0.
bigsize = pupBig.shape[0]
niter = 0 # desired iteration number
phaseBig[(bigsize-phsize)/2:(bigsize+phsize)/2,  (bigsize-phsize)/2:(bigsize+phsize)/2] = phase[niter,:,:]
pupBig[(bigsize-phsize)/2:(bigsize+phsize)/2, (bigsize-phsize)/2:(bigsize+phsize)/2] = pup
PSFSE = np.fft.fftshift(np.abs( np.fft.fft2( pupBig*np.exp(1j*phaseBig*2*np.pi/lambdaTargetInMicrons) ))**2 / (np.sum(pupBig)**2))

plt.matshow(np.log(PSFSE)) # plot of the PSF in log scale...
plt.show()
print np.max(PSFSE) #Strehl ratio of the PSF...
