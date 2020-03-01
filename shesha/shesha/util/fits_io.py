from astropy.io.fits import open as fitsopen
from astropy.io.fits import writeto as fitswrite
def fitsread(filename):
    return (fitsopen(filename)[0].data)

