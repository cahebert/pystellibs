import urllib.request
from astropy.io import fits
import numpy as np

remote_repo = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/'

params = []
spectra = []

# open catalog first and get the set of files (eg no logg variation)
with urllib.request.urlopen(remote_repo + 'catalog.fits') as f:
    catalog = fits.open(f)

files = list(set([d.split('[')[0] for d in catalog[1].data['FILENAME']]))

catalog.close()

for fname in files[:4]:
    # recover logZ, Teff, and logT
    temp = fname.split('/')[0][-3:]
    if temp[0] == 'p':
        logZ = float(temp[1:]) / 10
    elif temp[0] == 'm':
        logZ = -float(temp[1:]) / 10
    else:
        raise ValueError('Something went wrong recovering log Z!')

    Teff = int(fname.split('_')[-1].split('.')[0])
    logT = np.log10(Teff)

    with urllib.request.urlopen(remote_repo + fname) as f:
        hdulist = fits.open(f)

    if fname == files[0]:
        wavelengths = hdulist[1].data['WAVELENGTH']

    for logg_str in [c.name for c in hdulist[1].columns if 'g' in c.name][:2]:
        spectra.append(hdulist[1].data[logg_str])

        logg = float(logg_str.strip('g')) / 10
        params.append([logZ, logg, logT, Teff])

    hdulist.close()

params = np.array(np.array(params),
                  dtype=[('logZ', float),
                         ('logg', float),
                         ('logT', float),
                         ('Teff', float)])

spectra = np.array(spectra)
wavelengths = np.array(wavelengths)

# chuck params and spectra into an HDU
