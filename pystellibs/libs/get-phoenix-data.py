import urllib.request
from astropy.io import fits
import numpy as np
import astropy.table

remote_repo = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/'
ZSUN = .0134

params = []
spectra = []

# open catalog first and get the set of files (eg no logg variation)
with urllib.request.urlopen(remote_repo + 'catalog.fits') as f:
    catalog = fits.open(f)

files = list(set([d.split('[')[0] for d in catalog[1].data['FILENAME']]))

catalog.close()

ref_wavelengths = None

for fname in files:
    # recover logZ, Z, Teff, and logT
    temp = fname.split('/')[0][-3:]
    if temp[0] == 'p':
        logZ = float(temp[1:]) / 10
    elif temp[0] == 'm':
        logZ = -float(temp[1:]) / 10
    else:
        raise ValueError('Something went wrong recovering log Z!')
    Z = 10**(logZ + np.log10(ZSUN))

    Teff = int(fname.split('_')[-1].split('.')[0])
    logT = np.log10(Teff)

    with urllib.request.urlopen(remote_repo + fname) as f:
        hdulist = fits.open(f)

    wavelengths = hdulist[1].data['WAVELENGTH']
    
    # check whether wavelengths match the reference 
    # or do we prefer just going with a nearest neighbor in the grid instead of combining?
    if len(wavelengths) != 5000:
        if ref_wavelengths is None:
            raise ValueError('Reference wavelengths not yet defined!')
    elif ref_wavelengths is None:
        ref_wavelengths = wavelengths

    for logg_str in [c.name for c in hdulist[1].columns if 'g' in c.name]:
        logg = float(logg_str.strip('g')) / 10
        params.append([logZ, logg, logT, Teff, Z])

        spec = hdulist[1].data[logg_str]
        # resample the spectra to the reference so that the lengths all match
        if len(wavelengths) != 5000:
            resampled_spec = np.interp(ref_wavelengths, wavelengths, spec)
            spectra.append(resampled_spec)
        else:
            spectra.append(spec)

    hdulist.close()

params = astropy.table.Table(np.array(params, dtype=np.float32),
                             names=['logZ', 'logg', 'logT', 'Teff', 'Z'])

# wavelengths are the last row of the array
spectra = np.vstack([spectra, wavelengths])

hdul = fits.HDUList()
hdul.append(fits.PrimaryHDU(data=spectra))
hdul.append(fits.BinTableHDU(params, name='PHOENIX'))
hdul.writeto('phoenix.grid.fits', overwrite=True)
