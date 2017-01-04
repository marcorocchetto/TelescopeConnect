from library import *
import argparse
import pytz
from astropy.io import fits
import sys

#loading parameter file parser
parser = argparse.ArgumentParser()
parser.add_argument('--fits_fname',
                    dest='fits_fname',
                    type=str,
                    default=False,
                    )
parser.add_argument('--time_zone',
                    dest='time_zone',
                    type=str,
                    default='Europe/London',
                    )

options = parser.parse_args()

if not options.fits_fname:
    RaiseError('File name not specified')

if not options.time_zone:
    RaiseError('Time zone not specified')
else:
    try:
        local = pytz.timezone(options.time_zone)
    except pytz.exceptions.UnknownTimeZoneError:
        RaiseError('Unknown timezone')


# check file exists
if not os.path.isfile(options.fits_fname):
    RaiseError('File not found')

# open fits file hdu
try:
    hdulist = fits.open(options.fits_fname)
except:
    RaiseError('Unexpected error:', sys.exc_info()[0])

# get header & data
header = hdulist[0].header
data = hdulist[0].data

# check that some required FITS headers are present, otherwise raise error
header_keys = ['DATE-OBS', 'IMAGETYP', 'CCD-TEMP', 'EXPTIME']
for key in header_keys:
    if not key in hdulist[0].header:
        RaiseError('Header key `%s` missing' % key)

# determine Image type from header IMAGETYP
imagetype = get_ImageType(header['IMAGETYP'])
if not imagetype:
    RaiseError('Cannot determine ImageType')

# Date and time of observation: keyword DATE-OBS
dateobs_utc_str, dateobs_utc_datetime = get_DateObs(header['DATE-OBS'])


# Determine OBSNIGHT, i.e. night of observation
obsnight_str = get_ObsNight(dateobs_utc_datetime, local)


# RA and DEC, header keys: OBJCTRA, OBJCTDEC or RA, DEC
if 'OBJCTRA' in header:
    ra_str = header['OBJCTRA']
    ra = SexToDeg(ra_str, 'ra')
elif 'RA' in header:
    ra_str = header['RA']
    ra = SexToDeg(ra_str, 'ra')
else:
    ra_str = ''
    ra = None
if 'OBJCTDEC' in header:
    dec_str = header['OBJCTDEC']
    dec = SexToDeg(dec_str, 'dec')
elif 'DEC' in header:
    dec_str = header['OBJCTDEC']
    dec = SexToDeg(dec_str, 'dec')
else:
    dec_str = ''
    dec = None

# optional headers: TELESCOP, INSTRUME, OBSERVER, OBJECT
if 'TELESCOP' in header:
    telescope = header['TELESCOP']
else:
    telescope = ''
if 'INSTRUME' in header:
    instrument = header['INSTRUME']
else:
    instrument = ''
if 'OBSERVER' in header:
    observer = header['OBSERVER']
else:
    observer = ''
if 'OBJECT' in header:
    object = header['OBJECT']
else:
    object = ''

# exposure time
exptime = header['EXPTIME']

# filter
if 'FILTER' in header and imagetype != 'BIAS':
    filter = header['FILTER']
else:
    filter = ''

# ccd temperature
ccdtemp = header['CCD-TEMP']

# image size
width_px = len(data[0,:])
height_px = len(data[:,0])

# set filename
filename = dateobs_utc_datetime.strftime('%Y-%m-%dT%H-%M-%S')
if imagetype == 'LIGHT':
    if 'OBJECT' in header:
        filename += '-' + header['OBJECT']
    filename += header['FILTER']
if imagetype == 'BIAS':
    filename += '-Bias'
if imagetype == 'DARK':
    filename += '-Dark'
if imagetype == 'FLAT':
    filename += '-Flat'
    filename += header['FILTER']
filename += '.fits'
filename = filename.replace(' ', '_')

output = {
    'dateobs_utc_str': dateobs_utc_str,
    'obsnight_str': obsnight_str,
    'imagetype': imagetype,
    'object': object,
    'ra_float': ra,
    'dec_float': dec,
    'ra_str': ra_str,
    'dec_str': dec_str,
    'telescope': telescope,
    'instrument': instrument,
    'observer': observer,
    'exptime': exptime,
    'filter': filter,
    'width_px': width_px,
    'height_px': height_px,
    'ccdtemp': ccdtemp,
    'indexing_version': get_Version(),
    'filename': filename,
}

print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))