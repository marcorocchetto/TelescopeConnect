'''

Indexing of frames

The input to the script is a json file with the following format. See input_examples

Example:

python /home/telescope/TelescopeConnect/index_fits.py --json=/home/telescope/TelescopeConnect/input_examples/index_fits.json

'''


import argparse
import pytz
from astropy.io import fits
import sys
import json
import os

import library.general
from library.general import *

# #loading parameter file parser
# parser = argparse.ArgumentParser()
# parser.add_argument('--fits_fname',
#                     dest='fits_fname',
#                     type=str,
#                     default=False,
#                     )
# parser.add_argument('--time_zone',
#                     dest='time_zone',
#                     type=str,
#                     default='Europe/London',
#                     )
# 
# options = parser.parse_args()


#loading parameter file parser
parser = argparse.ArgumentParser()
parser.add_argument('--json',
                    dest='json_filename',
                    type=str,
                    default=False,
                    )

options = parser.parse_args()

if not options.json_filename:
    RaiseError('You need an inpunt JSON file')

json_data = json.loads(open(options.json_filename).read())


if not json_data['fits_fname']:
    RaiseError('File name not specified')

if not json_data['time_zone']:
    RaiseError('Time zone not specified')
else:
    try:
        local = pytz.timezone(json_data['time_zone'])
    except pytz.exceptions.UnknownTimeZoneError:
        RaiseError('Unknown timezone')


# check file exists
if not os.path.isfile(json_data['fits_fname']):
    RaiseError('File not found')

# open fits file hdu
try:
    hdulist = fits.open(json_data['fits_fname'])
except:
    RaiseError('Unexpected error:', sys.exc_info()[0])

# get header & data
header = hdulist[0].header
data = hdulist[0].data

# temporary fix for PinPoint simulated images
if not 'IMAGETYP' in hdulist[0].header:
    header['IMAGETYP'] = 'Light'

# check that some required FITS headers are present, otherwise raise error
#header_keys = ['DATE-OBS', 'IMAGETYP', 'CCD-TEMP', 'EXPTIME']
header_keys = ['DATE-OBS', 'CCD-TEMP', 'EXPTIME'] # temp fix
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
if not 'FILTER' in header:
    header['FILTER'] = 'Clear'

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
    filename += '-Flat-'
    filename += header['FILTER']
filename = filename.replace(' ', '_')

filename_fits = filename + '.fits'

filename_jpg_large = filename + '.jpg'
filename_jpg_thumb = filename + '_thumb.jpg'

path_jpg_large = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_large))
path_jpg_thumb = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_thumb))

os.system("/usr/bin/convert '" + os.path.abspath(json_data['fits_fname']) + "' -contrast-stretch 3%  -resize 1024x '" + path_jpg_large + "'")
os.system("/usr/bin/convert '" + os.path.abspath(json_data['fits_fname']) + "' -contrast-stretch 3%  -resize 100x '" + path_jpg_thumb + "'")

output = {
    'result': 'SUCCESS',
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
    'fwhm': 0,
    'width_px': width_px,
    'height_px': height_px,
    'ccdtemp': ccdtemp,
    'indexing_version': get_Version(),
    'filename': filename,
    'jpg_large': path_jpg_large,
    'jpg_thumb': path_jpg_thumb
}

print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))