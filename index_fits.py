'''

Indexing of frames

The input to the script is a json file with the following format. See input_examples

Example:

python /home/telescope/TelescopeConnect/index_fits.py --json=/home/telescope/TelescopeConnect/input_examples/index_fits.json

'''


import argparse
import pytz
from astropy.io import fits

import json

import traceback

import library.general
from library.general import *

import contextlib
import io
import sys

import sewpy
import numpy as np
import astropy

import logging
logging.basicConfig(level=logging.CRITICAL)

sew = sewpy.SEW(params=["FWHM_IMAGE"],
                config={"DETECT_MINAREA":200}
                )

#loading parameter file parser
parser = argparse.ArgumentParser()
parser.add_argument('--json',
                    dest='json_filename',
                    type=str,
                    default=False,
                    )

options = parser.parse_args()

def return_error(error):

    output = {
        'result': 'Failed',
        'message': error
    }

    print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))
    
    sys.exit()

try:

    if not options.json_filename:
        return_error('You need an inpunt JSON file')

    json_data = json.loads(open(options.json_filename).read())


    if not json_data['fits_fname']:
        return_error('File name not specified')

    if not json_data['time_zone']:
        local = pytz.timezone('Etc/GMT')
    else:
        try:
            local = pytz.timezone(json_data['time_zone'])
        except pytz.exceptions.UnknownTimeZoneError:
            local = pytz.timezone('Etc/GMT')

    # check file exists
    if not os.path.isfile(json_data['fits_fname']):
        return_error('File not found')

    # open fits file hdu
    try:
        hdulist = fits.open(json_data['fits_fname'], ignore_missing_end=True)
    except Exception as e:
        return_error(str(e) + "Traceback: " + traceback.format_exc())


    # get header & data
    header = hdulist[0].header

    try:
        data = hdulist[0].data
    except Exception as e:
        return_error(str(e) + "Traceback: " + traceback.format_exc())

    # temporary fix for PinPoint simulated images
    if not 'IMAGETYP' in hdulist[0].header:
        header['IMAGETYP'] = 'Light'

    # check that some required FITS headers are present, otherwise raise error
    #header_keys = ['DATE-OBS', 'IMAGETYP', 'CCD-TEMP', 'EXPTIME']
    header_keys = ['DATE-OBS', 'CCD-TEMP', 'EXPTIME'] # temp fix
    for key in header_keys:
        if not key in hdulist[0].header:
            return_error('Header key `%s` missing' % key)


    # determine Image type from header IMAGETYP
    imagetype = get_ImageType(header['IMAGETYP'])
    if not imagetype:
        return_error('Cannot determine ImageType')

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
    if not 'CCD-TEMP' in header:
        ccdtemp = 0
    else:
        ccdtemp = header['CCD-TEMP']

    # image size
    width_px = len(data[0,:])
    height_px = len(data[:,0])

    # set filename
    filename = dateobs_utc_datetime.strftime('%Y-%m-%dT%H-%M-%S')
    if imagetype == 'LIGHT':
        if 'OBJECT' in header:
            filename += '-' + header['OBJECT']
        filename += '-' + header['FILTER'].replace("'", "prime")
    if imagetype == 'BIAS':
        filename += '-Bias'
    if imagetype == 'DARK':
        filename += '-Dark'
    if imagetype == 'FLAT':
        filename += '-Flat-'
        filename += header['FILTER'].replace("'", "prime")
    filename = filename.replace(' ', '_')

    filename_fits = filename + '.fits'

    filename_jpg_large = filename + '.jpg'
    filename_jpg_thumb = filename + '_thumb.jpg'

    path_jpg_large = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_large))
    path_jpg_thumb = os.path.abspath(os.path.join(json_data['output_folder'], filename_jpg_thumb))

    os.system("/usr/bin/convert '" + os.path.abspath(json_data['fits_fname']) + "' -contrast-stretch 3%  -resize 1024x '" + path_jpg_large + "'")
    os.system("/usr/bin/convert '" + os.path.abspath(json_data['fits_fname']) + "' -contrast-stretch 3%  -resize 100x '" + path_jpg_thumb + "'")

    # get seeing
    if not 'pixel_scale' in json_data:
        seeing = 0
    else:

        PixelScale = json_data['pixel_scale']

        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                out = sew(json_data['fits_fname'])
                seeing = round(np.median(out["table"][0][:]) * PixelScale, 1)

    # clean fits headers
    if 'WXSENSOR' in header:
        del header['WXSENSOR']
    if 'SWSERIAL' in header:
        del header['SWSERIAL']
    if 'HISTORY' in header:
        del header['HISTORY']
    if 'COMMENT' in header:
        del header['COMMENT']
    if 'PRESSURE' in header:
        del header['PRESSURE']
    if 'SBSTDVER' in header:
        del header['SBSTDVER']
    if 'SWOWNER' in header:
        del header['SWOWNER']
    if 'PLTSOLVD' in header:
        del header['PLTSOLVD']

    if 'telescope_model' in json_data:
        header['TELESCOP'] = json_data['telescope_model']

    if 'imager_name' in json_data:
        header['INSTRUME'] = json_data['imager_name']

    header['OBSERVER'] = 'TelescopeLive'

    if 'telescope_guid' in json_data:
        header['TELID'] = json_data['telescope_guid']

    if 'imager_guid' in json_data:
        header['IMAGERID'] = json_data['imager_guid']

    if 'frame_guid' in json_data:
        header['FRAMEID'] = json_data['frame_guid']

    # save modified fits
    hdulist[0].header = header
    hdulist.writeto(json_data['fits_fname'], overwrite=True)

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
        'width_px': width_px,
        'height_px': height_px,
        'ccdtemp': ccdtemp,
        'indexing_version': get_Version(),
        'filename': filename,
        'jpg_large': path_jpg_large,
        'jpg_thumb': path_jpg_thumb,
        'qa_passed': False,
        'qa_pointing': False,
        'qa_seeing': False,
        'poor_fwhm': False,
        'background_gradient': False,
        'no_stars': False,
        'ellipticity': 0,
        'strikes': False,
        'star_number': 0,
        'seeing': seeing,
        'internal_reflection': False,
        'flat_field_residuals': False,
        'rbi_residuals': False,
        'saturated ': False
    }

    print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))

except Exception as e:

    return_error(str(e) + "Traceback: " + traceback.format_exc())
