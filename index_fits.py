'''

Indexing of frames

The input to the script is a json file with the following format. See input_examples

Example:

python /home/telescope/TelescopeConnect/index_fits.py --json=/home/telescope/TelescopeConnect/input_examples/index_fits.json

'''

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append('/home/telescope/anaconda3/lib/python3.6/site-packages')

import argparse
import pytz
from astropy import log
import traceback
from library.general import *
import io
import logging
logging.basicConfig(level=logging.CRITICAL)

import os
import shutil
import sewpy
import subprocess
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy import wcs
from astropy import units
from astroquery.xmatch import XMatch
import math

#loading parameter file parser
parser = argparse.ArgumentParser()
parser.add_argument('--json',
                    dest='json_filename',
                    type=str,
                    default=False,
                    )

options = parser.parse_args()

solvefield = '/usr/local/astrometry/bin/solve-field'

log.setLevel('ERROR')

save_stdout = sys.stdout
# sys.stdout = io.StringIO()


def return_error(error):

    output = {
        'result': 'Failed',
        'message': error
    }

    print(error)

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
        hdu = fits.open(json_data['fits_fname'], ignore_missing_end=True)
    except Exception as e:
        return_error(str(e) + "Traceback: " + traceback.format_exc())

    # init var
    qa_pointing = False
    qa_seeing = False
    qa_passed = False
    seeing = 0

    # get header & data
    header = hdu[0].header

    try:
        data = hdu[0].data
    except Exception as e:
        return_error(str(e) + "Traceback: " + traceback.format_exc())

    # temporary fix for PinPoint simulated images
    if not 'IMAGETYP' in hdu[0].header:
        header['IMAGETYP'] = 'Light'

    # check that some required FITS headers are present, otherwise raise error
    #header_keys = ['DATE-OBS', 'IMAGETYP', 'CCD-TEMP', 'EXPTIME']
    header_keys = ['DATE-OBS', 'CCD-TEMP', 'EXPTIME', 'XBINNING', 'YBINNING']
    for key in header_keys:
        if not key in hdu[0].header:
            return_error('Header key `%s` missing' % key)

    # determine Image type from header IMAGETYP
    imagetype = get_ImageType(header['IMAGETYP'])
    if not imagetype:
        return_error('Cannot determine ImageType')

    # Date and time of observation: keyword DATE-OBS
    dateobs_utc_str, dateobs_utc_datetime = get_DateObs(header['DATE-OBS'])

    if imagetype != 'LIGHT':

        # it's a calibration frame
        qa_passed = True
        qa_pointing = True
        seeing = 0

    else:

        if 'max_allowed_seeing'in json_data:
            max_allowed_seeing = json_data['max_allowed_seeing']
        else:
            max_allowed_seeing = 99

        # get seeing
        seeing = 0
        if not 'pixel_scale' in json_data:
            seeing = 0
        else:

            PixelScale = json_data['pixel_scale']

            workenv = json_data['output_folder']

            margin = 200
            scalelow = PixelScale * 0.9
            scalehigh = PixelScale * 1.1
            hdu = fits.open(os.path.abspath(json_data['fits_fname']))
            xaxis = hdu[0].header['NAXIS1']
            yaxis = hdu[0].header['NAXIS2']
            binning = hdu[0].header['XBINNING']
            ra = hdu[0].header['OBJCTRA'].replace(' ', ':')
            dec = hdu[0].header['OBJCTDEC'].replace(' ', ':')
            radius = 5*round(np.sqrt(2) * xaxis * PixelScale / 2 / 60)  # radius of search in arcmin (5 field diagonals)

            args = '--verbose --use-sextractor --tweak-order 2 ' \
                   '--ra %s --dec %s --radius 2 ' \
                   '--dir %s --no-plots --overwrite ' \
                   ' --scale-units arcsecperpix  --scale-low %.2f --scale-high %.2f' \
                   ' %s' % (ra, dec, workenv, scalelow, scalehigh, os.path.abspath(json_data['fits_fname']))

            with open(os.devnull, 'wb') as devnull:
                try:
                    output = subprocess.check_output(solvefield + ' ' + args,
                                                     stderr=subprocess.STDOUT, shell=True)
                except subprocess.CalledProcessError as e:
                    return_error(e.output)

            # check if test.solved exists
            if not os.path.isfile(os.path.join(workenv, os.path.splitext(os.path.basename(os.path.abspath(json_data['fits_fname'])))[0] + '.solved')):

                seeing = 0
                platesolved = False

            else:

                platesolved = True

                filename_new = os.path.join(workenv, os.path.splitext(os.path.basename(os.path.abspath(json_data['fits_fname'])))[0] + '.new')
                filename_final = os.path.abspath(json_data['fits_fname'])

                shutil.move(filename_new, filename_final)

                hdu = fits.open(filename_final)
                image = hdu[0].data
                header = hdu[0].header

                try:

                    sew = sewpy.SEW(params=["FWHM_IMAGE", "X_WORLD", "Y_WORLD", "X_IMAGE", "Y_IMAGE", "BACKGROUND", "FLUX_ISOCOR"],
                                    config={"DETECT_MINAREA": 50, "DETECTION_TRESH": 2})
                    sewout = sew(filename_final)

                    sort = sewout['table'].argsort(['FLUX_ISOCOR'])[::-1]
                    sorted_table = sewout['table'][:][sort]
                    fieldstars_world = []
                    fieldstars_image = []
                    fieldstars_seeing = []
                    fieldstars_background = []

                    # remove stars with bad background
                    bgrm_table = sorted_table[sorted_table['BACKGROUND'] < np.median(image)]

                    # remove edge stars
                    for fieldstar in bgrm_table[:150]: # take top 150
                        if fieldstar["X_IMAGE"] > margin and fieldstar["X_IMAGE"] < (xaxis - margin) and \
                                fieldstar["Y_IMAGE"] > margin and fieldstar["Y_IMAGE"] < (yaxis - margin):
                            fieldstars_world.append((fieldstar["X_WORLD"], fieldstar["Y_WORLD"]))
                            fieldstars_image.append((fieldstar["X_IMAGE"], fieldstar["Y_IMAGE"]))
                            fieldstars_seeing.append(fieldstar["FWHM_IMAGE"])
                            fieldstars_background.append(fieldstar["BACKGROUND"])

                    ra = hdu[0].header['CRVAL1']
                    dec = hdu[0].header['CRVAL2']
                    c = SkyCoord(ra, dec, unit="deg")

                    fn_fieldstars = os.path.join(workenv, os.path.splitext(os.path.basename(filename_new))[0] + '.fieldstars')
                    f = open(fn_fieldstars, "w+")
                    f.write('ra, dec,idx \n')
                    for idx, result in enumerate(fieldstars_world):
                        f.write(str(result[0]) + ', ' + str(result[1]) + ',' + str(idx) + '\n')
                    f.close()

                    table = XMatch.query(cat1=open(fn_fieldstars),
                                         cat2='vizier:II/246/out',
                                         max_distance=5 * u.arcsec,
                                         colRA1='ra',
                                         colDec1='dec')

                    seeing_values = []
                    for idx in table['idx']:
                        seeing_values.append(fieldstars_seeing[idx])
                    seeing = np.median(seeing_values) * PixelScale * int(binning)


                except Exception as e:

                    seeing = 0

        if not seeing or math.isnan(seeing):
            seeing = 0

        if seeing > 0:
            if seeing > max_allowed_seeing:
                qa_seeing = False
                poor_fwhm = True
            else:
                qa_seeing = True
        else:
            qa_seeing = False



        # check pointing
        if platesolved:
            imgwcs = wcs.WCS(hdu[0].header)
            coord1 = wcs.utils.pixel_to_skycoord(xaxis / 2, yaxis / 2, imgwcs)

            coord2 = SkyCoord(hdu[0].header['OBJCTRA'], hdu[0].header['OBJCTDEC'],
                              unit=(units.hourangle, units.degree)).transform_to('fk5')
            separation = coord1.separation(coord2).arcsecond

            if separation > 0.1 * PixelScale * xaxis:
                qa_pointing = False
            else:
                qa_pointing = True

        else:
            qa_pointing = False

        if not qa_seeing or not qa_pointing:
            qa_passed = False
        else:
            qa_passed = True

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

    # # filter
    # if not 'FILTER' in header:
    #     header['FILTER'] = 'Clear'

    if 'FILTER' in header and imagetype != 'BIAS':
        filter = header['FILTER']
    else:
        filter = ''

    # fix filter names (chilescope)
    if filter == 'Lum':
        filter = 'Luminance'
    elif filter == 'H-alpha':
        filter = 'Halpha'
    elif filter == 'Oiii':
        filter = 'OIII'
    elif filter == 'Sii':
        filter = 'SII'
    elif filter == 'SLOAN':
        filter = 'Sloan r'
    elif filter == 'r-Sloan':
        filter = 'Sloan r'

    header['FILTER'] = filter

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
            object = header['OBJECT']
            object = object.replace('\'', '-')
            object = object.replace('/', '-')
            object = object.replace('(', '_')
            object = object = object.replace(')', '_')
            object = object.replace(' ', '_')
            filename += '-' + object
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

    # get binning
    binning = header['XBINNING'] # assume YBINNING is the same

    # clean fits headers
    if 'WXSENSOR' in header:
        del header['WXSENSOR']
    if 'SWSERIAL' in header:
        del header['SWSERIAL']
    if 'SWCREATE' in header:
        del header['SWCREATE']
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

    observationBlockGuid = None
    if 'PLAN' in header:
        if len(header['PLAN'].split('_')) > 1:
            observationBlockGuid = header['PLAN'].split('_')[1]

    # save modified fits
    hdu[0].header = header


    hdu.writeto(json_data['fits_fname'], overwrite=True)

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
        'qa_passed': qa_passed,
        'qa_pointing': qa_pointing,
        'qa_seeing': qa_seeing,
        'background_gradient': False,
        'no_stars': False,
        'ellipticity': 0,
        'strikes': False,
        'star_number': 0,
        'seeing': round(seeing, 1),
        'internal_reflection': False,
        'flat_field_residuals': False,
        'rbi_residuals': False,
        'saturated ': False,
        'binning': binning,
    }
    if observationBlockGuid:
        output['observationBlockGuid'] = observationBlockGuid


    sys.stdout = save_stdout

    print(json.dumps(output, separators=(',',':'), sort_keys=True, indent=4))

except Exception as e:

    return_error(str(e) + "Traceback: " + traceback.format_exc())
