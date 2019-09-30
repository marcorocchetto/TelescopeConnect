import os
import sys
import datetime
import time
from astropy.io import fits
from astropy.time import Time
import numpy as np
import os
import json
import time
import datetime

def get_Version():
    version_file = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '_version'))
    version = version_file.read().strip()
    return version

def RaiseError(error_description):
    print(json.dumps({'result': 'ERROR', 'error_description': error_description}, separators=(',',':'), indent=4))
    sys.exit()

def get_FileName(header):

    dateobs_utc_str, dateobs_utc_datetime = get_DateObs(header['DATE-OBS'])
    filename = dateobs_utc_datetime.strftime('%Y-%m-%dT%H-%M-%S')
    ccdtemp = round(header['CCD-TEMP'])
    exptime = round(header['EXPTIME'])
    imagetype = get_ImageType(header['IMAGETYP'])
    if imagetype == 'LIGHT':
        if 'OBJECT' in header:
            object = header['OBJECT']
            object = object.replace('\'', '-')
            object = object.replace('/', '-')
            object = object.replace('(', '_')
            object = object.replace(')', '_')
            object = object.replace(' ', '_')
            filename += '_' + object

        filename += '_' + header['FILTER'].replace("'", "prime")
        filename += '_T'
        filename += str(ccdtemp)
        filename += '_'
        filename += str(exptime)
        filename += 's'
    elif imagetype == 'BIAS' and not 'NCOMBINE' in header:
        filename += '_Bias'
    elif imagetype == 'BIAS' and 'NCOMBINE' in header:
        filename += '_MasterBias'
    elif imagetype == 'DARK' and not 'NCOMBINE' in header:
        filename += '_Dark'
        filename += '_T'
        filename += str(ccdtemp)
        filename += '_'
        filename += str(exptime)
        filename += 's'
    elif imagetype == 'DARK' and 'NCOMBINE' in header:
        filename += '_MasterDark'
        filename += '_T'
        filename += str(ccdtemp)
        filename += '_'
        filename += str(exptime)
        filename += 's'
    elif imagetype == 'FLAT' and not 'NCOMBINE' in header:
        filename += '_Flat_'
        filename += header['FILTER'].replace("'", "prime")
        filename += '_T'
        filename += str(ccdtemp)
        filename += '_'
        filename += str(exptime)
        filename += 's'
    elif imagetype == 'FLAT' and 'NCOMBINE' in header:
        filename += '_MasterFlat_'
        filename += header['FILTER'].replace("'", "prime")
        filename += '_T'
        filename += str(ccdtemp)
        filename += '_'
        filename += str(exptime)
        filename += 's'
    filename = filename.replace(' ', '_')

    return filename


def get_ImageType(imagetyp_header):

    flat_values = ['FLAT', 'FLAT FRAME', 'FLAT FIELD']
    bias_values = ['BIAS', 'BIAS FRAME', 'ZERO']
    dark_values = ['DARK', 'DARK FRAME']
    light_values = ['LIGHT FRAME', 'LIGHT']



    imagetype = False
    for value in flat_values:
        if imagetyp_header.upper() == value:
            imagetype = 'FLAT'
    for value in dark_values:
        if imagetyp_header.upper() == value:
            imagetype = 'DARK'
    for value in bias_values:
        if imagetyp_header.upper() == value:
            imagetype = 'BIAS'
    for value in light_values:
        if imagetyp_header.upper() == value:
            imagetype = 'LIGHT'

    return imagetype


def get_ImageHeader(fname, index=0):
    hdulist = fits.open(fname, ignore_missing_end=True)
    header = hdulist[index].header
    return header

def get_ImageData(fname, index=0):
    hdulist = fits.open(fname, ignore_missing_end=True)
    imagedata = hdulist[index].data
    imagedata = np.asarray(imagedata, dtype=float)
    header = hdulist[index].header
    return imagedata, header

def get_DateObs(dateobs_utc_input):

    # Allowed formats are: '%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S',  '%Y-%m-%d', '%d/%m/%y'
    # Assume UTC
    # Returns  dateobs_utc_str, a date string of the form %Y-%m-%dT%H:%M:%S.%f

    # Allowed formats in DATE-OBS header key value
    dateformats = ['%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S', '%Y-%m-%d',
                   '%d/%m/%y']

    for dateformat in dateformats:
        try:
            dateobs_utc_structtime = time.strptime(dateobs_utc_input, dateformat)
        except ValueError:
            pass
        except OverflowError:
            pass

    if not dateobs_utc_structtime:
        RaiseError('Header key DATE-OBS is not in correct format')

    # convert struct_time to datetime
    dateobs_utc_datetime = datetime.datetime.fromtimestamp(time.mktime(dateobs_utc_structtime))
    dateobs_utc_str = dateobs_utc_datetime.strftime('%Y-%m-%dT%H:%M:%S.%f')
    return dateobs_utc_str, dateobs_utc_datetime

def get_ObsNight(dateobs_utc_datetime, local):

    # Definition: An observing night is a date defined at the midnight of the local time, starting from
    # the midday of the given date, to the midday of the day after.
    # Returns obsnight_str, a date string of the form %Y-%m-%dT

    # convert dateobs_utc_datetime to local timezone, then use pytz to get local time offset
    offset_seconds = local.localize(dateobs_utc_datetime).utcoffset().total_seconds()
    # get observing date in local time
    dateobs_lt_datetime = dateobs_utc_datetime + datetime.timedelta(seconds=offset_seconds)
    obsnight_datetime = dateobs_lt_datetime
    hh = dateobs_lt_datetime.strftime('%H')
    if int(hh) < 16:
        obsnight_datetime -= datetime.timedelta(seconds=86400)  # previous day if hour < 12
    obsnight_str = obsnight_datetime.strftime('%Y-%m-%d')
    return obsnight_str

def SexToDeg(value, coord):

    # convert OBJCTRA and OBJCDEC from sexagesimal to decimal degrees
    # assume RA: hh mm ss.sss
    # assume Dec: +/-deg arcmin arcsec.ss

    if coord.upper() == 'RA':
        ra = value.split(":")
        if len(ra) < 2: ra = value.split(" ")
        hh = float(ra[0]) * 15
        mm = (float(ra[1]) / 60) * 15
        ss = (float(ra[2]) / 3600) * 15
        return hh + mm + ss
    elif coord.upper() == 'DEC':
        dec = value.split(":")
        if len(dec) < 2: dec = value.split(" ")
        hh = abs(float(dec[0]))
        mm = float(dec[1]) / 60
        ss = float(dec[2]) / 3600
        if float(dec[0]) < 0:
            return (hh + mm + ss) * (-1)
        else:
            return (hh + mm + ss)


