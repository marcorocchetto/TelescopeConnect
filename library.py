import os
import json
import time
import datetime


def get_Version():
    version_file = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_version'))
    version = version_file.read().strip()
    return version

def RaiseError(error_description):
    print(json.dumps({'result': 'ERROR', 'error_description': error_description}, separators=(',',':'), indent=4))
    exit()

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
    offset_seconds = local.localize(dateobs_utc_datetime, is_dst=None).utcoffset().total_seconds()
    # get observing date in local time
    dateobs_lt_datetime = dateobs_utc_datetime + datetime.timedelta(seconds=offset_seconds)
    obsnight_datetime = dateobs_lt_datetime
    hh = dateobs_lt_datetime.strftime('%H')
    if int(hh) < 12:
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
