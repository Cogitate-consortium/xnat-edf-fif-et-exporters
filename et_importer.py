#!/usr/bin/env python
# author: Praveen Sripad <praveen.sripad@ae.mpg.de>
# author: Mohana Ramaratnam <mohanaramaratnam@flywheel.io>


"""Extract headers from eye tracker data files for import into XNAT.

--------
.. code-block:: console
    $ et_importer.py -f SF105_ET_V1_DurR1.csv -o headers.json


Eye tracking data uploaded to Cogitate XNAT may be of two types - csv or asc.
The asc file is the result of converting Eyelink Data Format (EDF)
produced by Eyelink eye tracker. The csv file is produced by the Tobii eye tracker.

1. Reads asc / csv file and extracts header information from each.
2. Provides an option to limit output to five important fields that will be modeled on XNAT.
3. ET data header extraction is highly device specific, so this is script tries
   to be as explicit and descriptive as possible.
"""

import sys
import os.path as op
import json

def generateETProppertiesFileForXnat(et_headers, out_fname):
        print('Generating Properties File for XNAT: ' + out_fname)

        out_file = open(r'%s'% out_fname,"w")
        out_file.write('gnmd_isAnon=%s\n'% et_headers['isAnon'])
        out_file.write('bids_Manufacturer=%s\n'% et_headers['Manufacturer'])
        out_file.write('bids_SamplingFrequency=%s\n'% et_headers['SamplingFrequency'])
        out_file.write('gnmd_DistanceToScreen=%s\n'% et_headers['DistanceToScreen'])
        out_file.write('gnmd_ScreenSize=%s\n'% et_headers['ScreenSize'])
        out_file.write('gnmd_RecordedEye=%s\n'% et_headers['RecordedEye'])
        out_file.close()


def extract_asc_file_headers(input_fname, out_fname, xnat_fields_only=True,
                             generate_properties_file_for_xnat=True):
    """Extract file headers from the given ascii (.asc) file for EyeLink eye tracker.

    The files are from proprietary Eyelink eye tracker devices converted from
    Eyelink Data Format to ascii, saved with the extension .asc.
    
    Parameters
    ----------
    input_fname : str
        Input file name of type .asc. Should be the output of the edf2asc tool provided by Eyelink. 
    out_fname : str
        Output file name
        relative paths are saved relative to parent dir of fif_fname
        None will save to parent dir of fif_fname with default prefix
    xnat_fields_only: bool
        Set to True to return only a small set of fields to be modeled on XNAT.
    generate_properties_file_for_xnat: bool
        Generates the properties file used by XNAT
    """

    # tons of explicitly declared headers to be extracted
    header_fields_to_extract = ['TYPE', 'VERSION', 'SOURCE', 'CAMERA',
                                'SERIAL NUMBER', 'CAMERA_CONFIG']
    msg_fields_to_extract = ['Calibration_area', 'Screen_size_mm',
                             'Screen_distance_mm', 'Camera_position_mm']
    config_fields_to_extract = ['RECCFG', 'ELCLCFG', 'GAZE_COORDS', 'THRESHOLDS',
                                'ELCL_WINDOW_SIZES', 'ELCL_PROC', 'PUPIL', 'PUPIL_DATA_TYPE',
                                'ELCL_PCR_PARAM']
    et_xnat_fields = ['Manufacturer', 'SamplingFrequency', 'DistanceToScreen',
                      'ScreenSize', 'RecordedEye']

    # read asc file and extract available header information
    with open(input_fname, 'r') as f:
        asc = [line.rstrip('\n') for line in f]

    et_config_dict = dict()

    # loop through the interested fields, loop through full data
    # extract field and make sure it does not change throughout the data
    et_config_dict['Manufacturer'] = 'EYELINK'

    for header_field in header_fields_to_extract:
        if not any([1 if myl.find(header_field + ':') != -1 else 0 for myl in asc]):
            et_config_dict[header_field] = 'N/A'

    # extract header fields
    for header_field in header_fields_to_extract:
        for myl in asc:
            if myl.find(header_field + ':') != -1:
                if et_config_dict.get(header_field):
                    assert et_config_dict[header_field] == myl.split(': ')[-1], 'Different values found for the same header - %s.' % header_field
                else:
                    et_config_dict[header_field] = myl.split(': ')[-1]

    for msg_field in msg_fields_to_extract:
        if not any([1 if myl.find(msg_field + ':') != -1 else 0 for myl in asc]):
            et_config_dict[msg_field] = 'N/A'

    # extract MSG fields
    for msg_field in msg_fields_to_extract:
        for myl in asc:
            if myl.find(msg_field) != -1:
                if et_config_dict.get(msg_field):
                    assert et_config_dict[msg_field] == myl.split(msg_field + ': ')[-1], 'Different values found for the same header - %s.' % msg_field
                else:
                    et_config_dict[msg_field] = myl.split(msg_field + ': ')[-1]
    
    for config_field in config_fields_to_extract:
        if not any([1 if myl.find(config_field) != -1 else 0 for myl in asc]):
            et_config_dict[config_field] = 'N/A'

    # extract config fields
    for config_field in config_fields_to_extract:
        for myl in asc:
            if myl.find(config_field) != -1:
                if config_field == 'PUPIL' and 'PUPIL_DATA_TYPE' in myl:
                    # special case to make sure it is not the PUPIL_DATA_TYPE and PUPIL only
                    # this is the problem with string matching
                    continue
                if et_config_dict.get(config_field):
                    assert et_config_dict[config_field] == myl.split(config_field)[-1].strip(), 'Different values found for the same header - %s.' % config_field
                else:
                    et_config_dict[config_field] = myl.split(config_field)[-1].strip()

    # post processing include additional information on the meaning of the fields

    # RECCFG <tracking mode> <sampling rate> <file sample filter> <link sample filter> <eye(s) recorded>
    # This specifies the tracking mode used (pupil-only vs. pupil-CR mode), 
    # sampling rate (250, 500, 1000, or 2000 hz), file sample filter 
    # (0 – filter off; 1 – standard filter; 2 – extra filter), 
    # link/analog filer (0 – filter off; 1 – standard filter; 2 – extra filter), 
    # and the eyes (L, R, or LR) recorded in the trial.
    reccfg_fields = ['tracking_mode', 'sampling_rate', 'file_sample_filter',
                     'link_sample_filter', 'eyes_recorded']
    for myfield, myval in zip(reccfg_fields, et_config_dict['RECCFG'].split(' ')):
        et_config_dict[myfield] = myval

    # ELCLCFG <mount configuration>
    # mount_configuration, mount_configuration_description
    mount_configuration_dict = {
        "MTABLER": "Desktop, Stabilized Head, Monocular",
        "BTABLER": "Desktop, Stabilized Head, Binocular/Monocular",
        "RTABLER": "Desktop (Remote mode), Target Sticker, Monocular",
        "RBTABLER": "Desktop (Remote mode), Target Sticker, Binocular/Monocular",
        "AMTABLER": "Arm Mount, Stabilized Head, Monocular",
        "ARTABLER": "Arm Mount (Remote mode), Target Sticker, Monocular",
        "BTOWER": "Binocular Tower Mount, Stabilized Head, Binocular/Monocular",
        "TOWER": "Tower Mount, Stabilized Head, Monocular",
        "MPRIM": "Primate Mount, Stabilized Head, Monocular",
        "BPRIM": "Primate Mount, Stabilized Head, Binocular/Monocular",
        "MLRR": "Long-Range Mount, Stabilized Head, Monocular, Camera Level",
        "BLRR": "Long-Range Mount, Stabilized Head, Binocular/Monocular, Camera Angled"
    }
    et_config_dict['mount_configuration'] = et_config_dict['ELCLCFG']
    et_config_dict['mount_configuration_description'] = mount_configuration_dict[et_config_dict['ELCLCFG']]

    # GAZE_COORDS <left> <top> <right> <bottom>
    # This reports the pixel resolution of the tracker recording. 
    # Left, top, right, and bottom refer to the x-y coordinates of the top-left and bottom-right corners of display.
    # et_config_dict['GAZE_COORDS_description'] = "GAZE_COORDS <left> <top> <right> <bottom>. 
    # This reports the pixel resolution of the tracker recording. Left, top, right, 
    # and bottom refer to the x-y coordinates of the top-left and bottom-right corners of display."

    # THRESHOLDS <eye> <pupil> <CR>
    # This reports the pupil and CR thresholds of the tracked eye(s).
    # et_config_dict['THRESHOLDS_description'] = "THRESHOLDS <eye> <pupil> <CR>. 
    # This reports the pupil and CR thresholds of the tracked eye(s)."

    # ELCL_PROC <pupil tracking algorithm>
    # This reports the pupil fitting processing type (i.e., ELLIPSE or CENTROID)
    et_config_dict['pupil_tracking_algorithm'] = et_config_dict['ELCL_PROC']

    # PUPIL <data type>
    # This specifies the type of pupil size data (AREA or DIAMETER) recorded.
    # et_config_dict['pupil_data_type_description'] = "PUPIL <data type>.
    # This specifies the type of pupil size data (AREA or DIAMETER) recorded in the trial."
    for myl in asc:
        if myl.find('PUPIL\t') != -1:
            et_config_dict['pupil_size_data'] = myl.split('\t')[-1]

    # add some bids fields
    et_config_dict['ManufacturersModelName'] = et_config_dict['VERSION']
    et_config_dict['Eye tracker camera'] = et_config_dict['CAMERA']
    et_config_dict['DeviceSerialNumber'] = et_config_dict['SERIAL NUMBER']
    et_config_dict['SamplingFrequency'] = et_config_dict['sampling_rate']

    if et_config_dict.get('Screen_size_mm') == 'N/A':
        et_config_dict['ScreenSize'] = 'N/A'
    else:
        # get the screen size in cms, some weirdfu here !
        et_config_dict['ScreenSize'] = \
            str(float(et_config_dict['Screen_size_mm'].split()[0].lstrip('-')) / 10) + ' ' + \
            str(float(et_config_dict['Screen_size_mm'].split()[1].lstrip('-')) / 10) + ' cms'

    if et_config_dict.get('Screen_distance_mm') == 'N/A':
        et_config_dict['DistanceToScreen'] = 'N/A'
    else:
        # get the distance to screen in cms
        # (pick one value only and convert from mm)
        et_config_dict['DistanceToScreen'] = \
            str(float(et_config_dict['Screen_distance_mm'].split(' ')[0].strip('-')) / 10) + ' cms'

    # convert eyes_recorded by Eyelink into BIDS friendly RecordedEye
    if et_config_dict['eyes_recorded'] == "LR":
        et_config_dict['RecordedEye'] = "Both"
    elif et_config_dict['eyes_recorded'] == "L":
        et_config_dict['RecordedEye'] = "Left"
    elif et_config_dict['eyes_recorded'] == "R":
        et_config_dict['RecordedEye'] = "Right"
    else:
        raise ValueError('eyes_recorded has to be L, R or LR.')

    if not all(etf in et_config_dict for etf in et_xnat_fields):
        raise ValueError('Mandatory fields are missing.')

    # check if 'CONVERTED FROM and DATE' strings are in the data
    if all(["CONVERTED FROM" not in line and "DATE" not in line for line in asc]):
        et_config_dict['isAnon'] = True

    if xnat_fields_only:
        et_output = {key: et_config_dict[key] for key in et_config_dict.keys() & et_xnat_fields}
    else:
        et_output = et_config_dict

    # out_fname = op.basename(input_file).split('.')[0] + '.json'
    with open(out_fname, "w") as f:
        json.dump(et_output, f, indent=4)

    if generate_properties_file_for_xnat:
        prop_file_name = out_fname + '.prop'
        generateETProppertiesFileForXnat(et_config_dict, prop_file_name)

    return et_output


def extract_csv_file_headers(input_fname, out_fname, xnat_fields_only=True,
                             generate_properties_file_for_xnat=True):
    """Extract file headers from the given asc file for Tobii eye tracker.

    Parameters
    ----------
    input_fname : str
        Input file name of type .csv. This is a csv file produced by other scripts inhouse
        specifically for the Tobii eye tracker.
    out_fname : str
        Output file name
        relative paths are saved relative to parent dir of fif_fname
        None will save to parent dir of fif_fname with default prefix
    xnat_fields_only: bool
        Set to True to return only a small set of fields to be modeled on XNAT.
    generate_properties_file_for_xnat: bool
        Generates the properties file used by XNAT
    """

    import pandas as pd
    import numpy as np

    def is_unique(s):
        """Test first row against all other rows of df."""
        a = s.to_numpy()  # s.values (pandas<0.24)
        return (a[0] == a).all()

    et_xnat_fields = ['Manufacturer', 'SamplingFrequency', 'DistanceToScreen',
                      'ScreenSize', 'RecordedEye']

    et_config_dict = dict()

    # read the csv file
    df = pd.read_csv(input_fname)

    # get list of interested headers
    tobii_headers = df.columns[1:6].to_list()

    et_config_dict['Manufacturer'] = 'TOBII'

    # extract the headers
    for header_field in tobii_headers:
        assert is_unique(df[header_field]), 'Different values found for the same header - %s.' % header_field
        # get the first value of the field
        et_config_dict[header_field] = str(df[header_field][0])

    # get the sampling rate from the device timestamp
    # get mean difference between consecutive timestamps
    # convert them to seconds from microseconds and get the sampling rate
    et_config_dict['SamplingFrequency'] = str(np.round(1 / (np.round(df['device_time_stamp'].diff(periods=1).mean(), 2) / 1e6), 2))

    # Tobii is always recorded with both eyes / binocular
    # TODO hardcoded value because there is no information in this csv
    et_config_dict['RecordedEye'] = "Both"

    et_config_dict['ScreenSize'] = et_config_dict['ScreenHeightCm'] + ' ' + \
        et_config_dict['ScreenWidthCm'] + ' cms'

    if not all(etf in et_config_dict for etf in et_xnat_fields):
        raise ValueError('Mandatory fields are missing.')

    # the csv file is a product of an inhouse script that strips everything
    # it should be anonymized by default
    # TODO hardcoded value because there is no information in this csv
    et_config_dict['isAnon'] = True

    if xnat_fields_only:
        et_output = {key: et_config_dict[key] for key in et_config_dict.keys() & et_xnat_fields}
    else:
        et_output = et_config_dict

    # out_fname = op.basename(input_file).split('.')[0] + '.json'
    with open(out_fname, "w") as f:
        json.dump(et_output, f, indent=4)

    if generate_properties_file_for_xnat:
        prop_file_name = out_fname + '.prop'
        generateETProppertiesFileForXnat(et_config_dict, prop_file_name)

    return et_output


def run():
    """Run ET files to extract headers."""
    from optparse import OptionParser
    parser = OptionParser()

    parser.add_option("-f", "--file", type="string", dest="file",
                      help="Name of raw data file.", metavar="FILE",
                      default=None)
    parser.add_option("-o", "--output", type="string", dest="output",
                      help="Name of extracted headers output file.",
                      metavar="OUTFILE", default=None)
    parser.add_option("-x", "--xProperties",  action="store_true", dest="xProperties",
                      help="Additionaly generate properties file required by XNAT.",
                      metavar="GENERATEPROPERTIESFILEFORXNAT", default=False)

    options, args = parser.parse_args()

    if options.file is None:
        print('Input file not provided.')
        parser.print_help()
        sys.exit(1)

    fname = options.file

    # check/set output header name
    if options.output is None:
        out_fname = op.basename(fname) + '.json'
    elif not options.output.endswith('.json'):
        out_fname = options.output + '.json'
    else:
        out_fname = options.output

    et_template_json_string = '{"InstitutionName": "", "InstitutionAddress": "", "Manufacturer": "", "ManufacturersModelName": "", "SoftwareVersions": "", "TaskDescription": "", "Instructions": "", "CogAtlasID": "", "CogPOID": "", "DeviceSerialNumber": "", "Firmware #": "", "SamplingFrequency": null, "SampleCoordinateUnit": "", "SampleCoordinateSystem": "", "EnvironmentCoordinates": [], "EventIdentifier": "", "RawSamples": "", "IncludedEyeMovementEvents": [], "DetectionAlgorithm": "", "Eye tracker camera": "", "Eye tracker lens": "", "StartMessage": "", "EndMessage": "", "KeyPressMessage": "", "CalibrationType": "", "CalibrationPosition": null, "CalibrationUnit": "", "MaximalCalibrationError": null, "AverageCalibrationError": null, "CalibrationList": [], "RecordedEye": "", "EyeCameraSettings": [], "FeatureDetectionSettings": [], "GazeMappingSettings": [], "DetectionAlgorithmSettings": [], "RawDataFilters": "", "ScreenSize": "", "ScreenResolution": "", "ScreenRefreshRate": "", "AOIDefinition": [], "PupilPositionType": "", "PupilFitMethod": ""}'

    generate_properties_file_for_xnat = False
    if options.xProperties:
        generate_properties_file_for_xnat = True

    xnat_fields_only = True

    # ET data can be .asc or .csv files
    if op.splitext(fname)[-1] == '.asc':
        et_config_dict = extract_asc_file_headers(fname, out_fname, xnat_fields_only,
                                                  generate_properties_file_for_xnat)
    elif op.splitext(fname)[-1] == '.csv':
        et_config_dict = extract_csv_file_headers(fname, out_fname, xnat_fields_only,
                                                  generate_properties_file_for_xnat)
    else:
        raise IOError('Input file should be of type asc or csv.')

    print('Eye Tracking ASC/CSV file header extraction complete.')


is_main = (__name__ == '__main__')
if is_main:
    run()
