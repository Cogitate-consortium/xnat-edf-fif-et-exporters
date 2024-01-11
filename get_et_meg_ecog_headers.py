#!/usr/bin/env python
# author: Mohana Ramaratnam <mohanaramaratnam@flywheel.io>
# author: Praveen Sripad <praveen.sripad@ae.mpg.de>

"""Extract FIF/EDF file header into a properties file.

--------
.. code-block:: console
    $ get_et_meg_ecog_headers.py -f sample_audvis_raw.fif [-x] [-v]
"""

import os.path as op
import sys
import json
import mne


def generateFIFProppertiesFileForXnat(meg_headers, out_fname):
    """Returns MEG / FIF headers for xnat."""
    print('Generating Properties File for XNAT ' + out_fname)

    isAnon = False
    if meg_headers['experimenter'] == 'XXX':
        isAnon = True

    out_file = open(r'%s'% out_fname,"w")
    out_file.write('gnmd_isAnon=%s\n'% isAnon)
    out_file.write('bids_sfrequency=%s\n'% meg_headers['SamplingFrequency'])
    out_file.write('bids_powerLineFrequency=%s\n'% meg_headers['PowerLineFrequency'])
    out_file.write('bids_dewarPosition=%s\n'% meg_headers['DewarPosition'])
    out_file.write('bids_softwareFilters_default_highpass=%s\n'% meg_headers['SoftwareFilters']['default']['highpass'])
    out_file.write('bids_softwareFilters_default_lowpass=%s\n'% meg_headers['SoftwareFilters']['default']['lowpass'])
    out_file.write('bids_digitizedLandmarks=%s\n'% meg_headers['DigitizedLandmarks'])
    out_file.write('bids_digitizedHeadPoints=%s\n'% meg_headers['DigitizedHeadPoints'])
    out_file.write('bids_megChannelCount=%s\n'% meg_headers['MEGChannelCount'])
    out_file.write('bids_recordingDuration=%s\n'% meg_headers['RecordingDuration'])
    out_file.write('bids_numberOfChannels=%s\n'% meg_headers['NumberOfChannels'])
    out_file.write('bids_megCoordinateSystem=%s\n'% meg_headers['MEGCoordinateSystem'])
    out_file.close()


def generateEDFProppertiesFileForXnat(ecog_headers, raw, out_fname):
    """Returns EDF headers for xnat."""
    print('Generating Properties File for XNAT' + out_fname)

    isAnon = False
    if ecog_headers['patient_id'] == 'X X X X_X':
        isAnon = True

        print(raw.info['subject_info'] )
        out_file = open(r'%s' % out_fname, "w")
        out_file.write('gnmd_isAnon=%s\n' % isAnon)
        out_file.write('bids_frequency=%s\n'% ecog_headers['SamplingFrequency'])
        out_file.write('bids_numberOfSignals=%s\n'% ecog_headers['NumberOfChannels'])
        out_file.write('bids_hardwareFilters_HighpassFilter_CutoffFrequency=%s\n'% ecog_headers['HardwareFilters']['HighpassFilter']['CutoffFrequency'])
        out_file.write('bids_hardwareFilters_LowpassFilter_CutoffFrequency=%s\n'% ecog_headers['HardwareFilters']['LowpassFilter']['CutoffFrequency'])
        out_file.write('bids_recordingDuration=%s\n'% ecog_headers['RecordingDuration'])
        out_file.close()


def generateETProppertiesFileForXnat(et_headers, out_fname):
    """Returns eye tracker headers for xnat."""
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
    """Extract file headers from the given csv file for Tobii eye tracker.

    The files are typically generated from the output of the Tobii eye tracker
    and hand engineered into csv files.

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


def fif_header_to_properties(fif_fname, json_template, out_fname, 
                             generate_properties_file_for_xnat):
    """Call *fif_header_to_properties* on fif file and save json output.

    Parameters
    ----------
    fif_fname : str
        Raw fif File
    json_template: str
        BIDS template used to populate various fields.
    out_fname : str
        Output file name
        relative paths are saved relative to parent dir of fif_fname
        None will save to parent dir of fif_fname with default prefix
    generate_properties_file_for_xnat : flag
        Generates the properties file used by XNAT
    """
    # load raw object
    raw = mne.io.read_raw_fif(fif_fname, verbose='ERROR')

    meg_headers = json.loads(json_template)

    if raw.info['device_info']:
        meg_headers['InstitutionName'] = raw.info['device_info']['site']
        meg_headers['Manufacturer'] = raw.info['device_info']['type']
        meg_headers['ManufacturersModelName'] = raw.info['device_info']['model']
        meg_headers['DeviceSerialNumber'] = raw.info['device_info']['serial']

    meg_headers['NumberOfChannels'] = raw.info['nchan']
    meg_headers['SamplingFrequency'] = raw.info['sfreq']
    meg_headers['PowerLineFrequency'] = raw.info['line_freq']
    meg_headers['DewarPosition'] = raw.info['gantry_angle']
    meg_headers['SoftwareFilters'] = {"default": {"highpass": raw.info['highpass'],
                                                   "lowpass": raw.info["lowpass"]}}

    if raw.info['dig'] is not None and len(raw.info['dig']) != 0:
        meg_headers['DigitizedLandmarks'] = True
        meg_headers['DigitizedHeadPoints'] = True

    # channel counts
    meg_headers['MEGChannelCount'] = len(mne.pick_types(raw.info, meg=True))
    meg_headers['MEGREFChannelCount'] = len(mne.pick_types(raw.info, ref_meg=True))
    meg_headers['EEGChannelCount'] = len(mne.pick_types(raw.info, eeg=True))
    meg_headers['ECOGChannelCount'] = len(mne.pick_types(raw.info, ecog=True))
    meg_headers['SEEGChannelCount'] = len(mne.pick_types(raw.info, seeg=True))
    meg_headers['EOGChannelCount'] = len(mne.pick_types(raw.info, eog=True))
    meg_headers['ECGChannelCount'] = len(mne.pick_types(raw.info, ecg=True))
    meg_headers['EMGChannelCount'] = len(mne.pick_types(raw.info, emg=True))
    meg_headers['MiscChannelCount'] = len(mne.pick_types(raw.info, misc=True))
    meg_headers['TriggerChannelCount'] = len(mne.pick_types(raw.info, stim=True))

    # recording duration = total samples / sampling frequency
    # in seconds
    meg_headers['RecordingDuration'] = round(raw.n_times / raw.info['sfreq'], 2)

    # get a list of Head Position Indicator coil frequencies
    try:
        chpi_freqs, ch_idx, chpi_codes = mne.chpi.get_chpi_info(info=raw.info)
        meg_headers['HeadCoilFrequency'] = str(chpi_freqs)
        meg_headers['ContinuousHeadLocalization'] = True

        chpi_amplitudes = mne.chpi.compute_chpi_amplitudes(raw)
        chpi_locs = mne.chpi.compute_chpi_locs(raw.info, chpi_amplitudes)
        head_pos = mne.chpi.compute_head_pos(raw.info, chpi_locs, verbose=True)

        # get the max movement among x, y and z positions in mm
        # originally in meters
        meg_headers['MaxMovement'] = max(max(head_pos[:, 4]),
                                         max(head_pos[:, 5]),
                                         max(head_pos[:, 6])) * 1000
    except ValueError:
        meg_headers['ContinuousHeadLocalization'] = False

    # channels information
    # name, type, units, sampling frequency

    # contains a list of channel names
    meg_headers['channel_names'] = raw.ch_names
    meg_headers['channel_types'] = [mne.channel_type(raw.info, idx)
                                    for idx in mne.pick_types(raw.info, include=raw.ch_names)]

    _unit2human = {'FIFF_UNIT_V': 'V',
                   'FIFF_UNIT_T': 'T',
                   'FIFF_UNIT_T_M': 'T/m',
                   'FIFF_UNIT_MOL': 'M',
                   'FIFF_UNIT_NONE': 'NA',
                   'FIFF_UNIT_CEL': 'C'}

    mne_units = [str(raw.info['chs'][i]['unit']).split(' ')[1].strip('()')
                 for i in range(len(raw.info['chs']))]
    meg_headers['channel_units'] = [_unit2human[i] for i in mne_units]

    # default mne values
    meg_headers['MEGCoordinateSystem'] = 'ElektaNeuromag'
    meg_headers['MEGCoordinateUnits'] = 'm'

    meg_headers['HeaderExtractionTool'] = 'mne-python'
    meg_headers['experimenter'] = raw.info['experimenter']

    with open(out_fname, "w") as f:
        json.dump(meg_headers, f, indent=4)
    if generate_properties_file_for_xnat:
        prop_file_name = out_fname + '.prop'
        generateFIFProppertiesFileForXnat(meg_headers, prop_file_name)

    return meg_headers


def edf_header_to_properties(edf_fname, json_template, out_fname,
                             generate_properties_file_for_xnat):
    """Call *fif_header_to_properties* on fif file and save json output.

    Parameters
    ----------
    fif_fname : str
        Raw fif File
    json_template: str
        BIDS template used to populate various fields.
    out_fname : str
        Output file name
        relative paths are saved relative to parent dir of fif_fname
        None will save to parent dir of fif_fname with default prefix
    generate_properties_file_for_xnat : flag
        Generates the properties file used by XNAT            
    """
    # load raw object
    raw = mne.io.read_raw_edf(edf_fname, verbose='ERROR')
    ecog_headers = json.loads(json_template)

    if raw.info['device_info']:
        ecog_headers['InstitutionName'] = raw.info['device_info']['site']
        ecog_headers['Manufacturer'] = raw.info['device_info']['type']
        ecog_headers['ManufacturersModelName'] = raw.info['device_info']['model']
        ecog_headers['DeviceSerialNumber'] = raw.info['device_info']['serial']

    ecog_headers['NumberOfChannels'] = raw.info['nchan']
    ecog_headers['SamplingFrequency'] = raw.info['sfreq']
    ecog_headers['PowerLineFrequency'] = raw.info['line_freq']
    ecog_headers['HardwareFilters'] = {
        'HighpassFilter': {'CutoffFrequency': [raw.info['highpass']]},
        'LowpassFilter': {'CutoffFrequency': [raw.info['lowpass']]}
        }

    # channel counts
    ecog_headers['EEGChannelCount'] = len(mne.pick_types(raw.info, eeg=True))
    ecog_headers['ECOGChannelCount'] = len(mne.pick_types(raw.info, ecog=True))
    ecog_headers['SEEGChannelCount'] = len(mne.pick_types(raw.info, seeg=True))
    ecog_headers['EOGChannelCount'] = len(mne.pick_types(raw.info, eog=True))
    ecog_headers['ECGChannelCount'] = len(mne.pick_types(raw.info, ecg=True))
    ecog_headers['EMGChannelCount'] = len(mne.pick_types(raw.info, emg=True))
    ecog_headers['MiscChannelCount'] = len(mne.pick_types(raw.info, misc=True))
    ecog_headers['TriggerChannelCount'] = len(mne.pick_types(raw.info, stim=True))

    # recording duration = total samples / sampling frequency in seconds
    ecog_headers['RecordingDuration'] = round(raw.n_times / raw.info['sfreq'], 2)

    # channels information
    # name, type, units, sampling frequency

    # contains a list of channel names
    ecog_headers['channel_names'] = raw.ch_names
    ecog_headers['channel_types'] = [mne.channel_type(raw.info, idx)
                                     for idx in mne.pick_types(raw.info,
                                     include=raw.ch_names)]

    _unit2human = {'FIFF_UNIT_V': 'V',
                   'FIFF_UNIT_T': 'T',
                   'FIFF_UNIT_T_M': 'T/m',
                   'FIFF_UNIT_MOL': 'M',
                   'FIFF_UNIT_NONE': 'NA',
                   'FIFF_UNIT_CEL': 'C'}

    mne_units = [str(raw.info['chs'][i]['unit']).split(' ')[1].strip('()')
                 for i in range(len(raw.info['chs']))]
    ecog_headers['channel_units'] = [_unit2human[i] for i in mne_units]

    ecog_headers['HeaderExtractionTool'] = 'mne-python'

    with open(edf_fname, 'rb') as fid:
        # directly read patient ID from the file
        ecog_headers['patient_id'] = fid.read(80).decode('latin-1').rstrip()[8:]

    with open(out_fname, "w") as f:
        json.dump(ecog_headers, f, indent=4)
    
    if generate_properties_file_for_xnat:
        prop_file_name = out_fname + '.prop'
        generateEDFProppertiesFileForXnat(ecog_headers, raw, prop_file_name)

    return ecog_headers


def run():
    """Run *et/fif/edf_header_to_properties* command."""
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
    parser.add_option("-v", "--verbose",  action="store_true", dest="verbose",
                      help="Print verbose messages.",
                      metavar="VERBOSE", default=False)
    
    options, args = parser.parse_args()

    if options.file is None:
        print('Raw file not provided.')
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

    generate_properties_file_for_xnat = False
    if options.xProperties:
        generate_properties_file_for_xnat = True

    if options.verbose:
        mne.set_log_level('DEBUG')
    else:
        mne.set_log_level('ERROR')

    ecog_template_json_string = '{"TaskName": "", "SamplingFrequency": "", "PowerLineFrequency": "", "SoftwareFilters": "", "HardwareFilters": {"HighpassFilter": {"CutoffFrequency": []}, "LowpassFilter": {"CutoffFrequency": []}}, "Manufacturer": "", "ManufacturersModelName": "", "TaskDescription": "", "Instructions": "", "CogAtlasID": "", "CogPOID": "", "InstitutionName": "", "InstitutionAddress": "", "DeviceSerialNumber": "", "ECOGChannelCount": "", "SEEGChannelCount": "", "EEGChannelCount": "", "EOGChannelCount": "", "ECGChannelCount": "", "EMGChannelCount": "", "MiscChannelCount": "", "TriggerChannelCount": "", "RecordingDuration": "", "RecordingType": "", "EpochLength": "", "SubjectArtefactDescription": "", "SoftwareVersions": "", "iEEGReference": "", "ElectrodeManufacturer": "", "ElectrodeManufacturersModelName": "", "iEEGGround": "", "iEEGPlacementScheme": "", "iEEGElectrodeGroups": "", "ElectricalStimulation": "", "ElectricalStimulationParameters": ""}'
    meg_template_json_string = '{"TaskName": "", "InstitutionName": "", "InstitutionAddress": "", "Manufacturer": "", "ManufacturersModelName": "", "SoftwareVersions": "", "TaskDescription": "", "Instructions": "", "CogPOID": "", "DeviceSerialNumber": "", "SamplingFrequency": "", "PowerLineFrequency": "", "DewarPosition": "", "SoftwareFilters": "", "DigitizedLandmarks": "", "DigitizedHeadPoints": "", "MEGChannelCount": "", "MEGREFChannelCount": "", "EEGChannelCount": "", "ECOGChannelCount": "", "SEEGChannelCount": "", "EOGChannelCount": "", "ECGChannelCount": "", "EMGChannelCount": "", "MiscChannelCount": "", "TriggerChannelCount": "", "RecordingDuration": "", "RecordingType": "", "EpochLength": "", "ContinuousHeadLocalization": "", "HeadCoilFrequency": "", "MaxMovement": "", "SubjectArtefactDescription": "", "AssociatedEmptyRoom": "", "EEGSamplingFrequency": "", "EEGPlacementScheme": "", "ManufacturersAmplifierModelName": "", "EEGReference": ""}'
    et_template_json_string = '{"InstitutionName": "", "InstitutionAddress": "", "Manufacturer": "", "ManufacturersModelName": "", "SoftwareVersions": "", "TaskDescription": "", "Instructions": "", "CogAtlasID": "", "CogPOID": "", "DeviceSerialNumber": "", "Firmware #": "", "SamplingFrequency": null, "SampleCoordinateUnit": "", "SampleCoordinateSystem": "", "EnvironmentCoordinates": [], "EventIdentifier": "", "RawSamples": "", "IncludedEyeMovementEvents": [], "DetectionAlgorithm": "", "Eye tracker camera": "", "Eye tracker lens": "", "StartMessage": "", "EndMessage": "", "KeyPressMessage": "", "CalibrationType": "", "CalibrationPosition": null, "CalibrationUnit": "", "MaximalCalibrationError": null, "AverageCalibrationError": null, "CalibrationList": [], "RecordedEye": "", "EyeCameraSettings": [], "FeatureDetectionSettings": [], "GazeMappingSettings": [], "DetectionAlgorithmSettings": [], "RawDataFilters": "", "ScreenSize": "", "ScreenResolution": "", "ScreenRefreshRate": "", "AOIDefinition": [], "PupilPositionType": "", "PupilFitMethod": ""}'

    # MEG data is always in FIF, ECOG data in EDF format
    if fname.split('.')[-1] in ['FIF', 'fif']:
        fif_headers = fif_header_to_properties(fname, meg_template_json_string, out_fname,
                                               generate_properties_file_for_xnat)
    elif fname.split('.')[-1] in ['EDF', 'edf']:
        edf_headers = edf_header_to_properties(fname, ecog_template_json_string, out_fname,
                                               generate_properties_file_for_xnat)
    elif fname.split('.')[-1] in ['asc', 'ASC']:
        et_config_dict = extract_asc_file_headers(fname, out_fname,
                                                  generate_properties_file_for_xnat=generate_properties_file_for_xnat)
    elif fname.split('.')[-1] in ['csv', 'CSV']:
        et_config_dict = extract_csv_file_headers(fname, out_fname,
                                                  generate_properties_file_for_xnat=generate_properties_file_for_xnat)
    else:
        raise ValueError('%s does not seem to be a FIF or an EDF file.' % fname)

    # print('ET/FIF/EDF File header extraction complete. Output written to %s' % out_fname)


is_main = (__name__ == '__main__')
if is_main:
    run()

