#!/usr/bin/env python
# author: Praveen Sripad <praveen.sripad@ae.mpg.de>

"""Extract FIF/EDF file header into a properties file.

--------
.. code-block:: console
    $ get_fif_bids_headers.py -f sample_audvis_raw.fif -b sub-01_task-FullExample_acq-CTF_run-1_proc-sss_meg.json -o headers.json [-x]
"""

import sys
import os.path as op
import json
import mne


def generateFIFProppertiesFileForXnat(meg_headers, out_fname):
        print('Generating Properties File for XNAT ' + out_fname)

        isAnon = False
        if meg_headers['experimenter'] == 'XXX':
            isAnon = True

        out_file = open(r'%s'% out_fname,"w")
        out_file.write('isAnon=%s\n'% isAnon)
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
        print('Generating Properties File for XNAT')

        isAnon = False
        if ecog_headers['patient_id'] == 'X X X X_X':
            isAnon = True

        print(raw.info['subject_info'] )
        out_file = open(r'%s'% out_fname,"w")
        out_file.write('isAnon=%s\n'% isAnon)
        out_file.write('bids_frequency=%s\n'% ecog_headers['SamplingFrequency'])
        out_file.write('bids_numberOfSignals=%s\n'% ecog_headers['NumberOfChannels'])
        #out_file.write('localSubjectIdentification=%s\n'% raw.info['patient']['id'])
        #out_file.write('identificationCode=%s\n'% raw.info['version'])
        #out_file.write('localRecordingIdentification=%s\n'% raw.info['recording'])
        #out_file.write('numberOfDataRecords=%s\n'% raw.info['datarecords'])
        #out_file.write('dataRecord_duration=%s\n'% raw.info['duration'])
        out_file.write('bids_hardwareFilters_HighpassFilter_CutoffFrequency=%s\n'% ecog_headers['HardwareFilters']['HighpassFilter']['CutoffFrequency'])
        out_file.write('bids_hardwareFilters_LowpassFilter_CutoffFrequency=%s\n'% ecog_headers['HardwareFilters']['LowpassFilter']['CutoffFrequency'])
        out_file.write('bids_recordingDuration=%s\n'% ecog_headers['RecordingDuration'])
        out_file.close()
    

def fif_header_to_properties(fif_fname, json_template, out_fname, generate_properties_file_for_xnat):
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

    print('Writing output to ')
    with open(out_fname, "w") as f:
        json.dump(meg_headers, f, indent=4)
    if generate_properties_file_for_xnat:
        prop_file_name = out_fname + '.prop'
        generateFIFProppertiesFileForXnat(meg_headers, prop_file_name)

def edf_header_to_properties(edf_fname, json_template, out_fname, generate_properties_file_for_xnat):
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
    ecog_headers['HardwareFilters'] = {'HighpassFilter': {'CutoffFrequency': [raw.info['highpass']]},
                                       'LowpassFilter': {'CutoffFrequency': [raw.info['lowpass']]}}

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
                                    for idx in mne.pick_types(raw.info, include=raw.ch_names)]

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

    print('Writing output to ')
    with open(out_fname, "w") as f:
        json.dump(ecog_headers, f, indent=4)
    
    if generate_properties_file_for_xnat:
        prop_file_name = out_fname+ '.prop'
        generateEDFProppertiesFileForXnat(ecog_headers, raw, prop_file_name)


def run():
    """Run *fif_header_to_properties* command."""
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
    
    ecog_template_json_string = '{"TaskName": "", "SamplingFrequency": "", "PowerLineFrequency": "", "SoftwareFilters": "", "HardwareFilters": {"HighpassFilter": {"CutoffFrequency": []}, "LowpassFilter": {"CutoffFrequency": []}}, "Manufacturer": "", "ManufacturersModelName": "", "TaskDescription": "", "Instructions": "", "CogAtlasID": "", "CogPOID": "", "InstitutionName": "", "InstitutionAddress": "", "DeviceSerialNumber": "", "ECOGChannelCount": "", "SEEGChannelCount": "", "EEGChannelCount": "", "EOGChannelCount": "", "ECGChannelCount": "", "EMGChannelCount": "", "MiscChannelCount": "", "TriggerChannelCount": "", "RecordingDuration": "", "RecordingType": "", "EpochLength": "", "SubjectArtefactDescription": "", "SoftwareVersions": "", "iEEGReference": "", "ElectrodeManufacturer": "", "ElectrodeManufacturersModelName": "", "iEEGGround": "", "iEEGPlacementScheme": "", "iEEGElectrodeGroups": "", "ElectricalStimulation": "", "ElectricalStimulationParameters": ""}'
    meg_template_json_string = '{"TaskName": "", "InstitutionName": "", "InstitutionAddress": "", "Manufacturer": "", "ManufacturersModelName": "", "SoftwareVersions": "", "TaskDescription": "", "Instructions": "", "CogPOID": "", "DeviceSerialNumber": "", "SamplingFrequency": "", "PowerLineFrequency": "", "DewarPosition": "", "SoftwareFilters": "", "DigitizedLandmarks": "", "DigitizedHeadPoints": "", "MEGChannelCount": "", "MEGREFChannelCount": "", "EEGChannelCount": "", "ECOGChannelCount": "", "SEEGChannelCount": "", "EOGChannelCount": "", "ECGChannelCount": "", "EMGChannelCount": "", "MiscChannelCount": "", "TriggerChannelCount": "", "RecordingDuration": "", "RecordingType": "", "EpochLength": "", "ContinuousHeadLocalization": "", "HeadCoilFrequency": "", "MaxMovement": "", "SubjectArtefactDescription": "", "AssociatedEmptyRoom": "", "EEGSamplingFrequency": "", "EEGPlacementScheme": "", "ManufacturersAmplifierModelName": "", "EEGReference": ""}'
    et_template_json_string = '{"InstitutionName": "", "InstitutionAddress": "", "Manufacturer": "", "ManufacturersModelName": "", "SoftwareVersions": "", "TaskDescription": "", "Instructions": "", "CogAtlasID": "", "CogPOID": "", "DeviceSerialNumber": "", "Firmware #": "", "SamplingFrequency": null, "SampleCoordinateUnit": "", "SampleCoordinateSystem": "", "EnvironmentCoordinates": [], "EventIdentifier": "", "RawSamples": "", "IncludedEyeMovementEvents": [], "DetectionAlgorithm": "", "Eye tracker camera": "", "Eye tracker lens": "", "StartMessage": "", "EndMessage": "", "KeyPressMessage": "", "CalibrationType": "", "CalibrationPosition": null, "CalibrationUnit": "", "MaximalCalibrationError": null, "AverageCalibrationError": null, "CalibrationList": [], "RecordedEye": "", "EyeCameraSettings": [], "FeatureDetectionSettings": [], "GazeMappingSettings": [], "DetectionAlgorithmSettings": [], "RawDataFilters": "", "ScreenSize": "", "ScreenResolution": "", "ScreenRefreshRate": "", "AOIDefinition": [], "PupilPositionType": "", "PupilFitMethod": ""}'

    # MEG data is always in FIF, ECOG data in EDF format
    if fname.split('.')[-1] in ['FIF', 'fif']:
        fif_header_to_properties(fname, meg_template_json_string, out_fname, generate_properties_file_for_xnat)
    elif fname.split('.')[-1] in ['EDF', 'edf']:
        edf_header_to_properties(fname, ecog_template_json_string, out_fname, generate_properties_file_for_xnat)
    else:
        raise ValueError('%s does not seem to be a FIF or an EDF file.' % fname)
    
    print('FIF/EDF File header extraction complete. Output written to %s' % out_fname)


is_main = (__name__ == '__main__')
if is_main:
    run()
