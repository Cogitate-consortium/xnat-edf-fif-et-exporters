# xnat-edf-fif-et-exporters
Backend scripts that export headers that will be used by the XNAT edf-fif-et importer plugin. 

The scripts have to be installed on the XNAT server such that they are available in the `PATH` of the XNAT user. The edf-fif-et-importer plugin installed on XNAT uses these scripts for extracting header information for the additional data types of MEEG (FIF files), ECOG (EDF Files) and Eye Tracking (ASC and CSV files produced by EYELINK and TOBII eye trackers).

For a fresh setup the below instructions may be followed,

1. Connect to the xnat server (via ssh) and switch to the xnat user.
2. Setup a working python environment with mne-python installed following the instructions provided [here](https://mne.tools/stable/install/manual_install.html). Ensure that the environment is correctly activated.
3. Make the scripts executable with `chmod +x et_importer.py` and `chmod +x get_meg_ecog_headers.py`.
4. Move them to `~/.local/bin` direectory and ensure that this directory is added to the PATH by adding this line `PATH="$HOME/.local/bin:$PATH"` to the `.bashrc` file.
5. Test them locally to check if properties files are being correctly generated. See below for more details.

## Extracting file headers
Scripts to extract header information from raw ET/MEEG/ECOG data and make them available to XNAT.

The scripts here will identify and parse out relevant headers to be modeled on XNAT. Mostly BIDS fields are chosen where available. While MEG/EDF headers are relatively standard, both raw file formats headers fields for eye tracking data vary significantly across devices. Here we use Eyelink data converted to asc files and and Tobii eye tracker data converted to csv.

Note: The scripts may have to be modified based on the headers of interest.

### Requirements
Requires [mne-python](https://mne.tools/stable/install/index.html) package.

`pip install mne`

### Usage

1. Install mne
2. Run get_et_meg_ecog_headers.py -f <file_name> -x
3. Should generate properties `*.prop` file with given header values.

### Examples

#### MEG

mpg:isAnon=True
bids:samplingFrequency=1000.0
bids:powerLineFrequency=50.0
bids:dewarPosition=68.0
bids:softwareFilters_default_highpass=0.10000000149011612
bids:softwareFilters_default_lowpass=330.0
bids:digitizedLandmarks=True
bids:digitizedHeadPoints=True
bids:megChannelCount=306
bids:recordingDuration=866.0
bids:numberOfChannels=341
bids:megCoordinateSystem=ElektaNeuromag

#### ECOG
mpg:isAnon=True
bids:samplingFrequency=1024.0
bids:numberOfChannels=144
bids:hardwareFilters_HighpassFilter_CutoffFrequency=[0.0]
bids:hardwareFilters_LowpassFilter_CutoffFrequency=[512.0]
bids:recordingDuration=938.0

#### ET

mpg:isAnon=True
bids:Manufacturer=EYELINK
bids:SamplingFrequency=1000
mpg:DistanceToScreen=140.0 cms
mpg:ScreenSize=39.0 22.0 cms
mpg:RecordedEye=Both

#### Abbreviations

- ET: Eye Tracking
- MEG: Magnetoencephalography
- ECOG: Electrocorticography
- FIF: Functional Imaging Format (MEG files are stored this way)
- EDF: European Data Format (ECOG files are stored this way)

## Notes
Names of the import handlers for the additional datatypes for XNAT REST calls via `/data/services/import` - **MEG/EEG** and **ECOG** and **EyeTracker** respectively.

## Repository Contributors
- [Adeel Ansari](https://github.com/adeel-ansari)
- [Mohana Ramaratnam](https://github.com/mohanakannan9)
- [Praveen Sripad](https://github.com/pravsripad)

## Acknowledgements
This project was made possible through the support of a grant from Templeton World Charity Foundation, Inc. The tool has been developed as a part of the [ARC-COGITATE](https://www.arc-cogitate.com/) project.
