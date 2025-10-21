# Reference

This analysis code is adapted from code used in: 

Banks A, Eldin Abdelaal A, Salcudean S. Head motion-corrected eye gaze tracking with the da Vinci surgical system. Int J Comput Assist Radiol Surg. 2024 Jul;19(7):1459-1467. doi: 10.1007/s11548-024-03173-4. Epub 2024 Jun 18. PMID: 38888820. 

[https://pubmed.ncbi.nlm.nih.gov/38888820/](https://pubmed.ncbi.nlm.nih.gov/38888820/)

# Eye Gaze Calibration Steps

**Done when collecting data**

- Do an 8-point calibration to fit the standard polynomial regressor.
- Do 5x1-point calibration to fit the corner-contingent head compensation model. During this calibration, get the user to look at a single dot in the center of the screen, then lift there head, and then return their head and look at the dot again. Do this five times.

# Installation
1. Python version 3.10.16
2. MATLAB version R2023b
3. Run "pip install -r cchc_requirements.txt"

# Data Structure After All Scripts are Run
- CCHC_Gaze/
- GazeData/ (or other root directory at same level as CCHC_Gaze code directory)
  - surgery1/
    - eyeVideo_MM-DD-YYYY_HR-MM-SC.avi
    - gazelog_MM-DD-YYYY_HR-MM-SC.txt
    - eyeCornerData_MM-DD-YYYY_HR-MM-SC.csv
    - mergedData_MM-DD-YYYY_HR-MM-SC.csv
    - CCHC_POG_MM-DD-YYYY_HR-MM-SC.csv
    - ...
    - calibration/
      - eyeVideo_MM-DD-YYYY_HR-MM-SC.avi
      - gazelog_MM-DD-YYYY_HR-MM-SC.txt
      - eyeCornerData_MM-DD-YYYY_HR-MM-SC.csv
      - mergedData_MM-DD-YYYY_HR-MM-SC.csv
      - CCHC_Calibration.mat
      - ...
  - surgery2/
  - surgery3/
  - ...

**Notes:
1. The "calibration/" subdirectory is where the eye recordings are for computing the calibration params. Outside of this subdirectory in the "surgery<X>/" subdirectory is the overall data collected (other than the calibration data) and is the data that the calibration parameters are applied to to estimate the POG.
2. "eyeVideo_MM-DD-YYYY_HR-MM-SC.avi" is used to determine the eye corners and create a corresponding "eyeCornerData_MM-DD-YYYY_HR-MM-SC.csv".
3. "gazelog_MM-DD-YYYY_HR-MM-SC.txt" is used to get the eye gaze tracking parameters and is then merged with the corresponding "eyeCornerData_MM-DD-YYYY_HR-MM-SC.csv" to give the "mergedData_MM-DD-YYYY_HR-MM-SC.csv"
4. CCHC_Calibration.mat is where the calibration parameters are stored.
5. CCHC_POG_MM-DD-YYYY_HR-MM-SC.csv is the estimated POG (with columns for each eye left/right in x/y screen coordinates, and columns indicating whether the CCHC parameters were used) from the data in the "surgery1/" subdirectory**

  
# Running the Code

**Note that this is development code and may need to be modified if you are using different paths.**

1. Run _EyeCornerDetector.py_ on all the data in the "calibration/" directories in each "surgery<x>/" directory. This creates the "eyeCornerData_MM-DD-YYYY_HR-MM-SC.csv" in the "surgery<x>/calibration/" directory. **Set the "IS_ON_CALIB" parameter=True on line 41 to run the eye corner detection on the calibration data**. Change "data_root" on line 39 if your root directory isn't "GazeData/"
2. Run _EyeCornerDetector.py_ on all the data to estimate the POG on in each "surgery<x>/" directory (not in the "calibration/" directories). This creates the "eyeCornerData_MM-DD-YYYY_HR-MM-SC.csv" in the "surgery<x>/" directory. **Set the "IS_ON_CALIB" parameter=False on line 41 to run the eye corner detection on the non-calibration (estimation) data**. Change "data_root" on line 39 if your root directory isn't "GazeData/". 
3. Run the _CCHCCalibration.m_ script. This runs through each "calibration/" subdirectory in each "surgery<x>/" directory, creates the "mergedData_MM-DD-YYYY_HR-MM-SC.csv" file (with eye gaze parameters and eye corners synchronized and merged).
4. Run the _CCHCGazeCompensation.m_ script. This runs through each "surgery<x>/" directory, uses the calibration parameters previously computed to estimate the POG for a given gazelog_MM-DD-YYYY_HR-MM-SC.txt. It creates the "CCHC_POG_MM-DD-YYYY_HR-MM-SC.csv" file (with columns for each eye left/right in x/y screen coordinates, and columns indicating whether the CCHC parameters were use).

# License

This work is protected under the [Apache-2.0](https://www.apache.org/licenses/LICENSE-2.0) license.
