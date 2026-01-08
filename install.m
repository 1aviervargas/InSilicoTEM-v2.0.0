%Add Matlab Path:
addpath('dipimage_2.9_lin64/dip/common/mlv7_6/diplib')
addpath('dipimage_2.9_lin64/dip/common/mlv7_6/dipimage_mex')
addpath('dipimage_2.9_lin64/dip/common/dipimage/aliases')
addpath('dipimage_2.9_lin64/dip/common/dipimage')

%With subfolders
addpath('/home/jvargas/Software/InSilicoTEM-v2.0.0');

dip_initialise
dipsetpref('imagefilepath','/herebeimages')

addpath(genpath('M4I-nanoscopy-InSilicoTEM-c3606da'))

