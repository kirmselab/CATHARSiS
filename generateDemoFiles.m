%% This script generates a set of demo files for testing CATHARSiS
% ("Calcium transient detection harnessing spatial similarity").
% ----------------------------------------------------------------------
% Reference:
% Jürgen Graf, Vahid Rahmati, Myrtill Majoros, Otto W. Witte, Christian Geis,
% Stefan J. Kiebel, Knut Holthoff, Knut Kirmse. (2021) Network instability
% dynamics drive a transient bursting period in the developing hippocampus
% in vivo. bioRxiv. doi.org/10.1101/2021.05.28.446133
% ======================================================================
% Use:
% (1) >> generateDemoFiles
% ======================================================================
% Help:
% Please see "CATHARSiS_User_Manual.pdf".
% ======================================================================
%% Clean-up
clear variables; close all;
%% Parameters (two cells will be simulated)
inR = 9; % inner radius of simulated cells [pixels]
outR = 15; % outer radius of simulated cells [pixels]
distrange = 13; % a scalar or vector specifying the offset(s) between cells [pixels]
% Note: If distrange = 11:10:31, then
% • three TIF files (simulated Ca2+ imaging data)
% • three TXT files (containing the corresponding ROI corrdinates)
% • one TXT containing drift intervals (here, an empty TXT file, i.e. no drift)
% • one TXT specifying file onset frames
% will be generated.
%% Simulate data & save demo files
disp('Generating simulated data. Please wait ...')
timestamp = datestr(now, 'yymmdd_HHMMSS');
dir_out = ['Demo_', timestamp];
mkdir(dir_out);
sim_rings(dir_out, inR, outR, distrange);
dlmwrite([dir_out, '\drift.txt'], [], '\t'); % no drift periods defined
dlmwrite([dir_out, '\onsets.txt'], 1, '\t'); % only the first frame represents a file onset
disp(['Sample files were generated and stored in the following folder: ', dir_out, '.'])