% runLUT.m
% runs LUT model script

addpath kdtree

% set no. 1
calibFile = 'E:\Fld_Res\Calib\BCLBL021.FLD';
lutFile = 'LUT2307.mat';
datFile = 'E:\Fldat\FLDTL021.FLD';
resFile = 'E:\Velfld\VELLUTL021.FLD';


% Options
opt = struct;

% General values:
opt.numArrays = 5; % number of arrays in the probe
opt.numWires  = 4; % number of wires per array

% Calibration options

opt.lenCalibFile = 490; % number of calibration points
opt.testMode = 1; % use 1 if the self-test is needed, with the same data as calibrated
opt.tosave = 1; % use 1 first time the lutFile is generated or need to be regenerated

opt.numelVel = 50;    % 50 different velocity interpolations
opt.dy = 1; opt.dz = 1; % in degrees, smallest resolved angle in LUT, smaller value means longer run and memory consumption
opt.du = .5; % m/s: 9 - 15 m/s = 6/.3 is about 20 different levels
opt.dy2 = .2; opt.dz2 = .2; % in degrees, smallest resolved angle in LUT, smaller value means longer run and memory consumption
opt.du2 = .2; % m/s: 9 - 15 m/s = 6/.3 is about 20 different levels



% RAW data reading options:
opt.skipLength =  0;
opt.readLength = 200; % 490*8 values

% Writing results options
opt.write = 0; % if you do not want to write the following resFile, = zero


% lut_calibration(calibFile,datFile,resFile,lutFile,opt)
% 
% set no 2

calibFile   = 'E:\Fld_Res\Calib\BCLBL117.FLD';
lutFile      = 'E:\Fld_Res\Calib\BCLBL117.MAT';
datFile     = 'E:\Fldat\Fldtl117.fld';
resFile     = 'E:\Velfld\VELLUTL117.FLD';

lut_calibration(calibFile,datFile,resFile,lutFile,opt)

% if it was successful run:
opt.tosave = 0;
opt.testMode = 0;
lut_calibration(calibFile,datFile,resFile,lutFile,opt)

% and if this one worked out, then the full run
opt.write = 1;
opt.readLength = 2^16;
lut_calibration(calibFile,datFile,resFile,lutFile,opt)


