%  eprload_jeol       Reads in data from JEOL files
%
%  [Data,Abscissa,Parameters] = eprload_jeol(FileName)
%
%  JEOL file formats for JES-FA and JES-X3 are supported.
%  This implementation is based on official documentation provided by JEOL.
%
%  Input:
%    FileName    file name of JEOL data file
%
%  Output:
%    Data        1D or 2D data, real or complex
%    Abscissa    vector of horizontal axis values
%    Parameters  structure containing all parameters
%                extracted from file
%
%  The file format is (C) 2015 by JEOL RESONANCE INC.
%  The implementation is (C) 2015 Stefan Stoll

function [Data,Abscissa,Parameters] = eprload_jeol(FileName)

h = fopen(FileName,'r','ieee-le');
if h<0, error('Could not open %s.',FileName); end

% [A] HEADER block (0x0000-0x00143, dataicon 0x00144-0x01143)
%====================================================================
jeol_header.processType = readcstring(h,16);
jeol_header.filename = readcstring(h,64);
fread(h,1,'int32'); % reserved
jeol_header.auxNum = fread(h,1,'int32');
jeol_header.dataSet = fread(h,1,'int32');
jeol_header.dataLength = fread(h,1,'int32');
jeol_header.dataKind = readcstring(h,15);
jeol_header.dispPart = fread(h,1,'*char');
jeol_header.xOffset = fread(h,1,'float32');
jeol_header.xRange = fread(h,1,'float32');
jeol_header.xUnit = readcstring(h,16);
jeol_header.yMin = fread(h,1,'float32');
jeol_header.yMax = fread(h,1,'float32');
jeol_header.yUnit = readcstring(h,16);
jeol_header.zPoint = fread(h,1,'float32');
jeol_header.zUnit = readcstring(h,16);
jeol_header.xStart = fread(h,1,'float32');
jeol_header.xStop = fread(h,1,'float32');
jeol_header.yStart = fread(h,1,'float32');
jeol_header.yStop = fread(h,1,'float32');
jeol_header.zStart = fread(h,1,'float32');
jeol_header.zStop = fread(h,1,'float32');
fread(h,88,'*char'); % reserved
jeol_header.waveIndex = fread(h,1,'int32');
fread(h,28,'*char'); % reserved
jeol_endor = strcmp(jeol_header.processType,'endor');
fread(h,4096,'int8'); % data icon bitmap

% [B] GENERAL block (0x01144-0x018AB)
%====================================================================
fread(h,16,'*char'); % reserved
jeol_general.system = readcstring(h,32);
jeol_general.version = readcstring(h,16);
jeol_general.date = readcstring(h,32);
jeol_general.title = readcstring(h,128);
jeol_general.sample = readcstring(h,128);
jeol_general.comment = readcstring(h,128);
fread(h,224,'*char'); % reserved

jeol_general.previousINFO_filename = readcstring(h,64);
jeol_general.previousINFO_dir = readcstring(h,512);
fread(h,1,'int32'); % not used
jeol_general.previousINFO_processType = readcstring(h,16);

jeol_general.originalINFO_filename = readcstring(h,64);
jeol_general.originalINFO_dir = readcstring(h,512);
fread(h,1,'int32'); % not used
jeol_general.originalINFO_processType = readcstring(h,16);

% [C] SPECTROMETER block (0x018AC-0x0205B)
%====================================================================
jeol_spectrometer.sType = readcstring(h,12);
jeol_spectrometer.SpectrometerFlag = fread(h,1,'int32');
jeol_spectrometer.magnet = readcstring(h,16);
jeol_spectrometer.water = readcstring(h,16);
jeol_spectrometer.magPower = readcstring(h,16);
jeol_spectrometer.magCurrent = readcstring(h,16);
jeol_spectrometer.centerField = readcstring(h,16);
jeol_spectrometer.sweepWidFin = readcstring(h,16);
jeol_spectrometer.sweepWidCor = readcstring(h,16);
jeol_spectrometer.sweepPulse = readcstring(h,16);
jeol_spectrometer.pulseNumber = readcstring(h,16);
jeol_spectrometer.sweepFTime = readcstring(h,16);
jeol_spectrometer.sweepBTime = readcstring(h,16);
jeol_spectrometer.scanTime = readcstring(h,16);
jeol_spectrometer.swCont = readcstring(h,16);
jeol_spectrometer.sPosition = readcstring(h,16);
jeol_spectrometer.ePosition = readcstring(h,16);
jeol_spectrometer.pPosition = readcstring(h,16);
jeol_spectrometer.cPosition = readcstring(h,16);
jeol_spectrometer.lfsBlock = readcstring(h,16);
jeol_spectrometer.xPosition = readcstring(h,16);
fread(h,112,'*char'); % reserved
jeol_spectrometer.modFreq = readcstring(h,16);
jeol_spectrometer.modWidFin = readcstring(h,16);
jeol_spectrometer.modWinCor = readcstring(h,16);
jeol_spectrometer.phase = readcstring(h,16);
jeol_spectrometer.Amp1_rcvMode = readcstring(h,16);
jeol_spectrometer.Amp1_phasePlus = readcstring(h,16);
jeol_spectrometer.Amp1_amplitudeFin = readcstring(h,16);
jeol_spectrometer.Amp1_amplitudeCor = readcstring(h,16);
jeol_spectrometer.Amp1_timeConstant = readcstring(h,16);
jeol_spectrometer.Amp1_zero = readcstring(h,16);
jeol_spectrometer.Amp2_rcvMode = readcstring(h,16);
jeol_spectrometer.Amp2_phasePlus = readcstring(h,16);
jeol_spectrometer.Amp2_amplitudeFin = readcstring(h,16);
jeol_spectrometer.Amp2_amplitudeCor = readcstring(h,16);
jeol_spectrometer.Amp2_timeConstant = readcstring(h,16);
jeol_spectrometer.Amp2_zero = readcstring(h,16);
jeol_spectrometer.amplock = readcstring(h,16);
jeol_spectrometer.sourceCH2 = readcstring(h,8);
fread(h,120,'*char'); % reserved
jeol_spectrometer.uType = readcstring(h,8);
jeol_spectrometer.uFreq = readcstring(h,16);
jeol_spectrometer.uFreqUnit = readcstring(h,8);
jeol_spectrometer.uPower = readcstring(h,16);
jeol_spectrometer.uPowerUnit = readcstring(h,8);
jeol_spectrometer.uPhase = readcstring(h,16);
jeol_spectrometer.uCoupling = readcstring(h,16);
jeol_spectrometer.slavePower = readcstring(h,16);
jeol_spectrometer.slavePhase = readcstring(h,16);
jeol_spectrometer.uRef = readcstring(h,16);
jeol_spectrometer.autoTuneMode = readcstring(h,16);
jeol_spectrometer.autoTune = readcstring(h,16);
jeol_spectrometer.afc = readcstring(h,8);
jeol_spectrometer.afcclk = readcstring(h,16);
jeol_spectrometer.mod = readcstring(h,8);
jeol_spectrometer.gunPower = readcstring(h,8);
jeol_spectrometer.ref = readcstring(h,8);
jeol_spectrometer.att_30dB = readcstring(h,8);
jeol_spectrometer.afcbalance = readcstring(h,8);
jeol_spectrometer.uDetM = readcstring(h,8);
jeol_spectrometer.uBalM = readcstring(h,8);
jeol_spectrometer.afcphase = readcstring(h,8);
jeol_spectrometer.uStab = readcstring(h,8);
jeol_spectrometer.shflock = readcstring(h,16);
jeol_spectrometer.uFreqCal = readcstring(h,16);
jeol_spectrometer.uFreqCalUnits = readcstring(h,8);
jeol_spectrometer.RotationStepAngle = readcstring(h,10);
jeol_spectrometer.RotationAngle = readcstring(h,10);
jeol_spectrometer.RotationZero = readcstring(h,10);
jeol_spectrometer.motorLow = readcstring(h,40);
jeol_spectrometer.motorUp = readcstring(h,40);
fread(h,18,'*char'); % reserved
jeol_spectrometer.acqPoint = readcstring(h,8);
jeol_spectrometer.dataLength = readcstring(h,8);
jeol_spectrometer.acqPulsMode = readcstring(h,8);
jeol_spectrometer.enType = readcstring(h,16);
jeol_spectrometer.ch1_eLeftFreq = readcstring(h,16);
jeol_spectrometer.ch1_eRightFreq = readcstring(h,16);
jeol_spectrometer.ch1_eCurrentFreq = readcstring(h,16);
jeol_spectrometer.ch1_eFunit = readcstring(h,8);
jeol_spectrometer.ch1_ePower = readcstring(h,16);
jeol_spectrometer.ch1_eDB = readcstring(h,12);
jeol_spectrometer.ch1_ePunit = readcstring(h,8);
jeol_spectrometer.ch1_eSweepTime = readcstring(h,12);
jeol_spectrometer.ch1_eModWidth = readcstring(h,8);
jeol_spectrometer.ch2_eLeftFreq = readcstring(h,16);
jeol_spectrometer.ch2_eRightFreq = readcstring(h,16);
jeol_spectrometer.ch2_eCurrentFreq = readcstring(h,16);
jeol_spectrometer.ch2_eFunit = readcstring(h,8);
jeol_spectrometer.ch2_ePower = readcstring(h,16);
jeol_spectrometer.ch2_eDB = readcstring(h,12);
jeol_spectrometer.ch2_ePunit = readcstring(h,8);
jeol_spectrometer.ch2_eSweepTime = readcstring(h,12);
jeol_spectrometer.ch2_eModWidth = readcstring(h,8);
fread(h,8,'*char'); % reserved
jeol_spectrometer.vtType = readcstring(h,16);
jeol_spectrometer.temperature = readcstring(h,16);
jeol_spectrometer.tempStart = readcstring(h,16);
jeol_spectrometer.tempEnd = readcstring(h,16);
jeol_spectrometer.tempStep = readcstring(h,16);
jeol_spectrometer.tempUnit = readcstring(h,8);
jeol_spectrometer.tempStatus = readcstring(h,8);
jeol_spectrometer.tempError = readcstring(h,16);
jeol_spectrometer.readyTime = readcstring(h,16);
jeol_spectrometer.okTempRange = readcstring(h,16);
jeol_spectrometer.tempControl = readcstring(h,8);
jeol_spectrometer.tempLock = readcstring(h,16);
jeol_spectrometer.endorLock = readcstring(h,16);
jeol_spectrometer.selfCheck = readcstring(h,8);
jeol_spectrometer.eCenterFreq = readcstring(h,16);
jeol_spectrometer.eSweepFreq = readcstring(h,16);
fread(h,48,'*char');
jeol_spectrometer.MarkerPos1 = readcstring(h,16);
jeol_spectrometer.MarkerPos2 = readcstring(h,16);
jeol_spectrometer.controlFA = readcstring(h,16);
fread(h,112,'*char');

% [D] GENERATOR block (CW: 0x0205C-0x21D7, ENDOR: 0x0205C-0x0250F)
%====================================================================
jeol_generator = struct;
if ~jeol_endor
  jeol_generator.date = readcstring(h,32);
  jeol_generator.accumulation = fread(h,1,'int32');
  jeol_generator.preAccumulation = fread(h,1,'int32');
  jeol_generator.delayTime = fread(h,1,'int32');
  jeol_generator.intervalTime = fread(h,1,'int32');
  jeol_generator.index = fread(h,1,'int32');
  jeol_generator.repetition = fread(h,1,'int32');
  fread(h,4,'int32');
  jeol_generator.sampleTime = fread(h,1,'int32');
  fread(h,4,'int32');
  jeol_generator.accumWave = readcstring(h,8);
  jeol_generator.baselineMode = readcstring(h,8);
  jeol_generator.sampleMode = readcstring(h,8);
  jeol_generator.sigTrig = readcstring(h,8);
  jeol_generator.refFile = readcstring(h,64);
  jeol_generator.paraFile = readcstring(h,64);
  jeol_generator.chMode = readcstring(h,16);
  fread(h,112,'*char');
elseif jeol_endor
  jeol_generator.date = readcstring(h,32);
  jeol_generator.accumulation = fread(h,1,'int32');
  jeol_generator.accumWave = readcstring(h,8);
  jeol_generator.baselineMode = readcstring(h,8);
  jeol_generator.chMode = readcstring(h,16);
  jeol_generator.MagnetPosition = readcstring(h,16);
  jeol_generator.ch1_eLeftFreq = readcstring(h,16);
  jeol_generator.ch1_eRightFreq = readcstring(h,16);
  jeol_generator.ch1_eCurrentFreq = readcstring(h,16);
  jeol_generator.ch1_eFunit = readcstring(h,8);
  jeol_generator.ch1_ePower = readcstring(h,16);
  jeol_generator.ch1_eDB = readcstring(h,12);
  jeol_generator.ch1_ePunit = readcstring(h,8);
  jeol_generator.ch1_eSweepTime = readcstring(h,12);
  jeol_generator.ch1_eModWidth = readcstring(h,16);
  fread(h,248,'*char');
  jeol_generator.ch2_eLeftFreq = readcstring(h,16);
  jeol_generator.ch2_eRightFreq = readcstring(h,16);
  jeol_generator.ch2_eCurrentFreq = readcstring(h,16);
  jeol_generator.ch2_eFunit = readcstring(h,8);
  jeol_generator.ch2_ePower = readcstring(h,16);
  jeol_generator.ch2_eDB = readcstring(h,12);
  jeol_generator.ch2_ePunit = readcstring(h,8);
  jeol_generator.ch2_eSweepTime = readcstring(h,12);
  jeol_generator.ch2_eModWidth = readcstring(h,16);
  fread(h,248,'*char');
  jeol_generator.eCenterFreq = readcstring(h,16);
  jeol_generator.eSweepFreq =  readcstring(h,16);
  fread(h,352,'*char');
end

% [E] PROCESS PARA block
%====================================================================
% (CW start: 0x021D8; ENDOR start: 0x02510; length: 0x0000-0x0343)
jeol_processpara = struct;
jeol_processpara.res_q = fread(h,1,'float32');
jeol_processpara.micro_hz = fread(h,1,'float32');
jeol_processpara.efect_pls = fread(h,1,'float32');
jeol_processpara.nike_hz = fread(h,1,'float32');
jeol_processpara.pulser = fread(h,1,'float32');
jeol_processpara.digitzer = fread(h,1,'float32');
jeol_processpara.resonator = fread(h,1,'float32');
jeol_processpara.fftCount = fread(h,1,'int32');
jeol_processpara.lastProcessTyp = readcstring(h,16);
jeol_processpara.std_spinFlag = fread(h,1,'int32');
jeol_processpara.std_marleft = fread(h,1,'float32');
jeol_processpara.std_marright = fread(h,1,'float32');
jeol_processpara.std_marintegval = fread(h,1,'float32');
jeol_processpara.std_raw_marintegval = fread(h,1,'float32');
jeol_processpara.std_sigleft = fread(h,1,'float32');
jeol_processpara.std_sigright = fread(h,1,'float32');
jeol_processpara.std_sigintegval = fread(h,1,'float32');
jeol_processpara.std_raw_sigintegval = fread(h,1,'float32');
jeol_processpara.std_spinnum = fread(h,1,'float32');
fread(h,4,'float32'); % reserved
jeol_processpara.spinFlag = fread(h,1,'int32');
jeol_processpara.marleft = fread(h,1,'float32');
jeol_processpara.marright = fread(h,1,'float32');
jeol_processpara.marintegval = fread(h,1,'float32');
jeol_processpara.raw_marintegval = fread(h,1,'float32');
jeol_processpara.sigleft = fread(h,1,'float32');
jeol_processpara.sigright = fread(h,1,'float32');
jeol_processpara.sigintegval = fread(h,1,'float32');
jeol_processpara.raw_sigintegval = fread(h,1,'float32');
jeol_processpara.spinnum = fread(h,1,'float32');
fread(h,4,'float32'); % reserved
jeol_processpara.xAreaM.xAreaFlag = fread(h,1,'int32');
jeol_processpara.xAreaM.left = fread(h,1,'float32');
jeol_processpara.xAreaM.right = fread(h,1,'float32');
jeol_processpara.xAreaM.peak = fread(h,1,'float32');
jeol_processpara.xAreaM.peakwidth = fread(h,1,'float32');
jeol_processpara.xAreaM.peakMin = fread(h,1,'float32');
jeol_processpara.xAreaM.peakMax = fread(h,1,'float32');
jeol_processpara.xAreaM.pleft = fread(h,1,'int32');
jeol_processpara.xAreaM.pright = fread(h,1,'int32');

jeol_processpara.xAreaS.xAreaFlag = fread(h,1,'int32');
jeol_processpara.xAreaS.left = fread(h,1,'float32');
jeol_processpara.xAreaS.right = fread(h,1,'float32');
jeol_processpara.xAreaS.peak = fread(h,1,'float32');
jeol_processpara.xAreaS.peakwidth = fread(h,1,'float32');
jeol_processpara.xAreaS.peakMin = fread(h,1,'float32');
jeol_processpara.xAreaS.peakMax = fread(h,1,'float32');
jeol_processpara.xAreaS.pleft = fread(h,1,'int32');
jeol_processpara.xAreaS.pright = fread(h,1,'int32');
for area = 1:8
  jeol_processpara.xArea(area).xAreaFlag = fread(h,1,'int32');
  jeol_processpara.xArea(area).left = fread(h,1,'float32');
  jeol_processpara.xArea(area).right = fread(h,1,'float32');
  jeol_processpara.xArea(area).peak = fread(h,1,'float32');
  jeol_processpara.xArea(area).peakwidth = fread(h,1,'float32');
  jeol_processpara.xArea(area).peakMin = fread(h,1,'float32');
  jeol_processpara.xArea(area).peakMax = fread(h,1,'float32');
  jeol_processpara.xArea(area).pleft = fread(h,1,'int32');
  jeol_processpara.xArea(area).pright = fread(h,1,'int32');
end
jeol_processpara.xLinem1.xLineFlag = fread(h,1,'int32');
jeol_processpara.xLinem1.pos = fread(h,1,'float32');
jeol_processpara.xLinem1.gvalue = fread(h,1,'float32');
jeol_processpara.xLinem1.point = fread(h,1,'int32');
jeol_processpara.xLinem2.xLineFlag = fread(h,1,'int32');
jeol_processpara.xLinem2.pos = fread(h,1,'float32');
jeol_processpara.xLinem2.gvalue = fread(h,1,'float32');
jeol_processpara.xLinem2.point = fread(h,1,'int32');
for m = 1:8
  jeol_processpara.xLine(m).xLineFlag = fread(h,1,'int32');
  jeol_processpara.xLine(m).pos = fread(h,1,'float32');
  jeol_processpara.xLine(m).gvalue = fread(h,1,'float32');
  jeol_processpara.xLine(m).point = fread(h,1,'int32');
end
jeol_processpara.aope = fread(h,4,'*char');
jeol_processpara.modFilename = readcstring(h,16);
jeol_processpara.yLine = fread(h,1,'int32');
jeol_processpara.calibFlag = fread(h,1,'int32');
fread(h,128,'*char');

% [F] CALCULATOR block
%====================================================================
jeol_calc = struct;
switch jeol_header.processType
  case 'yZero' % F1
    fread(h,160,'*char');
  case 'xZero' % F2
    jeol_calc.xzero = fread(h,1,'float32');
  case 'yGain' % F3
    jeol_calc.gain = fread(h,1,'float32');
  case 'xShift' % F4
    jeol_calc.shiftvalue = fread(h,1,'float32');
  case 'reverse' % F5
    jeol_calc.reversevalue = fread(h,1,'float32');
  case 'fill' % F6
    jeol_calc.fillingValue = fread(h,1,'float32');
    jeol_calc.oriEnd = fread(h,1,'int32');
    jeol_calc.length = fread(h,1,'int32');
  case 'window' % F7
    jeol_calc.parcent = fread(h,4,'float32');
    jeol_calc.xPoint = fread(h,4,'float32');
    jeol_calc.xValue = fread(h,4,'float32');
    jeol_calc.function = readcstring(h,32);
  case {'exp','log','power'} % F8
    % no data
  case {'fft','ifft'} % F9
    % no data
  case {'diff','inte'} % F10
    jeol_calc.inteFlag = fread(h,1,'int32');
    jeol_calc.order = fread(h,1,'int32');
    jeol_calc.points = fread(h,1,'float32');
    jeol_calc.pLeft = fread(h,1,'float32');
    jeol_calc.pRight = fread(h,1,'float32');
    jeol_calc.xLeft = fread(h,1,'float32');
    jeol_calc.xRight = fread(h,1,'float32');
    jeol_calc.baseline = fread(h,1,'float32');
  case 'smooth' % F11
    jeol_calc.ummy = readcstring(h,16);
  case 'spin' % F12
    % no data
  case 'fit' % F13
    fread(h,90120,'*char');
  case {'fRtime','tRtime'} % F14
    fread(h,4264,'*char');
  case {'add','sub','mul','div'} % F15
    % no data
  case 'phase' % F16
    jeol_calc.theta = fread(h,1,'float32');
    jeol_calc.radian = fread(h,1,'float32');
    jeol_calc.comment = readcstring(h,256);
    fread(h,544,'*char');
  case 'mT2MHz' % F17
    fread(h,4,'*char');
  case 'hw' % F18
    fread(h,40,'*char');
  case 'mem' % F19
    fread(h,2708,'*char');
  case 'combine' % F20
    % missing in JEOL documentation
  otherwise
    % For any other processType, there is no CALCULATOR block
end

% [G] DATA block
%====================================================================
% The official documentation is incomplete here. It only mentions auxNum=1
% and auxNum=12, whereas most JEOL files with 1D data have auxNum=0
auxNum = jeol_header.auxNum; % format identifier, can be 0, 1 or 12
dataSet = jeol_header.dataSet; % number of datasets
dataKind = jeol_header.dataKind; % 'real' or 'complex' or 'imagi'

jeol_data = struct;
Abscissa = [];
Data = [];

if (dataSet==1)
  
  % Type-1 data, 1D, real and complex; (G1) and (G2)
  %---------------------------------------------------------
  % SEQ_VALUE and SEQ_ITEM are documented, but appear to not be present in the
  % example files provided by JEOL
  %jeol_data.SEQ_VALUE = fread(h,1,'float32');
  %jeol_data.SEQ_ITEM = fread(h,8,'*char').';
  jeol_data.CH1_DATA = fread(h,jeol_header.dataLength,'float32');
  if strcmp(dataKind,'complex')
    jeol_data.CH2_DATA = fread(h,jeol_header.dataLength,'float32');
  end
  Abscissa = jeol_header.xOffset+linspace(0,jeol_header.xRange,jeol_header.dataLength);
  
elseif (auxNum==1) && strcmp(dataKind,'complex')
  
  % Type-1 data, 2D , complex; (G4)
  %---------------------------------------------------------
  for k=1:jeol_header.dataSet
    jeol_data.SEQ_VALUE(k,:) = fread(h,1,'float32');
    jeol_data.SEQ_ITEM(k,:) = fread(h,8,'*char').';
    jeol_data.CH1_DATA(k,:) = fread(h,jeol_header.dataLength,'float32');
    jeol_data.CH2_DATA(k,:) = fread(h,jeol_header.dataLength,'float32');
  end
  Abscissa = jeol_header.xOffset+linspace(0,jeol_header.xRange,jeol_header.dataLength);
  
elseif (auxNum==12) && strcmp(dataKind,'complex')
  
  % Type-2 data, 1D and 2D, complex; (G3) and (G5)
  %---------------------------------------------------------
  for k=1:jeol_header.dataSet
    jeol_data.SEQ_VALUE(k,:) = fread(h,12,'float32');
    jeol_data.SEQ_ITEM(k,:) = fread(h,96,'*char').';
    jeol_data.CH1_DATA(k,:) = fread(h,jeol_header.dataLength,'float32');
    jeol_data.CH2_DATA(k,:) = fread(h,jeol_header.dataLength,'float32');
  end
  Abscissa = jeol_header.xOffset+linspace(0,jeol_header.xRange,jeol_header.dataLength);
end

if strcmp(jeol_header.dataKind,'real')
  Data = jeol_data.CH1_DATA;
  %plot(Abscissa,Data);
elseif strcmp(jeol_header.dataKind,'complex')
  Data = complex(jeol_data.CH1_DATA,jeol_data.CH2_DATA);
  %plot(Abscissa,real(Data),Abscissa,imag(Data));
end

Parameters.Header = jeol_header;
Parameters.General = jeol_general;
Parameters.Spectrometer = jeol_spectrometer;
Parameters.Generator = jeol_generator;
Parameters.ProcessPara = jeol_processpara;
Parameters.Calculator = jeol_calc;
Parameters.Data = jeol_data;

fclose(h);


%-----------------------------------------------------------------
function ch = readcstring(ID,numchars)
% Read C string characters from file ID and convert to Matlab string
ch = fread(ID,numchars,'*char').';
ch = ch(1:find(ch==char(0))-1);
if isempty(ch), ch = ''; end
