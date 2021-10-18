%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_MagnettechBinary(FileName)
%-------------------------------------------------------------------------------
%   Binary file format of older Magnettech spectrometers (MS400 and prior)
%-------------------------------------------------------------------------------

hMagnettechFile = fopen(FileName,'r','ieee-le');
if hMagnettechFile<0
  error('Could not open Magnettech spectrometer file %s.',FileName);
end

% Check file size to determine whether this is a Magnettech data file
% - File overhead is 64 bytes, located at end of file
% - Each data point is 2 bytes
% - MiniScope MS400 benchtop spectrometer data have 4096 points
% - Data from other spectrometers (e.g MT500 L-band) can have 512,
%   1024, 2048, or 4096 points.
nPoints = [512 1024 2048 4096];
fileinfo = dir(FileName);
fileSize = fileinfo.bytes;
idx = find(fileSize==2*nPoints+64);
if isempty(idx)
  error('This file does not have the correct file size for a Magnettech SPE file.');
else
  nPoints = nPoints(idx);
end

% Read spectral data
Data = fread(hMagnettechFile,nPoints,'int16');

% Determine format version and other flags
fseek(hMagnettechFile,fileSize-4,'bof');
FileFlags = fread(hMagnettechFile,1,'uint8');
mwFreqAvailable = bitand(FileFlags,1)~=0;
oldSpeFormat = bitand(FileFlags,2)==0;
temperatureAvailable = bitand(FileFlags,4)~=0;
if oldSpeFormat
  fourbytesingle = @(x)x(1) + x(2)/100;
else
  fourbytesingle = @(x)x(1) + x(2)/1000;
end
readfbs = @()fourbytesingle(fread(hMagnettechFile,2,'int16'));

%if oldSpeFormat
%  Data = Data - 16384;
%end

% Read parameters
fseek(hMagnettechFile,2*nPoints,'bof');
Parameters.B0_Field = readfbs()/10; % G -> mT
Parameters.B0_Scan = readfbs()/10; % G -> mT
Parameters.Modulation = readfbs()/10000; % mT
Parameters.MW_Attenuation = readfbs(); % dB
Parameters.ScanTime = readfbs(); % s
GainMantissa = readfbs();
GainExponent = readfbs();
Parameters.Gain = GainMantissa*10^round(GainExponent);
Parameters.Number = readfbs();
reserve = readfbs();
Parameters.Time_const = readfbs(); % s
reserve = readfbs();
reserve = readfbs();
Parameters.NumberSamples = readfbs();
if temperatureAvailable
  Parameters.Temperature = fread(hMagnettechFile,1,'int32'); % degree C
else
  reserve = fread(hMagnettechFile,1,'int32');
  Parameters.Temperature = [];
end
reserve = readfbs();
Parameters.FileFlags = fread(hMagnettechFile,1,'uint8');
if mwFreqAvailable
  mwf = fread(hMagnettechFile,3,'uint8');
  Parameters.mwFreq = mwf(3) + 256*mwf(2) + 256^2*mwf(1);
  Parameters.mwFreq = Parameters.mwFreq/1e6; % kHz->GHz
else
  Parameters.mwFreq = [];
end
fclose(hMagnettechFile);

Abscissa = Parameters.B0_Field + linspace(-1/2,1/2,numel(Data))*Parameters.B0_Scan;
Abscissa = Abscissa(:);

return
%-------------------------------------------------------------------------------
