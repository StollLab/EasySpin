function [Data,Abscissa,Parameters] = eprload_CIQTEK(fileName)
%--------------------------------------------------------------------------
% JSON format of CIQTEK spectrometers
%--------------------------------------------------------------------------
% (implementation based on example files from the EasySpin user forum)
% Uses Matlab-internal JSON parser  matlab.internal.webservices.fromJSON
% This does not work prior to R2016b.

% Matlab version check
if verLessThan('matlab','9.1.0')
  error('Reading CIQTEK files requires Matlab R2016b (9.1.0) or later.');
end

% Read JSON string from file
fid = fopen(fileName);
if fid<0
  error('Could not open CIQTEK file %s.',fileName);
end
jsonstring = fread(fid, "*char").';
fclose(fid);

% Parse JSON string into MATLAB structure
try
  data = jsondecode(jsonstring);
catch
  error('Could not parse CIQTEK file %s.',fileName);
end

% Extract field axis and spectrum
lineData = data.dataStore.lineDataList;
nSpectra = numel(lineData);
for iSpectrum = 1:nSpectra
  B(:,iSpectrum) = lineData(iSpectrum).ReData(:,1);
  spcRe(:,iSpectrum) = lineData(iSpectrum).ReData(:,2);
  spcIm(:,iSpectrum) = lineData(iSpectrum).ImData(:,2);
end

% Output
Data = spcRe;
Abscissa = B;
Parameters = data;

end
