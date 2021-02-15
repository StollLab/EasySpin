function [Data,Abscissa,Parameters] = eprload_AdaniJSON(FileName)
%--------------------------------------------------------------------------
% JSON format of Adani spectrometers SPINSCAN etc. (e-Spinoza software)
%--------------------------------------------------------------------------
% (implementation based on official documentation from Adani)
% Uses Matlab-internal JSON parser  matlab.internal.webservices.fromJSON
% This does not work prior to R2016b.

% Matlab version check
if verLessThan('matlab','9.1.0')
  error('Reading Adani e-Spinoza JSON files requires Matlab R2016b (9.1.0) or later.');
end

% Read JSON string from file
fid = fopen(FileName);
if fid<0
  error('Could not open %s.',FileName);
end
jsonstring = textscan(fid,'%s');
jsonstring = jsonstring{1};
fclose(fid);

% Parse JSON string
try
  data = matlab.internal.webservices.fromJSON(jsonstring);
  data = data{1};
catch
  error('Could not parse JSON-format data from %s.',FileName);
end

% Get basic information about data
nSpectra = numel(data.Values);
nPoints = numel(data.Values(1).Values);
nPhases = numel(data.Values(1).Values(1).Points);

% Read 1D or 2D data
if isfield(data,'Phase')
  phase0 = data.Phase;
else
  phase0 = 0;
end
s = sin(2*pi*(0:nPhases-1).'/nPhases + phase0);
spc = zeros(nPoints,nSpectra);
for iSpectrum = 1:nSpectra
  acqdata = data.Values(iSpectrum).Values;
  for iPoint = 1:nPoints
    spc(iPoint,iSpectrum) = sum(acqdata(iPoint).Points.*s);
  end
end

% Construct field axis
d_ = data.ExperimentOptions;
CenterField = d_.CommonOptions.CenterMagneticField;
SweepWidth = d_.CommonOptions.SweepWidth;
B = linspace(-1,1,nPoints)*SweepWidth/2 + CenterField;

% Construct second axis
if nSpectra>1
  axis2 = linspace(d_.InitialValue2D,d_.FinalValue2D,nSpectra);
end

% Output
Data = spc;
if nSpectra==1
  Abscissa = B;
else
  Abscissa = {B,axis2};
end
Parameters = data;

return
