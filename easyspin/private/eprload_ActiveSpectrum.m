%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_ActiveSpectrum(FileName)
%-------------------------------------------------------------------------------
%   ESR file format of Active Spectrum spectrometers
%-------------------------------------------------------------------------------
fh = fopen(FileName);
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
allLines = allLines{1};
fclose(fh);
nLines = numel(allLines);

% Find start of data
for idx = 1:nLines
  found = strncmp(allLines{idx},'FIELD (G)',9);
  if found; break; end
end
if (~found)
  error('Could not find start of data in file %s',FileName);
end

dataLines = allLines(idx+1:end);
nPoints = numel(dataLines);

for idx = 1:nPoints
  data(idx,:) = sscanf(dataLines{idx},'%f',2);
end
Abscissa = data(:,1);
Data = data(:,2);
Parameters = [];
return
%-------------------------------------------------------------------------------

