%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_AdaniDAT(FileName)
%-------------------------------------------------------------------------------
%   Text-based file format of Adani spectrometers
%-------------------------------------------------------------------------------
fh = fopen(FileName);
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
allLines = allLines{1};
fclose(fh);
nLines = numel(allLines);

Line1 = '======================== Parameters: ========================';
if ~strncmp(allLines{1},Line1,length(Line1))
  error('The file %s is not an Adani spectrometer file',FileName);
end
Line2 = '=============================================================';
for idx = 1:nLines
  found = strncmp(allLines{idx},Line2,length(Line2));
  if found; break; end
end
if (~found)
  error('Could not find start of data in file %s',FileName);
end
nPoints = nLines - idx;
data = zeros(nPoints,3);
for p=1:nPoints
  L_ = allLines{idx+p};
  L_(L_==',') = '.';
  data(p,:) = sscanf(L_,'%f',3);
end
Abscissa = data(:,2);
Data = data(:,3);
Parameters = [];

return
%-------------------------------------------------------------------------------
