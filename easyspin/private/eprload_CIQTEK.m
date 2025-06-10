function [data,abscissa,parameters] = eprload_CIQTEK(fileName)
%--------------------------------------------------------------------------
% JSON format of CIQTEK spectrometers
%--------------------------------------------------------------------------
% (implementation based on example files from the EasySpin user forum)
% This does not work prior to R2016b because of jsondecode().

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
  filecontents = jsondecode(jsonstring);
catch
  error('Could not parse CIQTEK file %s.',fileName);
end

% Extract abscissa and in-phase/quadrature data
lineDataList = filecontents.dataStore.lineDataList;
nTraces = numel(lineDataList);
x(:) = lineDataList(1).ReData(:,1);
for iTrace = nTraces:-1:1
  dataRe(:,iTrace) = lineDataList(iTrace).ReData(:,2);
  dataIm(:,iTrace) = lineDataList(iTrace).ImData(:,2);
end

% Set up 2D abscissae
% For 2D data, lineDataList().params contains a field
%   delay:  2D field/delay sweep
%   power:  2D field/power sweep
%   time2:  pulse experiments
% If no specific field is present, then the second axis is simply indexed
if nTraces>1
  lineParams = lineDataList(1).params;
  y = 1:nTraces;
  fieldNames = {'delay','power','time2'};
  for f = 1:numel(fieldNames)
    if ~isfield(lineParams,fieldNames{f}), continue; end
    for iTrace = 1:nTraces
      y(iTrace) = sscanf(lineDataList(iTrace).params.(fieldNames{f}),'%f');
    end
    break
  end
end

% Output
if nTraces>1
  abscissa = {x,y};
else
  abscissa = x;
end
data = complex(dataRe,dataIm);
parameters = filecontents;

end
