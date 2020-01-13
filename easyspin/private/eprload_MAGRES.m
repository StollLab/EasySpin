%-------------------------------------------------------------------------------
function [Data, Abscissa, Parameters] = eprload_MAGRES(FileName)
%-------------------------------------------------------------------------------
% PLT file processing
%   MAGRES  Nijmegen EPR/ENDOR simulation program
%-------------------------------------------------------------------------------

[Line,found] = findtagsMAGRES(FileName,{'DATA'});
if found(1), nx = str2double(Line{1}); else, nx=0; end
if ~nx
  error('Unable to determine number of x points in PLT file.');
end

fid = fopen(FileName,'r');
if (fid<0), error(['Could not open ' FileName]); end

for k=1:3, fgetl(fid); end

% read data
ny = 1;
[Data,N] = fscanf(fid,'%f',[nx,ny]);
if N<nx*ny
  warning('Could not read entire data set from PLT file.');
end

% close file
St = fclose(fid);
if St<0, error('Unable to close PLT file.'); end

Abscissa = [];
Parameters = [];
return
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
function [out,found] = findtagsMAGRES(FileName,TagList)

% open file
fid = fopen(FileName,'r');
if fid<0, error(['Could not open ' FileName]); end

found = zeros(1,length(TagList));
out = cell(1,length(TagList));
while ~feof(fid)
  Line = fgetl(fid);
  whitespace = find(isspace(Line)); % space or tab
  if ~isempty(whitespace)
    endTag = whitespace(1)-1;
    if endTag>0
      I = strcmp(Line(1:endTag),TagList);
      if ~isempty(I)
        out{I} = fliplr(deblank(Line(end:-1:endTag+1)));
        found(I) = 1;
      end
    end
  end
end

% close file
St = fclose(fid);
if St<0, error('Unable to close data file.'); end

return
%-------------------------------------------------------------------------------
