% nuclist2string  Generate isotope list string
%
%   str = nuclist2string(Nucs)
%
%   Generates a string containing a commma-
%   separated list of isotope symbols given
%   in the cell array Nucs.

function NucString = nuclist2string(Nucs)

% Generate list of nuclei from a cell array of strings,
% e.g. {'63Cu','14N','14N','1H'} -> '63Cu,14N,14N,1H'

if ~iscell(Nucs)
  error('Nuclear isotope list must be a cell array!');
end

if isempty(Nucs)
  NucString = '';
  return;
end

% Append comma to each isotope symbol
Nucs = strcat(Nucs,',');

% Concatenate
NucString = strcat(Nucs{:});

% Remove trailing comma
NucString(end) = [];

return
