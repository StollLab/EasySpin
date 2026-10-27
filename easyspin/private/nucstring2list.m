% Converts a nuclear spin string to a cell array

% Input
%  NucString: string or character array containing comma-separated list
%             of nuclear spins
%    Syntax: 'X,Y,Z' where X, Y, Z can be either
%    - single isotopes: 13C, 14N, 1H, 63Cu
%    - natural mixtures: C, N, H, Cu
%    - custom isotope mixtures: (63,65)Cu, (12,13)C
%  ParseMode: parse mode
%    - 's' for single isotope per list entry
%    - 'm' for multiple isotopes per list entry
%
% Output: cell array of nuclear spin list, size (1,nNuclei); empty entries
% are dropped, and an empty or all-whitespace input gives an empty list
%
% Splits list of nuclei into a cell array of strings,
% e.g. '63Cu,14N,14N,1H'  -> {'63Cu','14N','14N','1H'}
%      '(63,65)Cu,N,35Cl' -> {'(63,65)Cu','N','35Cl'}

function Nucs = nucstring2list(NucString,ParseMode)

if iscell(NucString)
  Nucs = NucString;
  return
end

if nargin<2
  ParseMode = 's';
end

if ParseMode=='s'
  SingleIsotopologueMode = true;
elseif ParseMode=='m'
  SingleIsotopologueMode = false;
else
  error('Unknown parse mode - should be ''s'' or ''m''.');
end

str = NucString;

if isstring(str)
  if numel(str)~=1
    error('Isotope list must be a single string, not a string array with %d elements.',numel(str));
  end
  str = char(str);
end

if isempty(str)
  Nucs = cell(1,0);
  return
end

if ~ischar(str)
  error('Isotope list must be a string, not numeric.');
end

if size(str,1)~=1
  error('Isotope list must be a single string, not an array with %d rows.',size(str,1));
end

if any(str==';')
  error('Use , and not ; in nuclear spin string. You gave ''%s''.',str);
end

str(isspace(str)) = []; % remove whitespace
N = length(str);

if SingleIsotopologueMode
  commaidx = str==',';
else
  % Determine level (1 inside parentheses, 0 outside)
  level = zeros(1,N);
  level(str=='(') = +1;
  level(str==')') = -1;
  % Find positions of commas outside parentheses
  commaidx = (str==',') & (cumsum(level)==0);
end
comma = [0 find(commaidx) N+1];

nNuclei = numel(comma)-1;

Nucs = cell(1,nNuclei);
for k = 1:nNuclei
  Nucs{k} = str(comma(k)+1:comma(k+1)-1);
end
Nucs = Nucs(1,~cellfun('isempty',Nucs)); % drop empty entries, keep size (1,n)
