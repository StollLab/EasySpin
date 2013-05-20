function Nucs = nucstringparse(NucString,Mode)

% Mode:
% - 's' for single isotope per list entry
% - 'm' for multiple isotopes per list entry

% Splits list of nuclei into a cell array of strings,
% e.g. '63Cu,14N,14N,1H'  -> {'63Cu','14N','14N','1H'}
%      '(63,65)Cu,N,35Cl' -> {'(63,54)Cu','N','35Cl'}

if (nargin==0)
  NucString = '(63,65)Cu,N,14N,13C';
end

if iscell(NucString)
  Nucs = NucString;
  return;
end

if isempty(NucString)
  Nucs = [];
  return;
end

if (nargin==1)
  Mode = 's';
end

if (Mode=='s')
  SingleIsotopologue = 1;
elseif (Mode=='m')
  SingleIsotopologue = 0;
else
  error('Unknown parse mode.');
end

str = NucString;

str(isspace(str)) = []; % remove whitespace
if str(end)==',', str(end) = []; end % remove trailing comma if present
N = length(str);

if (SingleIsotopologue)
  commaidx = (str==',');
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
