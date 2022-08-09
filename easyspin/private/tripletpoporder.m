function pops = tripletpoporder(varargin)
% Take triplet zero-field populations (px, py, pz) and output energy-ordered populations
% (from lowest to highest energy) given the signs of D and E
%
% pops = tripletpoporder(Dprincipal,xyzpops)
% pops = tripletpoporder(D,E,xyzpops)
%
% Input: D          - zero-field splitting D value (with sign!)
%        E          - zero-field splitting E value
%        or
%        Dprincipal - principal values for zero-field splitting interaction
%        and
%        xyzpops    - zero-field triplet sublevel populations, [px py pz]
%
% Output: pops      - zero-field triplet sublevel populations ordered from 
%                     lowest to highest energy level
%

% Input parsing and definition of energies of the triplet state sublevels
% at zero-field
if nargin==2
  Dprincipal = varargin{1};
  xyzpops = varargin{2};
  % Energies
  En = - Dprincipal;
elseif nargin==3
  D = varargin{1};
  E = varargin{2};
  xyzpops = varargin{3};
  % Energies
  En(1) = (1/3)*D - E; % X
  En(2) = (1/3)*D + E; % Y
  En(3) = -(2/3)*D;    % Z
end

% Get indices for ordering from lowest to highest energy level
[~,ind] = sort(En,'ascend');

% Energy-ordered populations
pops = xyzpops(ind);