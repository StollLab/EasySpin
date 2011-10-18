% nmass   neutron mass 
%
%   mn = nmass
%   [mn,sigma] = nmass
%
%   Returns the mass of the neutron
%   in kg. sigma is the standard
%   uncertainty (2006 CODATA).

function [mn,sigma] = nmass

% 2006 CODATA value

% Value 1.674 927 211 x 10-27 kg
% Standard uncertainty 0.000 000 084 x 10-27 kg
% Relative standard uncertainty 5.0 x 10-8
% Concise form 1.674 927 211(84) x 10-27 kg

mn =    1.674927211e-27;
sigma = 0.000000084e-27;
