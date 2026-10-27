% p_resfieldspan  Estimate lowest and highest resonance field of an S=1/2 system
%
%   [minB,maxB] = p_resfieldspan(Sys,mwFreq)
%
%   Input:
%     Sys     spin system, processed by validatespinsys
%     mwFreq  microwave frequency, in GHz
%
%   Output:
%     minB    lowest resonance field, in mT
%     maxB    highest resonance field, in mT
%
%   The estimate uses the principal g values and, for each nucleus, the
%   largest principal hyperfine coupling, to first order. It is only
%   available for a single electron spin S=1/2 without orbital angular
%   momentum. Otherwise, an error is issued.

function [minB,maxB] = p_resfieldspan(Sys,mwFreq)

if ~isscalar(Sys.S) || Sys.S~=1/2 || any(Sys.L(:))
  error('Cannot automatically determine the field range for this spin system. Provide Exp.Range or Exp.CenterSweep.');
end

% Magnitudes of the principal g values
if Sys.fullg
  g = svd(Sys.g);
else
  g = abs(Sys.g(:));
end

% Largest principal hyperfine coupling of each nucleus
if Sys.nNuclei>0
  I = nucspin(Sys.Nucs);
  if Sys.fullA
    Amax = zeros(1,Sys.nNuclei);
    for iNuc = 1:Sys.nNuclei
      Amax(iNuc) = norm(Sys.A(3*iNuc-2:3*iNuc,:));
    end
  else
    Amax = max(abs(Sys.A),[],2).';
  end
  hf = sum(I.*Amax)*1e6; % MHz -> Hz
else
  hf = 0;
end

minB = planck*(mwFreq*1e9 - hf)/bmagn/max(g)/1e-3; % mT
maxB = planck*(mwFreq*1e9 + hf)/bmagn/min(g)/1e-3; % mT

end
