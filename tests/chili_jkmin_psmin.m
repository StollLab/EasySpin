function [err,data] = test(opt,olddata)

%===============================================================================
% Test whether Freed and general code are internally and among them consistent
% with basis set truncations using Opt.jKmin and Opt.pSmin.
%===============================================================================

% Set up experiments, system, options
%-------------------------------------------------------------------------------
Exp.mwFreq = 9.5;
Exp.Range = [326 344];

Sys.g = [2.05 2.03 2.00];
Sys.tcorr = [1 2 3]*1e-9; % rhombic diffusion
Sys.lambda = [1 1 1 1]; % complex potential

Opt.pqOrder = true;
Opt.LLKM = [6 3 2 2];

jKmin = [-1 1];
pSmin = [-1 0 1];

% Simulate spectra using all combinations of jKmin and pSmin
%-------------------------------------------------------------------------------
idx = 0;
for j = 1:numel(jKmin)
  for p = 1:numel(pSmin)
    Opt.jKmin = jKmin(j);
    Opt.pSmin = pSmin(p);
    
    idx = idx+1;
    Opt.LiouvMethod = 'Freed';
    [x,y1(:,idx)] = chili(Sys,Exp,Opt);
    
    Opt.LiouvMethod = 'general';
    [x,y2(:,idx)] = chili(Sys,Exp,Opt);
  end
end

% Determine whether spectra are identical
%-------------------------------------------------------------------------------
err = false;
for k = 1:size(y1,2)
  err = err || ~areequal(y1(:,k),y2(:,k),1e-2*max(y1(:,k)));
end

% Plotting
%-------------------------------------------------------------------------------
if opt.Display
  cla
  sh = 5;
  for k = 1:6
    plot(x,y1(:,k)+sh*(k-1),x,y2(:,k)+sh*(k-1));
  end
end


data = [];
