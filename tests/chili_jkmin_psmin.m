function [err,data] = test(opt,olddata)

%===============================================================================
% Test whether fast and general code are internally and among them consistent
% with basis set truncations using Opt.jKmin and Opt.pSmin.
%===============================================================================

% Set up experiments, system, options
%-------------------------------------------------------------------------------
Exp.mwFreq = 9.5;
Exp.Range = [326 344];

Sys.g = [2.05 2.03 2.00];
Sys.tcorr = [1 2 3]*1e-9; % rhombic diffusion
Sys.Potential = [2 0 0 1; 2 0 2 1; 4 0 0 1; 4 0 2 1];

Opt.pqOrder = true;
Opt.LLMK = [6 3 2 2];
Opt.nKnots = 5;

jKmin = [-1 1];
pSmin = [-1 1];

% Simulate spectra using all combinations of jKmin and pSmin
%-------------------------------------------------------------------------------
err = false;
idx = 0;
for p = 1:numel(pSmin)
  for j = 1:numel(jKmin)
    Opt.jKmin = jKmin(j);
    Opt.pSmin = pSmin(p);
    
    idx = idx+1;
    Opt.LiouvMethod = 'fast';
    [B,y_fast] = chili(Sys,Exp,Opt);
    
    Opt.LiouvMethod = 'general';
    [B,y_general] = chili(Sys,Exp,Opt);
    
    % Plotting
    if opt.Display
      subplot(2,2,idx);
      plot(B,y_fast,B,y_general);
      axis tight
      legend('fast','general');
      title(sprintf('jKmin = %d, pSmin = %d',Opt.jKmin,Opt.pSmin));
    end

    % Determine whether spectra are identical
    err(idx) = ~areequal(y_fast,y_general,1e-4,'rel');

  end
end
err = any(err);

data = [];
