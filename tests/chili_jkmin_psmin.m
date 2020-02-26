function ok = test(opt)

%===============================================================================
% Test whether fast and general code are internally and among them consistent
% with basis set truncations using Opt.jKmin and Opt.highField.
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
highField = [true false];

% Simulate spectra using all combinations of jKmin and highField
%-------------------------------------------------------------------------------
ok = false;
idx = 0;
for p = 1:numel(highField)
  for j = 1:numel(jKmin)
    Opt.jKmin = jKmin(j);
    Opt.highField = highField(p);
    
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
      title(sprintf('jKmin = %d, highField = %d',Opt.jKmin,Opt.highField));
    end

    % Determine whether spectra are identical
    ok(idx) = areequal(y_fast,y_general,1e-4,'rel');

  end
end
