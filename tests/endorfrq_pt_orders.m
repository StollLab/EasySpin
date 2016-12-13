function [err,data] = test(opt,olddata)

%-------------------------------------------------------------
% Make sure 2nd-order PT is always better than 1st-order PT
%-------------------------------------------------------------
clear
Exp.Field = 350;

% Loop over + and - signs for both A and gn
A = [+120 -120];
Nucs = {'1H','15N'};

for iA = 1:2
  for iNucs = 1:2
    Sys.Nucs = Nucs{iNucs};
    Sys.A = A(iA);
    
    freq_m = endorfrq(Sys,Exp);
    freq_m = sort(freq_m);
    
    Opt.PerturbOrder = 1;
    freq_p1 = endorfrq_perturb(Sys,Exp,Opt);
    freq_p1 = sort(freq_p1);
    
    Opt.PerturbOrder = 2;
    freq_p2 = endorfrq_perturb(Sys,Exp,Opt);
    freq_p2 = sort(freq_p2);
    dfreq1 = abs(freq_p1-freq_m);
    dfreq2 = abs(freq_p2-freq_m);
    ok(iA,iNucs) = all(dfreq2<dfreq1*1.01);
  end
end

err = any(~ok(:));

data = [];
