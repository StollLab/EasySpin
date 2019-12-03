function [err,data] = eprload_esp_xaxis(opt,olddata)

data = [];
err = [];

BaseDir = 'eprfiles/';

%-------------------------------------------------
% Check whether ESP abscissa vectors are correct
%-------------------------------------------------

Files{1} = '00011201.par';
Files{2} = 'esp.par';
Files{3} = 'frem_gly.par';
Files{4} = 'simfonia125.par';
Files{5} = 'strong1esp.par';
Files{6} = 'win.par';
Files{7} = 'inga17MAY3.PAR';
Files{8} = 'inga2ABR01.PAR';
Files{9} = 'inga_mar06_endor.PAR';
Files{10} = 'pandelia_mar06.PAR';
Files{11} = 'PNOX120DB.PAR';
Files{12} = 'PNOX220DB.PAR';
Files{13} = 'loncke_endor.par';

%Abscissa{1} = [200 4200]; %HCF/HSW
Abscissa{1} = [199.512 4199.512]; %GST/GSI
Abscissa{2} = [3463 3513];
%Abscissa{3} = [3443.8 3503.8]; % HCF/HSW
Abscissa{3} = [3449.71 3509.73]; % GST/GSI
%Abscissa{4} = 3357.213731+[-1 1]/2*141.121367; % HCF/HSW
Abscissa{4} = 3286.653048+[0 1]*141.121367; % GST/GSI
Abscissa{5} = [3375 3425];
Abscissa{6} = [3463 3513];
Abscissa{7} = [0 7632];
Abscissa{8} = [0 7632];
Abscissa{9} = [0.8212386 26.0212386];
Abscissa{10} = [3412.4 3462.4];
Abscissa{11} = [3421.1 3471.1];
Abscissa{12} = [3420.65 3470.65];
Abscissa{13} = [1 81];

%Univ. Jena ESP300E, cw EPR
Files{14} = 'jena_mi191000.par';
Files{15} = 'jena_se120802.par';
Files{16} = 'jena_ap151000.par';
Abscissa{14} = 3456.444 + [0  80.13]; % GST/GSI (differs from HCF/HCW)
Abscissa{15} = 3464.134 + [0  70.16]; % GST/GSI (differs from HCF/HCW)
Abscissa{16} = 3444.364 + [0 100.26]; % GST/GSI (differs from HCF/HCW)

xax_err = zeros(1,numel(Files));
for iFile = 1:numel(Files)
  if (opt.Display)
    disp(Files{iFile});
  end
  [xax,data] = eprload([BaseDir Files{iFile}]);
  xax0 = Abscissa{iFile};
  if ~areequal(xax(1),xax0(1),1e-4,'rel') || ~areequal(xax(end),xax0(2),1e-4,'rel')
    disp([' ' Files{iFile}]);
    fprintf('   read:      %0.9f %0.9f\n',xax([1 end]));
    fprintf('   should be: %0.9f %0.9f\n',xax0);
    xax_err(iFile) = 1;
  end
end

err = xax_err;

data = [];
