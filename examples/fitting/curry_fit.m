% Combined fit of magnetic susceptibility and magnetic moment

clear all
close all

% generate 'data'
DataSys.S = 3/2;
DataSys.D = -20*clight*1e-4;
DataSys.g = 2.3;

% Temperatures and field for magnetic susceptibilty:
Exp.chiTemperature =[2:0.5:20 25:5:300];
Exp.chiField =[100,500];

% Temperatures and field for magnetic moment
Exp.mTemperature =[2,5,10];
Exp.mField = 0:100:7000;

OptData.Output = 'ChiTCGs MvsB'; 

% generate 'data' and add a bit of noise
[ChiTD,MvsBD] = curry(DataSys,Exp,OptData);
ChiTData = addnoise(ChiTD,30);
MvsBData = addnoise(MvsBD,30);

figure(1)
subplot(2,1,1)
plot(Exp.chiTemperature,ChiTData,'+')
title('\chi T for S = 3 /2 with D = -20 cm^{-1} and g = 2.3');
subplot(2,1,2)
plot(Exp.mField,MvsBData,'+')
title('M(B)for S = 3 /2 with D = -20 cm^{-1} and g = 2.3');
ylim([0 2.2])

%% now we have 'data', let's do a least square fit
% make a first guess for the fit system
FitSys.S = 3/2;
FitSys.D = 20*clight*1e-4;
FitSys.g = 2.5;

%what should be fitted
Vary.D = 40*clight*1e-4;
Vary.g = 0.5;

% no scaling! Assuming correct diamagnetic corrections and 
% sample weight measure the absolute value contain valuable information!
FitOpt.Scaling = 'none';

% fit chiT and mag seperately
ChiOpt.Output = 'ChiTCGs OneColoumn'; 
bestsysChi = esfit('curry',ChiTData(:),FitSys,Vary,Exp,ChiOpt,FitOpt);
ChiOptData.Output = 'ChiTCGs';
ChiTSim = curry(bestsysChi,Exp,ChiOptData);

MOpt.Output = 'MvsB OneColoumn'; 
bestsysM = esfit('curry',ChiTData(:),FitSys,Vary,Exp,ChiOpt,FitOpt);
MOptData.Output = 'MvsB';
MvsBSim = curry(bestsysM,Exp,MOptData);
figure(2)
subplot(2,1,1)
plot(Exp.chiTemperature,ChiTData,'+', Exp.chiTemperature,ChiTSim)
title(['\chi T for S = 3 /2 with D = -20 cm^{-1} and g = 2.3 fitted with D = ', ...
  num2str(bestsysChi.D/clight*1e4,'%2.1f'), ' cm^{-1} and g = ', num2str(bestsysChi.g,'%1.2f')]);
subplot(2,1,2)
plot(Exp.mField,MvsBData,'+',Exp.mField,MvsBSim)
title(['M(B) for S = 3 /2 with D = -20 cm^{-1} and g = 2.3 fitted with D = ', ...
  num2str(bestsysM.D/clight*1e4,'%2.1f'), ' cm^{-1} and g = ', num2str(bestsysM.g,'%1.2f')]);
ylim([0 2.2])

%not too bad, but we can do better!

%%

% combine all data in 1-row array
fitdata = [ChiTData(:); MvsBData(:)];

% set output to ChiT and MvsB in 1-coloumn array
Opt.Output = 'ChiTCGs MvsB OneColoumn'; 

% no scaling! Assuming correct diamagnetic corrections and 
% sample weight measure the absolute value contain valuable information!
FitOpt.Scaling = 'none';

% Let's go and hope for the best!
[bestsys,bestspc] = esfit('curry',fitdata,FitSys,Vary,Exp,Opt,FitOpt);

[ChiTSim,MvsBSim] = curry(bestsys,Exp,OptData);
figure(3)
subplot(2,1,1)
plot(Exp.chiTemperature,ChiTData,'+', Exp.chiTemperature,ChiTSim)
title(['\chi T for S = 3 /2 with D = -20 cm^{-1} and g = 2.3 fitted with D = ', ...
  num2str(bestsys.D/clight*1e4,'%2.1f'), ' cm^{-1} and g = ', num2str(bestsys.g,'%1.2f')]);
subplot(2,1,2)
plot(Exp.mField,MvsBData,'+',Exp.mField,MvsBSim)
title(['M(B) for S = 3 /2 with D = -20 cm^{-1} and g = 2.3 fitted with D = ', ...
  num2str(bestsys.D/clight*1e4,'%2.1f'), ' cm^{-1} and g = ', num2str(bestsys.g,'%1.2f')]);
ylim([0 2.2])

% the fit look okay, doesn't it? let's try negative D

%%
FitSys2.S = 3/2;
FitSys2.D = -10*clight*1e-4;
FitSys2.g = 2.0;

% Let's go and hope for the better!
[bestsys2,bestspc2] = esfit('curry',fitdata,FitSys2,Vary,Exp,Opt,FitOpt);

[ChiTSim2,MvsBSim2] = curry(bestsys2,Exp,OptData);
figure(4)
subplot(2,1,1)
plot(Exp.chiTemperature,ChiTData,'+', Exp.chiTemperature,ChiTSim2)
title(['\chi T for S = 3 /2 with D = -20 cm^{-1} and g = 2.3 fitted with D = ', ...
  num2str(bestsys2.D/clight*1e4,'%2.1f'), ' cm^{-1} and g = ', num2str(bestsys2.g,'%1.2f')]);
subplot(2,1,2)
plot(Exp.mField,MvsBData,'+',Exp.mField,MvsBSim2)
title(['M(B) for S = 3 /2 with D = -20 cm^{-1} and g = 2.3 fitted with D = ', ...
  num2str(bestsys2.D/clight*1e4,'%2.1f'), ' cm^{-1} and g = ', num2str(bestsys2.g,'%1.2f')]);
ylim([0 2.2]);

% now we are close, right?
% BTW, if you are patient enough write down all the best fit values and
% rerun the script. Interesting, right? 


