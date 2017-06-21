clear Params Events system Vary
sqn = 1/2;
% Det = {spops(sqn,'x') spops(sqn,'z') spops(sqn,'p')/2};
Det = {spops(1/2,'z')};
% 
system.sqn = sqn;
system.interactions = {1,0,'z','e',1.5};
system.T1 = 1;
system.T2 = 5;
options.silent = 0;
options.relaxation = 1;
[system,state] = setup(system,options);
Ham = system.ham*1000/2/pi;
Sigma = state;
Relaxation.equilibriumState = system.eq;
Relaxation.Gamma = system.gamma;

Params.tp = 0.100; % us
% Params.Type = 'rectangular/linear';
Params.Frequency = [-100 100] + 1500; % MHz
Params.Flip = pi;
Params.TimeStep = 0.0001; % us
Params.Type = 'quartersin/linear';
Params.trise = 0.015; % us


for i = 1 : 4
    if mod(i,2)
        Events{i}.type = 'pulse';
        [t,IQ,modulation] = pulse(Params);
        IQ = IQ;
        Events{i}.PhaseCycle = [0 1];
        IQ(2,:) = IQ;
        Events{i}.PhaseCycle = [0 1; pi 1];
        Events{i}.t = t;
        Events{i}.IQ = IQ;
        Events{i}.xOp = spops(sqn,'x');
        Events{i}.ComplexExcitation = 0;
%         
    else
        Events{i}.type = 'free evolution';
        Events{i}.t = 0.2;
%         Events{i}.t = 0.2;
    end
    Events{i}.storeDensityMatrix = 1;
    Events{i}.Detection = 1;
    Events{i}.Propagation = [];
    Events{i}.Relaxation = false;
%     Events{i}.Relaxation = true;
end
% 
Vary.Table = [1, 1; 1, 2; 2 1; 2 2];
Vary.Events = {[1] 2};
Vary.Dimension{1}.IQs{1} = {IQ IQ/2};
Vary.Dimension{1}.IQs{3} = {IQ/4 IQ};
Vary.Dimension{1}.ts{1} = {t t};
Vary.Dimension{1}.ts{3} = {t t};

Vary.Dimension{2}.ts{2} = {0:Params.TimeStep:0.2; 0:Params.TimeStep:0.2};



% Events{1}.Relaxation = true;
% Events{3}.Relaxation = true;
Events{1}.Detection = 1;
Events{3}.Detection = 1;

% tic
[t, signal,state,sigmas,Eventsnew]=evolve2(Sigma, Ham, Det, Events, Relaxation, Vary);
% toc
disp(state)

try
a = squeeze(signal(:,:,1));
b = squeeze(signal(:,:,2));
c = squeeze(signal(:,:,3));
d = squeeze(signal(:,:,4));
figure(1); clf
plot(t,real(a))
hold on
plot(t,real(b))
plot(t,real(c))
plot(t,real(d))

catch
  figure(3);clf
  hold on
  for i = 1: length(signal)
    plot(t{i},real(signal{i}))
    
  end
end

% 
% Events{1}.Detection = 0;
% Events{3}.Detection = 0;
% Events{2}.Detection = 0;
% Events{4}.Detection = 0;
%         
% tic
% [t, signal,state,sigma,Eventnew]=evolve2(Sigma,Ham, Det, Events, Relaxation);
% toc
% % profile on
% tic
% [t, signal,state,sigmas]=evolve2(Sigma,Ham, Det, Eventsnew, Relaxation, Vary);
% % profile off
% toc
% % profile report
% disp(state)
% % figure(2); clf
% % plot(t,real(signal))
% % 
% Events{1}.Detection = 0;
% Events{3}.Detection = 0;
% Events{2}.Detection = 0;
% Events{4}.Detection = 0;
%         
% tic
% [t, signal,state,sigma,Eventsnew]=evolve2(Sigma,Ham, Det, Events, Relaxation,Vary);
% toc
% disp(state)
% 
% tic
% [t, signal,state,sigma,Eventsnew]=evolve2(Sigma,Ham, Det, Eventsnew, Relaxation,Vary);
% toc
% disp(state)