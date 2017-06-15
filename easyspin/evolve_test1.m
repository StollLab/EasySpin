clear Params Events system
sqn = 1/2;
Det = {spops(sqn,'x') spops(sqn,'z') spops(sqn,'p')/2};
% Det = {spops(1/2,'z')};
% 
system.sqn = sqn;
system.interactions = {1,0,'z','e',1.5};
system.T1 = 2;
system.T2 = 10;
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
        [t,IQ] = pulse(Params);
        IQ = IQ;
%         IQ(2,:) = IQ;
        Events{i}.t = t;
        Events{i}.IQ = IQ;
        Events{i}.xOp = spops(sqn,'x');
        Events{i}.ComplexExcitation = 0;
        Events{i}.PhaseCycle = [0 1];
%                                 pi -1];
%         
    else
        Events{i}.type = 'free evolution';
        Events{i}.t = 0:Params.TimeStep:0.2;
%         Events{i}.t = 0.2;
    end
    Events{i}.storeDensityMatrix = 1;
    Events{i}.Detection = 1;
    Events{i}.propagators = [];
end

% Events{3}.Detection = 0;
% Events{4}.Detection = 1;

tic
[t, signal,state,sigmas,Eventnew]=evolve2(Sigma,Ham, Det, Events, Relaxation);
toc
disp(state)
figure(1); clf
plot(t,real(signal))

Events{3}.ComplexExcitation = 1;
Events{1}.ComplexExcitation = 1;
        
tic
% profile on
[t, signal,state,sigmas,Eventnew]=evolve2(Sigma,Ham, Det, Events, Relaxation);
% profile off
toc
% profile report
disp(state)
figure(2); clf
plot(t,real(signal))