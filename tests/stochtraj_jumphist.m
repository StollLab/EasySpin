function [err,data] = test(opt,olddata)
% Check that using stochtraj_jump generates proper state probability
% distribution

Par.nTraj = 500;
Sys.TransRates = 1e9*[-0.5,  0.5;
                       0.5, -0.5];
Sys.Orientations = [0,  0, 0;
                    0, pi, 0];

tau = -1/(2*Sys.TransRates(1,1));  % mean residence time
            
Par.dt = tau/5;
Par.nSteps = ceil(200*tau/Par.dt);

[~,~,~,stateTraj] = stochtraj_jump(Sys,Par);

[n,~,~] = histcounts(stateTraj,'Normalization','Probability');

if abs(n(1)-n(2))>0.05
  err = 1;
  histogram(stateTraj)
else  
  err = 0;
end

data = [];

end
