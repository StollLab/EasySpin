function [FitOpt,info] = esfit_algdefaults(algorithm)

% Option parameters common to all algorithms
info.TolFun = 'termination tolerance for error function value';
FitOpt.TolFun = 1e-5;

info.maxTime = 'time, in minutes, after which esfit will terminate even if the fitting has not yet converged';
FitOpt.maxTime = inf;
FitOpt.IterFcn = [];
FitOpt.IterationPrintFunction = @(str)str;
FitOpt.InfoPrintFunction = @(str)str;

% Algorithm-specific parameters
switch algorithm
  case 'Nelder-Mead simplex'
    
    info.delta = 'edge length for initial simplex';
    FitOpt.delta = 0.1;

    info.TolEdgeLength = 'termination tolerance for the length of the parameter step';
    FitOpt.TolEdgeLength = 1e-4;

    info.SimplexPars = ['[rho, chi, psi, sigma]' ...
                        '(reflection coefficient, expansion coefficient,' ...
                        'contraction coefficient, reduction coefficient)'];
    FitOpt.SimplexPars = []; % [1, 2, 0.5, 0.5]; if nParams>2, FitOpt.SimplexPars = [1, 1+2/nParams, 0.75-1/(2*nParams), 1-1/nParams]; end;

    info.ScaleParams = 'rescale fitting parameters to within [-1,1] interval';
    FitOpt.ScaleParams = false;
    
    FitOpt.RandomStart = 0;

  case 'Levenberg-Marquardt'
    
    info.lambda = 'starting value of Marquardt parameter';
    FitOpt.lambda = 1e-3;
    
    info.delta = 'relative step for difference approximation';
    FitOpt.delta = 1e-3;

    info.TolStep = 'termation tolerance for parameter step (small step stops)';
    FitOpt.TolStep = 1e-6;
    
    info.Gradient = 'termation tolerance for gradient (small gradient stops)';
    FitOpt.Gradient = FitOpt.TolFun;
    
    info.ScaleParams = 'rescale fitting parameters to within [-1,1] interval';
    FitOpt.ScaleParams = false;

    FitOpt.RandomStart = 0;

  case 'Monte Carlo'

    info.nTrials = 'number of random trial simulations before termination';
    FitOpt.nTrials = 20000;

  case 'genetic algorithm'

    info.PopulationSize = 'size of the population, i.e. number of parameter sets and simulations in one generation';
    FitOpt.PopulationSize = 20;

    info.EliteCount = 'number of populations (ordered in terms of decreasing score) carried over to the next generation without recombination and mutation'; 
    FitOpt.EliteCount = []; %'max(2,ceil(0.1*FitOpt.PopulationSize))';

    info.maxGenerations = 'maximum number of generations the algorithm should run';
    FitOpt.maxGenerations = 10000;

  case 'grid search'

    info.GridSize = 'number or an array that specifies how many grid points there should be for each parameter';
    FitOpt.GridSize = 7;

    info.randGrid = 'randomize order of gridpoints';
    FitOpt.randGrid = true;

    FitOpt.maxGridPoints = 1e5;

  case 'particle swarm'

    info.nParticles = 'number of particles in the particle swarm'; 
    FitOpt.nParticles = []; % '20 + nParams*10;';

    info.SwarmParams = ['[k, w, c1, c2]' ...
                        '(velocity limit, inertial coefficient,' ...
                        'cognitive coefficient, social coefficient)'];
    FitOpt.SwarmParams = [0.2 0.5 2 1];

    info.TolStallIter = 'maximum number of consecutive iterations over which the best function value doesn''t change';
    FitOpt.TolStallIter = 6;

end



