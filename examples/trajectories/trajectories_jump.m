% trajectories of Markovian jumps between discrete states
%==========================================================================
% This script simulates jumps between 2 states using the provided
% transition probability matrix and time step. By calling stochtraj_jump
% without any output arguments, a single trajectory of states will be
% plotted, along with a histogram of state counts.

clear

Sys.TransProb = [ 0.1, 0.9;
                  0.3, 0.7 ];  % transition probability matrix

Par.dt = 1e-9;                 % time step, s

stochtraj_jump(Sys,Par)