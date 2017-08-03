isPulse  = [1 0  1  0  1 0 1 0];
EventLengths = [10 70 10 -30 10 100 10 40];

% close all
% clc, clear

% isPulse  = [1 0  1  0  1 0 1 0 1 0];
% EventLengths = [10 100 10 30 10 40 10 40 10 20];

% [Sequence, T] = reorder_events(EventLengths,isPulse);

% dt1 = -180;
% EventLengths(6) = EventLengths(6) + dt1;
% EventLengths(8) = EventLengths(8) - dt1;

% dt2 = 30;
% EventLengths(4) = EventLengths(4) + dt2;
% EventLengths(6) = EventLengths(6) - dt2;


[Sequence, T] = reorder_events(EventLengths,isPulse);