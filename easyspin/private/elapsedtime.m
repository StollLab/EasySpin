% Computes elapsed time between StartTime and EndTime (in seconds)
% in hours, minutes and seconds.

function [Hours,Minutes,Seconds] = elapsedtime(StartTime,EndTime)

Time = etime(EndTime,StartTime)/3600; % convert seconds to hours
Hours = fix(Time);
RemainingMinutes = (Time - Hours)*60;
Minutes = fix(RemainingMinutes);
Seconds = (RemainingMinutes - Minutes)*60;

return
