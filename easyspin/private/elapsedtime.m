% Computes elapsed time between StartTime and EndTime (in seconds)
% in hours, minutes and seconds.

function [hmsString,Hours,Minutes,Seconds] = elapsedtime(startTime,endTime)

Time = etime(endTime,startTime)/3600; % convert seconds to hours
Hours = fix(Time);
remainingMinutes = (Time - Hours)*60;
Minutes = fix(remainingMinutes);
Seconds = (remainingMinutes - Minutes)*60;

if Hours>0
  hmsString = sprintf('%dh%dm%0.3fs',Hours,Minutes,Seconds);
elseif Minutes>0
  hmsString = sprintf('%dm%0.3fs',Minutes,Seconds);
else
  hmsString = sprintf('%0.3fs',Seconds);
end

end
