function [NewSequence, NewEventLengths] = s_reorder_events(EventLengths,PulseList)

nEvents = numel(EventLengths);

% Transfers intervals to positions - absolute times when an event is to start
t = [0 cumsum(EventLengths)];

%-------------------------------------------------------------------------------
% Plots of how the time periods are defined before reordering (for debugging)
%-------------------------------------------------------------------------------
plotTimeLine = false;
if plotTimeLine
  dy = 1;
  figure
  for k = 1:nEvents
    if isPulse(k), col = 'r'; else, col = 'k'; end
    line([1 1]*t(k),[1 nEvents],'Color',[0 0.5 0]);
    h = line([t(k) t(k+1)],k*dy*[1 1],'Color',col,'LineWidth',3);
    if k<nEvents
      line([1 1]*t(k+1),[k,k+1]*dy,'Color',[1 1 1]*0);
    end
  end
end
%-------------------------------------------------------------------------------

% Sort the events according to their starting position
[TimeIntervals,NewSequence] = sort(t);
NewSequence = NewSequence(1:end-1);

% Recalculate intervals/event durations from the new starting positions
NewEventLengths = zeros(size(EventLengths));
for k = 1:nEvents
  NewEventLengths(k) = TimeIntervals(k+1) - TimeIntervals(k) ;
end

% Rounding is necessary because MATLAB does loses precision during the
% previous step, on the order of numeric precision (10^-16)
num_dig = 10;
NewEventLengths = round(NewEventLengths*(10^num_dig))/(10^num_dig);

% Creates an event matrix that correlates events as they were indixed when
% input to the sequence, values of two correspond to pulses, one to free
% evolution events
EventMatrix = zeros(nEvents,nEvents);
for iEvent = 1 : nEvents
  t_event = [t(iEvent) t(iEvent)+EventLengths(iEvent)];
  tEventStart = min(t_event);
  tEventEnd = max(t_event);
  idx = TimeIntervals>=tEventStart & TimeIntervals<tEventEnd;
  if PulseList(iEvent)==1
    EventMatrix(iEvent,idx) = 2;
  else
    EventMatrix(iEvent,idx) = 1;
  end
end

% Assert absence of pulse overlap, if pulse overlap is detected, an error
% is returned which tells the user the index of the overlapping events.
if any(sum(EventMatrix==2,1)>1)
  Overlaps = find(sum(EventMatrix==2,1)>1);
  OverlappingEvents = sum(EventMatrix(:,Overlaps)==2,2);
  OverlappingEvents = find(OverlappingEvents>0);
  ErrorMessage = 'Pulse Overlap! Events: ';
  ErrorMessage = [ErrorMessage num2str(OverlappingEvents(1)) ', ' num2str(OverlappingEvents(2))];
  if length(Overlaps) > 2
    for i = 3 : length(Overlaps)
      ErrorMessage = [ErrorMessage ', ' num2str(OverlappingEvents(i))];
    end
  end
  error(ErrorMessage);
end

% When a pulse and an event start at the same (absolute) time, the
% reordering fails. This loop takes care that in such a case the pulse is
% exchanged with the next free evolution event. Repeated for all pulses
% with zero pulse length.
for iEvent = 1 : nEvents
  if PulseList(NewSequence(iEvent)) && NewEventLengths(iEvent) == 0
    ShortenedPulse = NewSequence(iEvent);
    NewSequence(iEvent) = NewSequence(iEvent+1);
    NewSequence(iEvent+1) = ShortenedPulse;
  end
end

%-------------------------------------------------------------------------------
% Plots of how the time periods are defined after reordering (for debugging)
%-------------------------------------------------------------------------------
if plotTimeLine
  dy = 1;
  figure
  for k = 1:nEvents
    if PulseList(k), col = 'r'; else, col = 'k'; end
    line([1 1]*TimeIntervals(k),[1 nEvents],'Color',[0 0.5 0]);
    h = line([TimeIntervals(k) TimeIntervals(k+1)],NewSequence(k)*dy*[1 1],'Color',col,'LineWidth',3);
    if k<nEvents
      line([1 1]*TimeIntervals(k+1),[k,k+1]*dy,'Color',[1 1 1]*0);
    end
  end
end
%-------------------------------------------------------------------------------

end
