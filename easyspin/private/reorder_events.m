function [NewSequence, NewEventLenghts] = reorder_events(EventLengths,PulseList)

nEvents = numel(EventLengths);
t = [0 cumsum(EventLengths)];

% dy = 1;
% figure
% for k = 1:nEvents
%   if isPulse(k), col = 'r'; else, col = 'k'; end
%   line([1 1]*t(k),[1 nEvents],'Color',[0 0.5 0]);
%   h = line([t(k) t(k+1)],k*dy*[1 1],'Color',col,'LineWidth',3);
%   if k<nEvents
%     line([1 1]*t(k+1),[k,k+1]*dy,'Color',[1 1 1]*0);
%   end
% end

%%
[TimeIntervals,NewSequence] = sort(t);
NewSequence = NewSequence(1:end-1);

NewEventLenghts = zeros(size(EventLengths));
for k = 1:nEvents
NewEventLenghts(k) = TimeIntervals(k+1) - TimeIntervals(k) ;
end



%%

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

% NewSequence = zeros(size(isPulse));
% 
% for iSequence = 1 : size(M,2)
%   for iOldEvent = 1 : size(M,1)
%     if M(iOldEvent,iSequence) == 2
%       NewSequence(iSequence) = iOldEvent;
%       continue
%     elseif M(iOldEvent,iSequence) == 1 && NewSequence(iSequence) == 0 && ~any(NewSequence==iOldEvent)
%       NewSequence(iSequence) = iOldEvent;
%     end
%   end
% end

% if NewSequence(end) == 0
%   if any(NewSequence == nEvents)
%     error('Something went wrong.')
%   else
%     NewSequence(end) = nEvents;
%   end
% end
% 
% dy = 1;
% figure; clf
% for k = 1:nEvents
%   if PulseList(k), col = 'r'; else, col = 'k'; end
%   line([1 1]*TimeIntervals(k),[1 nEvents],'Color',[0 0.5 0]);
%   h = line([TimeIntervals(k) TimeIntervals(k+1)],NewSequence(k)*dy*[1 1],'Color',col,'LineWidth',3);
%   if k<nEvents
%     line([1 1]*TimeIntervals(k+1),[k,k+1]*dy,'Color',[1 1 1]*0);
%   end
% end

% Assert absence of pulse overlap
if any(sum(EventMatrix==2,1)>1)
%    message = 'Pulse Overlap!'
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

% isEvolveInterval = any(M==1) & all(M~=2);

% for iInterval = 1:nEvents
%   if ~isEvolveInterval, continue; end
%   % Check for conflicts between free-evolution events for this interval
%   for iEvent = 1:nEvents
%     if M(iEvent,iInterval)~=1, continue; end
%     
%   end
% end

end

