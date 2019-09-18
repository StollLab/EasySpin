function [err,data] = test(opt,olddata)

EventLengths1 = [10 50 10 50 10 50];
PulseList = [1 0 1 0 1 0 1];
Sequence1 = 1:6;

[NewSequence1, NewEventLengths1] = runprivate('s_reorder_events',EventLengths1,PulseList);


EventLengths2 = [10 50 10 -30 10 50];
Sequence2 = [1 2 5 6 3 4];

[NewSequence2, NewEventLengths2] = runprivate('s_reorder_events',EventLengths2,PulseList);
CorrectEventLengths = [10 30 10 10 10 30];


if any([~isequal(EventLengths1,NewEventLengths1) ~isequal(Sequence1,NewSequence1) ...
    ~isequal(CorrectEventLengths,NewEventLengths2) ~isequal(Sequence2,NewSequence2)])
  err = 1;
else
  err = 0;
end

data = [];

