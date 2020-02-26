function ok = test()

PulseList = [1 0 1 0 1 0 1];

EventLengths = [10 50 10 -20 10 50];
CorrectSequence = [1 2 5 6 3 4];

[NewSequence, NewEventLengths] = runprivate('s_reorder_events',EventLengths,PulseList);
CorrectEventLengths = [10 40 10 0 10 40];

ok = all([isequal(NewEventLengths,CorrectEventLengths) isequal(NewSequence,CorrectSequence)]);


