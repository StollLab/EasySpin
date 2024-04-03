% rcfilt  Apply an RC filter to a signal
%
%  yFilt = rcfilt(y,SampTime,TimeConstant)
%  yFilt = rcfilt(y,SampTime,TimeConstant,'up')
%  yFilt = rcfilt(y,SampTime,TimeConstant,'down')
%
%  Filters a spectrum using a RC low-pass
%  filter as built into cw spectrometers
%  to remove high-frequency noise.
%
%  Input:
%  - y              unfiltered spectrum
%  - SampleTime     sampling time
%  - TimeConstant   time constant of the filter
%  - 'up' or 'down' defines the direction of the
%       field sweep. If omitted, 'up' is assumed.
%
%    SampleTime and TimeConstant must have the same
%    units (s, ms, µs, ...)
%
%  Output:
%  - yFilt          filtered spectrum
%
%  For matrices, rcfilt operates along columns.

function varargout = rcfilt(y,SampleTime,TimeConstant,UpDown)

if nargin==0, help(mfilename); return; end
if nargin<4, UpDown = 'up'; end

% Check input arguments
err = '';
if ~isreal(TimeConstant), err = 'Time constant must be real!'; end
if ~isreal(SampleTime), err = 'Sampling time must be real!'; end
if numel(TimeConstant)~=1, err = 'Time constant must be a scalar!'; end
if numel(SampleTime)~=1, err = 'Sampling time must be a scalar!'; end
if TimeConstant<0, err = 'Time constant must be positive or zero!'; end
if SampleTime<=0, err = 'Sampling time must be positive!'; end
if all(~strcmp(UpDown,{'up','down','dn'})), err = 'Specify either ''up'' or ''down''.'; end
error(err);

RowVector = isrow(y);
if RowVector
  y = y.';
end

Invert = UpDown=="dn" || UpDown=="down";
if Invert
  y = y(end:-1:1,:);
end

% Apply low-pass RC filter
if TimeConstant==0
  yFiltered = y;
else
  e = exp(-SampleTime/TimeConstant);
  for iCol = size(y,2):-1:1
    yFiltered(:,iCol) = filter(1-e,[1 -e],y(:,iCol));
  end
end

if Invert
  y = y(end:-1:1,:);
  yFiltered = yFiltered(end:-1:1,:);
end
if RowVector
  y = y.';
  yFiltered  = yFiltered.';
end

switch nargout
  case 0
    x = 1:size(y,1);
    plot(x,y,'b',x,yFiltered,'g');
    legend('raw','filtered');
    xlabel('data point #');
    ylabel('value');
    title('Effect of RC filter on data');
    varargout = {};
  case 1
    varargout = {yFiltered};
  otherwise
    error('Wrong number of output arguments.');
end

end
