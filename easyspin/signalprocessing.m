% signalprocessing    Signal translation and clean up for spidyan
%
%     x = signalprocessing(TimeAxis,RawSignal,FreqTranslation)
% 
%     TimeAxis   ... time axis of the input signal in microseconds
%     RawSignal   ... original signal, as returned by saffron and spidyan
%     FreqTranslation   ... vector with frequencies in GHz for translation
%
%     out:
%       x          ... signal(s) after frequency translation

function [ProcessedSignal] = signalprocessing(TimeAxis,RawSignal,FreqTranslations)

nDetectionOperators = length(FreqTranslations);

% Tries down conversion, if it fails for any reason, the raw signal is
% returned. Errors are usually when the down conversion frequency is very
% wrong, or a non oscillating signal is to be down converted (Sz)
try
  % Recognizing the type of the input, whether it is a cell or numeric 
	% array and gets the number of acquisition points
  % For the cell array each element has to processed individually which is
  % much more straightforward
  if iscell(RawSignal)
    ProcessedSignal = cell(size(RawSignal));
    nPoints = numel(RawSignal);
  else
    % for a numeric array, several tests need to be done
    SignalSize = size(RawSignal);
    % First a check of the dimensionality of the input, necessary since
    % signalprocessing could be called by a user (and not from within
    % spidyan) where the data might have already been reduced from a
    % multi-dimensional array to a two dimensional array (e.g. only one
    % acquistion point)
    if ndims(RawSignal) == 2 %#ok<ISMAT>
      nPoints = 1;
      % turn it into a 3D array, because thats how they are processed
      if SignalSize(end) == nDetectionOperators
        RawSignal = reshape(RawSignal,[1 1 SignalSize]);
      else
        RawSignal = reshape(RawSignal,[1 SignalSize(end-1) SignalSize(end)]);
      end
    else
      nPoints = prod(SignalSize(1:end-2));
      % reshape the n-dimensional array into a 3-dimensional array, that
      % can be looped over linearly along the first dimension (which
      % corresponds to all the acquistion points)
      if SignalSize(end) == nDetectionOperators
        RawSignal = reshape(RawSignal,[nPoints SignalSize(end) SignalSize(end-1)]);
      else
        RawSignal = reshape(RawSignal,[nPoints SignalSize(end-1) SignalSize(end)]);
      end
    end
    
    ProcessedSignal = zeros(size(RawSignal));
    
    % reshape time axis if simulation had more than one indirect dimension
    % and/or time axes are not identical
    if ndims(TimeAxis) > 2 %#ok<ISMAT>
      TimeAxis = reshape(TimeAxis,[nPoints,SignalSize(end)]);
    end
  end
  
  % loop over all acquired data points
  for iPoint = 1 : nPoints
    if iscell(RawSignal)
      % Gets the size of the traces for the current data point and the time
      % axis if the RawSignal is a CellArray
      SignalSize = size(RawSignal{iPoint});
      if SignalSize(2) == nDetectionOperators
        RawSignal{iPoint} = reshape(RawSignal,[SignalSize(2) SignalSize(1)]);
      end
      Traces = zeros(size(RawSignal{iPoint}));
      DCTimeAxis = TimeAxis{iPoint};
    else
      % If a time axis for each acquistion point was provided, the correct
      % time axis is loaded for the downconversion
      if size(TimeAxis,1) > 1
        DCTimeAxis = TimeAxis(iPoint,:);
      else
        % if all signals have the same time axis
        DCTimeAxis = TimeAxis;
      end
    end
    % Loop over the Number of Detection Operators
    for iTrace = 1 : nDetectionOperators
      % Load RawSignal of the current Detection Operator
      if iscell(RawSignal)
        RFSignal = RawSignal{iPoint}(iTrace,:);
      else
        RFSignal = squeeze(RawSignal(iPoint,iTrace,:));
      end
      
      % Ensures that signal is purely real or imaginary, if it is supposed 
      % to be purely real/imaginary
      if max(abs(RFSignal))/max(imag(RFSignal)) > 1e8
        RFSignal = real(RFSignal);
        purelyImag = false;
      elseif max(abs(RFSignal))/max(real(RFSignal)) > 1e8
        RFSignal = imag(RFSignal);
        purelyImag = true;
      else
        % If signal is complex (S+ or S-), the normalization is wrong by a
        % factor of two and hence needs normalization, but only do it, if
        % the signal is also being demodulated. Else, if signalprocessing
        % is called manually at a later point, the signal gets multiplied
        % too many times
        if FreqTranslations(iTrace) ~= 0
          RFSignal = 2*RFSignal;
        end
      end
      % Break signal into smaller signals if a nondetected event is present
      % (non linear time axis)
      % store indices
          
      % Looks for nonlinear steps in the time axis (with rounding), if
      % jumps in the time axis are detected, the signal is processed
      % separatel for each continuuous section. This is the case when
      % detection several events are detected, with one in the middle that
      % is not, eg. Det.Events = [1 0 1]
      num_dig = 10;
      diffTime = round(diff(DCTimeAxis)*(10^num_dig))/(10^num_dig);
      
      % Extracting the time step
      dt = min(diffTime);
      Opt.dt = dt;
      
      % Indices for splitting the signal
      BreakIndices = [1  find(diffTime ~= dt)+1 length(DCTimeAxis)+1];
      
      
      % Does the downconversion, if no down conversion is requested for the
      % current detection operator, the cleaned up signal is written to
      % ProcessedSignal
      
      if FreqTranslations(iTrace) ~= 0
        DCSignal = zeros(1,length(RFSignal));
        % Depending on the type of trace, different down conversion types
        % are required (see rfmixer for details)
        for j = 1 : length(BreakIndices)-1
          Elements = BreakIndices(j):(BreakIndices(j+1)-1);
          if isreal(RFSignal)
            % Mixing if signal is real or imaginary
            [~, DCSignal(Elements)] = rfmixer(DCTimeAxis(Elements),RFSignal(Elements),FreqTranslations(iTrace),'IQdemod',Opt);
            
            if purelyImag
              DCSignal = real(DCSignal)*1i;
            else
              DCSignal = real(DCSignal);
            end
            
          else
            % Mixing if signal is complex          
            [~, DCSignal(Elements)] = rfmixer(DCTimeAxis(Elements),RFSignal(Elements),FreqTranslations(iTrace),'IQshift',Opt);
            
          end
        end
        
        if iscell(RawSignal)
          Traces(iTrace,:) = DCSignal;
        else
          ProcessedSignal(iPoint,iTrace,:) = DCSignal;
        end
        
      else
        % stores original cleaned up signal
        if iscell(RawSignal)
          Traces(iTrace,:) = RFSignal;
        else
          ProcessedSignal(iPoint,iTrace,:) = RFSignal;
        end
        
      end
    end
    if iscell(RawSignal)
      ProcessedSignal{iPoint} = Traces.';
    end
  end
  
  if ~iscell(RawSignal)
    ProcessedSignal = reshape(ProcessedSignal,SignalSize);
    nDimsSignal = length(SignalSize);
    ProcessedSignal = permute(ProcessedSignal,[1:nDimsSignal-2,nDimsSignal,nDimsSignal-1]);
  end
catch EM
  % If something goes wrong during down conversion, the RawSignal is
  % returned
  ProcessedSignal = RawSignal;
  message = ['Down conversion of the signal was not successful and created an error and the raw signal in the simulation frame was returned. The encountered error was: ' EM.message];
  warning(message)
end

% If only one acquisition point, singleton dimension is removed
if ~iscell(ProcessedSignal)
  if ndims(ProcessedSignal) == 3 && size(ProcessedSignal,1) == 1
    ProcessedSignal = squeeze(ProcessedSignal);
  end
  if ndims(ProcessedSignal) == 2 && size(ProcessedSignal,1) == 1 && nDetectionOperators == 1
    ProcessedSignal = permute(ProcessedSignal,[2 1]);
  end
end

