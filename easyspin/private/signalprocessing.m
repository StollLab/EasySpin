function [ ProcessedSignal ] = signalprocessing(TimeAxis,RawSignal,DetectionOperators,FreqTranslation)
% Number of detection operators
nDetectionOperators = length(DetectionOperators);

% Setup vector with down conversion frequencies
TranslationFrequencies = zeros(1,nDetectionOperators);

% Trys down conversion, if it fails for any reason, the raw signal is
% returned. Errors are usually when the down conversion frequency is very
% wrong, or a non oscillating signal is to be down converted (Sz)
try
  nDownConversionFrequencies = length(FreqTranslation);
  TranslationFrequencies(1:nDownConversionFrequencies) = FreqTranslation;
  
  % Recognizing the type of the input, wheter it is a cell or numeri array
  % and gets the number of acquisition points
  % For the cell array each element has to processed individually
  if iscell(RawSignal)
    ProcessedSignal = cell(size(RawSignal));
    nPoints = length(RawSignal);
  else
    ProcessedSignal = zeros(size(RawSignal));
    nPoints = size(RawSignal,1);
    DCTimeAxis = TimeAxis(1,:);
  end
  
  % loop over all acquired data points
  for iPoint = 1 : nPoints
    if iscell(RawSignal)
      % Gets the size of the traces for the current data point and the time
      % axis if the RawSignal is a CellArray
      Traces = zeros(size(RawSignal{iPoint}));
      DCTimeAxis = TimeAxis{1,iPoint};
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
        % factor of two and hence needs normalization
        RFSignal = 2*RFSignal;
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
      
      if TranslationFrequencies(iTrace) ~= 0
        DCSignal = zeros(1,length(RFSignal));
        % Depending on the type of trace, different down conversion types
        % are required (see rfmixer for details)
        for j = 1 : length(BreakIndices)-1
          Elements = BreakIndices(j):(BreakIndices(j+1)-1);
          if isreal(RFSignal)
            % Mixing if signal is real or imaginary
            [~, DCSignal(Elements)] = rfmixer(DCTimeAxis(Elements),RFSignal(Elements),TranslationFrequencies(iTrace),'IQdemod',Opt);
            
            if purelyImag
              DCSignal = real(DCSignal)*1i;
            else
              DCSignal = real(DCSignal);
            end
            
          else
            % Mixing if signal is complex          
            [~, DCSignal(Elements)] = rfmixer(DCTimeAxis(Elements),RFSignal(Elements),TranslationFrequencies(iTrace),'IQshift',Opt);
            
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
      ProcessedSignal{iPoint} = Traces;
    end
  end
  
catch EM
  % If something goes wrong during down conversion, the RawSignal is
  % returned
  ProcessedSignal = RawSignal;
  message = ['The down conversion of the signal was not successful and created an error. The raw signal in the simulation frame was returned. The error message was: ' EM.message];
  warning(message)
end

end

