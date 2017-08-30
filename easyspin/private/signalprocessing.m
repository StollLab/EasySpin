function [ ProcessedSignal ] = signalprocessing(TimeAxis,RawSignal,Opt)
% Number of detection operators
nDetectionOperators = length(Opt.DetectionOperators);

% Setup vector with down conversion frequencies
TranslationFrequencies = zeros(1,nDetectionOperators);

% Trys down conversion, if it fails for any reason, the raw signal is
% returned. Errors are usually when the down conversion frequency is very
% wrong, or a non oscillating signal is to be down converted (Sz)
try
  nDownConversionFrequencies = length(Opt.FreqTranslation);
  TranslationFrequencies(1:nDownConversionFrequencies) = Opt.FreqTranslation;
  
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
        RfSignal = RawSignal{iPoint}(iTrace,:);
      else
        RfSignal = squeeze(RawSignal(iPoint,iTrace,:));
      end
      
      % Ensures that signal is purely real or imaginary, if it is supposed 
      % to be purely real/imaginary
      if max(abs(RfSignal))/max(imag(RfSignal)) > 1e8
        RfSignal = real(RfSignal);
        purelyImag = false;
      elseif max(abs(RfSignal))/max(real(RfSignal)) > 1e8
        RfSignal = imag(RfSignal);
        purelyImag = true;
      else
        % If signal is complex (S+ or S-), the normalization is wrong by a
        % factor of two and hence needs normalization
        RfSignal = 2*RfSignal;
      end
      % Break signal into smaller signals if a nondetected event is present
      % (non linear time axis)
      % store indices
          
      num_dig = 10;

      diffTime = round(diff(DCTimeAxis)*(10^num_dig))/(10^num_dig);
      dt = min(diffTime);
      
      Opt.dt = dt;
      
      BreakIndices = [1  find(diffTime ~= dt)+1 length(DCTimeAxis)+1];
      % Does the downconversion, if no down conversion is requested for the
      % current detection operator, the cleaned up signal is written to
      % ProcessedSignal
      
      if TranslationFrequencies(iTrace) ~= 0
        DCSignal = zeros(1,length(RfSignal));
        % Depending on the type of trace, different down conversion types
        % are required (see rfmixer for details)
        for j = 1 : length(BreakIndices)-1
          Elements = BreakIndices(j):(BreakIndices(j+1)-1);
          if isreal(RfSignal)
            % Mixing if signal is real
            % loop over indices for broken (or not broken) down signal and
            % provide timestep
            
            [~, DCSignal(Elements)] = rfmixer(DCTimeAxis(Elements),RfSignal(Elements),TranslationFrequencies(iTrace),'IQdemod',Opt);
            
            if purelyImag
              DCSignal = real(DCSignal)*1i;
            else
              DCSignal = real(DCSignal);
            end
            
          else
            % Mixing if signal is complex
            
            % loop over indices for broken (or not broken) down signal and
            % provide timestep
            
            [~, DCSignal(Elements)] = rfmixer(DCTimeAxis(Elements),RfSignal(Elements),TranslationFrequencies(iTrace),'IQshift',Opt);
            
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
          Traces(iTrace,:) = RfSignal;
        else
          ProcessedSignal(iPoint,iTrace,:) = RfSignal;
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

