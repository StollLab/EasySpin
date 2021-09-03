%{
Simple wrapper class to run EasySpin's function-based tests through MATLAB's
testing framework.

Suggested by Andy Campbell, MathWorks, 9/2/2021

To run the tests, run the following commands from the easyspin root folder:

import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
import matlab.unittest.plugins.TestReportPlugin

runner = TestRunner.withTextOutput("OutputDetail", "Detailed")
runner.addPlugin(CodeCoveragePlugin.forFolder('easyspin',"Producing", CoverageReport))
runner.addPlugin(TestReportPlugin.producingPDF)

runner.run(testsuite("tests"))

%}

classdef LegacyTests < matlab.unittest.TestCase
  
  properties(TestParameter)
    fcn = getListOfTests;
  end
  
  methods(Test)
    
    function test(testCase, fcn)
      testFcn = str2func(fcn);
      nArgsOut = nargout(testFcn);
      nArgsIn = nargin(testFcn);
      
      % Apply options
      fcnInputs = {};
      if nArgsIn > 0
        opt.Display = false;
        opt.Regenerate = false;
        opt.Verbosity = false;
        fcnInputs{end+1} = opt;
      end
      
      % Get baseline data
      if nArgsIn==2 && nArgsOut==2
        dataLocation = fileparts(which(fcn));
        olddata = load(fullfile(dataLocation, "data", fcn + ".mat"));
        fcnInputs{end+1} = olddata.data;
      end
      
      % Run the function
      out = cell(1,nArgsOut);
      [out{:}] = testFcn(fcnInputs{:});
      testCase.verifyThat(EveryElementOf(out{1}), IsTrue, ...
        sprintf("Invoking the ""%s"" test should return true.", fcn));
    end
    
  end
  
end

function tests = getListOfTests()

% Get the test folder & its info
[folder, ~] = fileparts(which(mfilename));
d = dir(fullfile(folder,"*_*.m"));

% Grab all the files
tests = {d.name};

% Remove the extension
tests = erase(tests, '.m');

end

% "Import" functions
function p = EveryElementOf(varargin)
p = matlab.unittest.constraints.EveryElementOf(varargin{:});
end

function c = IsTrue(varargin)
c = matlab.unittest.constraints.IsTrue(varargin{:});
end
