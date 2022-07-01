function ok = test()
T = validateFunctionSignaturesJSON("..\easyspin\functionSignatures.json") ;
ok = isempty(T) ;
