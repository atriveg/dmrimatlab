function makefile_mexcode(varargin)
% This is a function that builds all (or only selected) modules of mex code
% within the toolbox. Use cases are:
%
%   + Print the list of all available modules that can be built:
%
%     >> makefile_mexcode list
%
%   + Check which modules need to be built:
%
%     >> makefile_mexcode missing
%
%   + Check which modules among a given list need to be built:
%
%     >> makefile_mexcode missing module1 module2 ... moduleN
%
%     where "modulei" are names within the list provided with the "list"
%     command (first use case).
%
%   + Build ALL modules, no matter whether they are already built or not:
%
%     >> makefile_mexcode
%     >> makefile_mexcode build
%
%   + Build only certain modules:
%
%     >> makefile_mexcode module1 module2 ... moduleN
%     >> makefile_mexcode build module1 module2 ... moduleN
%
%     where "modulei" are names within the list provided with the "list"
%     command (first use case).
%
%   + Delete ALL modules, no matter whether they are exist or not:
%
%     >> makefile_mexcode clean
%
%   + Delete only certain modules:
%
%     >> makefile_mexcode clean module1 module2 ... moduleN
%
%     where "modulei" are names within the list provided with the "list"
%     command (first use case).
%
% NOTE: if using windows, either MinGW or Cywin/gcc are strongly
%       recommended. We have not tested other compilers.

sf = check_software_platform;
if(sf==1)
    makefile_mexcode_matlab(varargin{:});
elseif(sf==2)
    makefile_mexcode_octave(varargin{:});
end

end
