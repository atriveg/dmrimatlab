function printBLASInfo
% function printBLASInfo
%
%   Prints the versions of BLAS and LAPACK Octave itself and the toolbox
%   mex-files are linked against.

[sf,sfname,vs] = check_software_platform;
if(sf~=2)
    error('This toolbox function works only in GNU Octave under Linux');
end

% Retrieve information on BLAS/LAPACK libraries:
fprintf(1,'BLAS and LAPACK used by Octave:\n');
fprintf(1,'   BLAS: %s\n',version('-blas'));
fprintf(1,'   LAPACK: %s\n',version('-lapack'));
if( exist ('mexLinkInfo.mex','file')==3 )
   fprintf(1,'BLAS and LAPACK used by mex files:\n');
   mexLinkInfo;
else
   warning('Mex-files are not built yet. Please run makefile_mexcode');
end

end
