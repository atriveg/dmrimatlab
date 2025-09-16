function [sf,sfname,vs] = check_software_platform
% function [sf,sfname,vs] = check_software_platform
%
%   Checks if the software is running under Mathwork's Matlab of GNU
%   Octave, and returns the version of the software.
%
%   OUTPUTS:
%
%      sf: either 1 for Matlab, 2 for Octave or -1 for other unknown
%          platforms.
%      sfname: either 'matlab', 'octave' or 'unknown'.
%      vs: a double with the version of the software, the integer part
%          containing the major version and the decimals the minor version
v = ver('Matlab');
if(isempty(v))
    v = ver('Octave');
    if(isempty(v))
        sf     = -1;
        vs     = NaN;
        sfname = 'unknown';
        return;
    else
        sf     = 2;
        sfname = 'octave';
    end
else
    sf     = 1;
    sfname = 'matlab';
end
vstr  = v.Version;
delim = strfind(vstr,'.');
if(length(delim)>1)
    vstr = vstr(1:delim(2)-1);
end
vs = str2double(vstr);
end
