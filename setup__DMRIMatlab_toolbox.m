function setup__DMRIMatlab_toolbox(varargin)
% function setup__DMRIMatlab_toolbox('opt1',value1,'opt2',value2,...)
%
%   Setup the dMRI toolbox for matlab by fixing the path and configuring
%   options. By now, options (given as regular matlab's 'option', value
%   pairs are:
%
%     usebroadcast: wether(true) or not (false) use Python-like broadcasted
%        operations instead of bsxfun for large matrixes. By default,
%        broadcast is used if available, i.e. if your Matlab version
%        supports it.
%     useparallel: wether (true) or not (false) use the parallel computing
%        toolbox to speed up certain computations, mainly voxel-by-voxel
%        matrix inversions or voxel-by-voxel eigenvalues computations
%        (default: false).

% First of all, set path:
if(~isdeployed)
    path0       = mfilename('fullpath');
    [path0,~,~] = fileparts(path0);
    forbidden = {};
    %forbidden   = { [path0,filesep,'bib'],...
    %    [path0,filesep,'MTIS',filesep,'ToMatlabBinWin32'] };
    parse__path_dmri(path0,forbidden);
end

% Check options:
% -------------------------------------------------------------------------
opt.usebroadcast = true; optchk.usebroadcast = [true,true]; % always 1x1 boolean
opt.useparallel = false; optchk.useparallel = [true,true];  % always 1x1 boolean
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------

global is_broadcast_available_test_var;
if(~opt.usebroadcast)
    is_broadcast_available_test_var = false;
else
    is_broadcast_available_test_var = [];
    is_broadcast_available_test;
end

global use_parallel_test_var;
if(opt.useparallel)
    % User requests using the parallel computing toolbox. Check if this is
    % possible:
    version = ver('parallel');
    if(~isempty(version))
        if(license('test','Distrib_Computing_Toolbox'))
            % Installed and licensed
            use_parallel_test_var = true;
            % Init the pool if it does not exist yet
            gcp;
        else
            % Not licensed
            use_parallel_test_var = false;
        end
    else
        % Not installed
        use_parallel_test_var = false;
    end
else
    use_parallel_test_var = false;
end

% -------------------------------------------------------------------------
function parse__path_dmri(path,forbidden)
addpath(path);
list = dir(path);
for n=1:length(list)
    if(list(n).isdir)
        if( (list(n).name(1)~='.') && ((list(n).name(1)~='@')) )
            nestedpath = [path,filesep,list(n).name];
            if(~any(strcmp(nestedpath,forbidden)))
                parse__path_dmri(nestedpath,forbidden);
            end
        end
    end
end
