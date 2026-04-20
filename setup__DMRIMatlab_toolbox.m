function setup__DMRIMatlab_toolbox(varargin)
% function setup__DMRIMatlab_toolbox('opt1',value1,'opt2',value2,...)
%
%   Setup the dMRI toolbox for matlab by fixing the path and configuring
%   options. By now, options (given as regular matlab's 'option', value
%   pairs are:
%
%     usebroadcast: whhether (true) or not (false) use Python-like
%        broadcasted operations instead of bsxfun for large matrixes. By
%        default, broadcast is used if available, i.e. if your Matlab
%        or Octave version supports it.
%     useparallel: wether (true) or not (false) use the parallel computing
%        toolbox to speed up certain computations, mainly voxel-by-voxel
%        matrix inversions or voxel-by-voxel eigenvalues computations
%        (default: false).

% First of all, check if we deal with Matlab or
% Octave, and set proper paths accordingly:
dmrilabvers = '1.0';
if(~isdeployed)
    % ---
    path0       = mfilename('fullpath');
    [path0,~,~] = fileparts(path0);
    % ---
    parse__path_dmri( sprintf('%s%sutils',path0,filesep), {} );
    [sf,sfname,vs] = check_software_platform;
    msg1 = sprintf('You''re running version %s of the toolbox under %s %1.2f\n',dmrilabvers,sfname,vs);
    % ---
    if(sf==1) % i.e. Matlab
        forbidden = {[path0,filesep,'graphics',filesep,'octave']};
    elseif(sf==2) % i.e Octave
        forbidden = {};
    else % i.e. Unknown
        error('Unknown software platform');
    end
    parse__path_dmri(path0,forbidden);
end


% Check options:
% -------------------------------------------------------------------------
opt.usebroadcast = true; optchk.usebroadcast = [true,true]; % always 1x1 boolean
opt.useparallel = false; optchk.useparallel = [true,true];  % always 1x1 boolean
opt.quiet = false;       optchk.quiet = [true,true];
opt.nographics = false;  optchk.nographics = [true,true];
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------

if(~isdeployed)
    if(~opt.nographics)
        if(sf==2)
            % Select graphics toolkit:
            tkts = available_graphics_toolkits;
            if( any(strcmp(tkts,'qt')) )
                try
                    graphics_toolkit('qt');
                    msg2 = sprintf('Using graphics toolkit: qt\n');
                catch
                    warning('Unable to set graphics toolkit to qt. Using default');
                end
            elseif( any(strcmp(tkts,'gnuplot')) )
                try
                    graphics_toolkit('gnuplot');
                    msg2 = sprintf('Using graphics toolkit: qt\n');
                catch
                    warning('Unable to set graphics toolkit to gnuplot. Using default');
                end
            elseif( any(strcmp(tkts,'fltk')) )
                try
                    graphics_toolkit('fltk');
                    msg2 = sprintf('Using graphics toolkit: fltk\n');
                catch
                    warning('Unable to set graphics toolkit to fltk. Using default');
                end
            else
                warning('Unable to auto-detect graphics toolkit');
            end
        end
    end
end

if(~opt.quiet)
    if(exist('msg1','var'))
        fprintf(1,msg1);
    end
    if(exist('msg2','var'))
        fprintf(1,msg2);
    end
end

global is_broadcast_available_test_var; %#ok<GVMIS>
if(~opt.usebroadcast)
    is_broadcast_available_test_var = false;
else
    is_broadcast_available_test_var = [];
    is_broadcast_available_test;
end

global use_parallel_test_var; %#ok<GVMIS>
if(opt.useparallel)
    if(sf==2) % parfors are not actually implemented in Octave
        warning('parfors are not currently supported in Octave. Ignoring ''useparallel''');
        use_parallel_test_var = false;
    else
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
