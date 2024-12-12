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

if(nargin>0)
    if(   strcmpi( varargin{1}, 'clean' )   )
        action = 'c';
    elseif(   strcmpi( varargin{1}, 'build' )   )
        action = 'b';
    elseif(   strcmpi( varargin{1}, 'list' )   )
        action = 'l';
    elseif(   strcmpi( varargin{1}, 'missing' )   )
        action = 'm';
    else
        action = 'x';
    end
else
    action = 'b';
end

% ------------------------------------------------------------------------------
mtlroot = matlabroot;

if(isunix)
    trlinks = {'-lpthread'};
    %%%
    libsdir = sprintf('%s/bin/glnxa64',mtlroot);
    mkllib  = sprintf('%s/mkl.so',libsdir);
    blaslinks = {'-lmwlapack','-lmwblas',mkllib};
    blasflags = {'CXXFLAGS="$CXXFLAGS -D_USE_MKL_THREAD_CONTROL"'};
end

if(ispc)
    trlinks = {};
    %%%
    blaslinks = {'-lmwlapack','-lmwblas'};
    blasflags = {'CXXFLAGS="$CXXFLAGS -fopenmp -D_USE_OMP_THREAD_CONTROL"','LDFLAGS="$LDFLAGS -fopenmp"'};
end

if( ismac )
    trlinks = {'-lpthread'};
    % Check if mkl.dylib is present. If it is the case, this is an old
    % Mac machine with an Intel processor, hence we can rely on MKL as
    % we do for linux (but we must explicitly tell the computer through
    % the -D_HAS_MKL_BLAS flag).
    libsdir = sprintf('%s/bin/maci64',mtlroot);
    mkllib  = sprintf('%s/mkl.dylib',libsdir);
    if( exist(mkllib,'file') ~= 0 )
        % MAC with intel processor
        blaslinks = {'-lmwlapack','-lmwblas',mkllib};
        blasflags = {'CXXFLAGS="$CXXFLAGS -D_USE_MKL_THREAD_CONTROL"'};
    else
        % MAC with ARM processor
        blaslinks = {'-lmwlapack','-lmwblas','-lmwopenblas'};
        blasflags = {'CXXFLAGS="$CXXFLAGS -D_USE_OPENBLAS_THREAD_CONTROL"'};
    end
end
% ------------------------------------------------------------------------------
oldpath = pwd;
% ------------------------------------------------------------------------------
cd( [fileparts(which('setup__DMRIMatlab_toolbox')),'/mexcode'] );
warning('off','MATLAB:mex:GccVersion_link');
mid = 1;
% -----------------
modules(mid).name = 'dmri_2F1_';
modules(mid).src = './misfit';
modules(mid).depends = {'../mathsmex/hypergeom2F1.cxx','../threads/threadHelper.cpp'};
modules(mid).links = trlinks;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/maths/legendre'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'mexGenerateSHMatrix';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'mexTripleSHProdReal';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'mexAllSquaredSHFactors';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'mexTestPosODFsGrads';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../mathsmex/posODFsMaths.cxx'};
modules(mid).links = blaslinks;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'signal2sqrtshodf';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../mathsmex/posODFsMaths.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MiSFIT/utils'];
mid = mid+1;
% -----------------
modules(mid).name = 'sh2squaredsh_';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/sh'];
mid = mid+1;
% -----------------
modules(mid).name = 'shodf2samples';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = trlinks;
modules(mid).flags = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/sh'];
mid = mid+1;
% -----------------
modules(mid).name = 'sh2hot_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/sh2hothardcodes.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
mid = mid+1;
% -----------------
modules(mid).name = 'hot2sh_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/sh2hothardcodes.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
mid = mid+1;
% -----------------
modules(mid).name = 'hot2signal_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
mid = mid+1;
% -----------------
modules(mid).name = 'atti2hydidsi_';
modules(mid).src = './hydidsi';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../gcv/compute_gcv.cxx','../threads/threadHelper.cpp','../mathsmex/sanityCheckDTI.cxx','../quadprog/dmriquadprog.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/HYDI-DSI'];
mid = mid+1;
% -----------------
modules(mid).name = 'hydidsiQIntegrals_';
modules(mid).src = './hydidsi';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/HYDI-DSI'];
mid = mid+1;
% -----------------
modules(mid).name = 'atti2dti_';
modules(mid).src = './tensor';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tensor'];
mid = mid+1;
% -----------------
modules(mid).name = 'dti2spectrum_';
modules(mid).src = './tensor';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tensor'];
mid = mid+1;
% -----------------
modules(mid).name = 'atti2mapl_';
modules(mid).src = './mapl';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','./hermitePols.cxx','./maplMaths.cxx','../quadprog/dmriquadprog.cxx','../gcv/compute_gcv.cxx','../mathsmex/sanityCheckDTI.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MAPL'];
mid = mid+1;
% -----------------
modules(mid).name = 'mapl2atti_';
modules(mid).src = './mapl';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','./hermitePols.cxx','./maplMaths.cxx','../mathsmex/sanityCheckDTI.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MAPL'];
mid = mid+1;
% -----------------
modules(mid).name = 'mapl2index_';
modules(mid).src = './mapl';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','./hermitePols.cxx','./maplMaths.cxx','../mathsmex/sanityCheckDTI.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MAPL'];
mid = mid+1;
% -----------------
modules(mid).name = 'mapl2eap_';
modules(mid).src = './mapl';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','./hermitePols.cxx','./maplMaths.cxx','../mathsmex/sanityCheckDTI.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MAPL'];
mid = mid+1;
% -----------------
modules(mid).name = 'mapl2odf_';
modules(mid).src = './mapl';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','./hermitePols.cxx','./maplMaths.cxx','../mathsmex/sanityCheckDTI.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MAPL'];
mid = mid+1;
% -----------------
modules(mid).name = 'mexDMRIquadprog';
modules(mid).src = './quadprog';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','dmriquadprog.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
mid = mid+1;
% -----------------
modules(mid).name = 'mapl2samples';
modules(mid).src = './mapl';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','./hermitePols.cxx','./maplMaths.cxx','../mathsmex/sanityCheckDTI.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MAPL'];
mid = mid+1;
% -----------------

if(nargin>1)
    if( action=='x' )
        smodules = varargin(1:end);
        action = 'b';
    else
        smodules = varargin(2:end);
    end
elseif(nargin==1)
    if( action=='x' )
        smodules = varargin(1);
        action = 'b';
    else
        smodules = {modules(:).name}; % Select all
    end
else
    smodules = {modules(:).name}; % Select all
end

% Check if all requested modules correspond to elegible modules:
for n=1:length(smodules)
    if(   ~any( ismember({modules(:).name},smodules{n}) )   )
        error('<%s> doesn''t sound like a buildable module',smodules{n});
    end
end

switch(action)
    case 'l'
        names = {modules(:).name};
        fprintf(1,'These are the %d modules you can build:\n\n',mid-1);
        for n=1:length(names)
            fprintf(1,'   %s\n',names{n});
        end
        return;
    case 'b'
        for n=1:length(smodules)
            idx = find(ismember({modules(:).name},smodules{n}));
            build_my_mex_module(modules(idx)); %#ok<FNDSB>
        end
    case 'c'
        for n=1:length(smodules)
            idx = find(ismember({modules(:).name},smodules{n}));
            clean_my_mex_module(modules(idx)); %#ok<FNDSB>
        end
    case 'm'
        for n=1:length(smodules)
            idx = find(ismember({modules(:).name},smodules{n}));
            missing_my_mex_module(modules(idx)); %#ok<FNDSB>
        end
end

cd(oldpath);
end

function missing_my_mex_module(module)
totest = fullfile(module.dest,[module.name,'.',mexext]);
if(~exist(totest,'file'))
    fprintf(1,'Module <%s> in <%s> has not been built yet\n',[module.name,'.',mexext],module.dest);
end
end

function clean_my_mex_module(module)
todelete = fullfile(module.dest,[module.name,'.',mexext]);
if(exist(todelete,'file'))
    fprintf(1,'Removing <%s> from <%s>\n',[module.name,'.',mexext],module.dest);
    delete(todelete);
else
    fprintf(1,'Couldn''t find <%s> in <%s> [SKIP DELETION]\n',[module.name,'.',mexext],module.dest);
end
end

function build_my_mex_module(module)
cwd = pwd;
cd(module.src);
if(~ispc)
    mex( '-R2018a', ...
        module.flags{:}, ...
        [module.name,'.cpp'], ...
        module.depends{:}, ...
        module.links{:} ...
        );
else
    LPATH = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft');
    mex( '-R2018a', ...
        module.flags{:}, ...
        ['-L',LPATH], ...
        [module.name,'.cpp'], ...
        module.depends{:}, ...
        module.links{:} ...
        );
end
fprintf(1,'Moving <%s> to <%s>\n',[module.name,'.',mexext],module.dest);
movefile( [module.name,'.',mexext], module.dest );
cd(cwd);
end
