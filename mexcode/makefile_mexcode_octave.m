function makefile_mexcode_octave(varargin)
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
%
% ADDITIONAL NOTES FOR OCTAVE (Linux-based):
%
%   Please, refer to the README.octave file within the same subfolder
%   this script is located at ('mexcode').
%

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
path0       = mfilename('fullpath');
[path0,~,~] = fileparts(path0);

BLAS_MODE = get_BLAS_mode(path0);

if(isunix)
    setenv( 'CPPFLAGS', '-DOCTAVE_BUILD' );
    setenv( 'CXXFLAGS', '-DMX_COMPAT_64  -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -pthread -fwrapv -O3 -DNDEBUG' );
    trlinks   = {'-lpthread'};
    switch(BLAS_MODE)
        case 1,
            % Netlib's BLAS, needs cblas and lapacke packages
            blaslinks = {'-lcblas','-llapack'};
            blasflags = {'-D_SYSTEM_BLAS_BUILD_'};
        case 2,
            % System-wide OpenBLAS, needs openblas package
            blaslinks = {'-lopenblas'};
            blasflags = {'-D_SYSTEM_OPENBLAS_BUILD_','-D_USE_OPENBLAS_THREAD_CONTROL'};
        case 3,
            % Locally compiled OpenBLAS, single-thread with locks
            suffix = check_local_openblas_available(path0,'single-thread');
            blaslinks = { sprintf('-lopenblas_%s',suffix), sprintf('-L%s/openblas-%s/lib',path0,suffix), sprintf('-Wl,-rpath=%s/openblas-%s/lib',path0,suffix) };
            blasflags = {'-D_LOCAL_OPENBLAS_BUILD_', sprintf('-I %s/openblas-%s/include',path0,suffix) };
        case 4,
            % System-wide OpenBLAS, avoid direct calls to BLAS functions (very inefficient)
            blaslinks = {'-lopenblas'};
            blasflags = {'-D_SYSTEM_OPENBLAS_BUILD_','-D_NO_BLAS_CALLS'};
        case 5,
            % Use Intel MKL
            mklroot = '/opt/intel/oneapi/mkl/2025.0';
            others = '/opt/intel/oneapi/redist/lib';
            blaslinks = { sprintf('-L%s/lib',mklroot), sprintf('-L%s',others), sprintf('-Wl,-rpath=%s/lib',mklroot), sprintf('-Wl,-rpath=%s',others), ...
                '-lmkl_intel_lp64', '-lmkl_intel_thread', '-lmkl_core', '-liomp5', '-lpthread', '-lm', '-ldl' };
            blasflags = {'-D_MKL_BLAS_BUILD_','-D_USE_MKL_THREAD_CONTROL', '-m64', '-Wl,--no-as-needed', sprintf('-I"%s/include"',mklroot),  };
        otherwise,
            error('Unable to determine the BLAS implementation to use');
    end
end

if(ispc)
    error('mexcode not implemented yet for Windows'' Octave');
end

if( ismac )
    error('mexcode not implemented yet for MAC OSX''s Octave');
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
modules(mid).flags = blasflags;
mid = mid+1;
% -----------------
modules(mid).name = 'mexGenerateSHMatrix';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests/bin'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'mexTripleSHProdReal';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests/bin'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'mexAllSquaredSHFactors';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests/bin'];
modules(mid).flags = {};
mid = mid+1;
% -----------------
modules(mid).name = 'mexTestPosODFsGrads';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../mathsmex/posODFsMaths.cxx'};
modules(mid).links = blaslinks;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests/bin'];
modules(mid).flags = blasflags;
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
modules(mid).flags = blasflags;
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
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests/bin'];
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
modules(mid).name = 'fastSweepingStep';
modules(mid).src = './finslertractography';
modules(mid).depends = {'../pixelIterators/iterators.cxx'};
modules(mid).links  = {};
modules(mid).flags = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/finslertractography'];
mid = mid+1;
% -----------------
modules(mid).name = 'backTracing_';
modules(mid).src = './finslertractography';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','../pixelIterators/iterators.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/finslertractography'];
mid = mid+1;
% -----------------
modules(mid).name = 'fiberTracking_';
modules(mid).src = './dtiTractography';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp','../pixelIterators/iterators.cxx'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = blasflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tractography'];
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
    baseflags = getenv('CXXFLAGS');
    newflags  = '';
    for k=1:length(module.flags)
        newflags = sprintf('%s %s',newflags,module.flags{k});
    end
    newflags = sprintf('%s %s',newflags,baseflags);
    setenv( 'CXXFLAGS', newflags );
    mex( [module.name,'.cpp'], ...
        module.depends{:}, ...
        module.links{:} ...
        );
    setenv( 'CXXFLAGS', sprintf('%s ',baseflags) );
else
    error('Window''s Octave is not able to compile our mex files by the moment');
end
fprintf(1,'Moving <%s> to <%s>\n',[module.name,'.',mexext],module.dest);
movefile( [module.name,'.',mexext], module.dest );
cd(cwd);
clear(module.name);
end


function suffix = check_local_openblas_available(path0,opts)
suffix    = get_BLAS_suffix(path0);
available = true;
available = available && (   exist( sprintf('%s/openblas-%s/lib/libopenblas_%s.so',path0,suffix,suffix), 'file' )   ~=   0   );
available = available && (   exist( sprintf('%s/openblas-%s/include/cblas.h',path0,suffix), 'file' )   ~=   0   );
available = available && (   exist( sprintf('%s/openblas-%s/include/lapacke.h',path0,suffix), 'file' )   ~=   0   );

if(~available)
    clc;
    fprintf(1,'The build option you have chosen via config.octave means\n');
    fprintf(1,'that the mex files will be linked against a local, non-\n');
    fprintf(1,'threaded version of OpenBLAS. However, this local version is\n');
    fprintf(1,'not available. I will try now to download, compile, and locally\n');
    fprintf(1,'install it for you. This process should be transparent for you,\n');
    fprintf(1,'but it might take a while. Please make sure you have the\n');
    fprintf(1,'following software installed:\n');
    fprintf(1,'   - git \n');
    fprintf(1,'   - gcc and g++ \n');
    fprintf(1,'   - gfortan or alike\n');
    fprintf(1,'In case this fails, please set some other value of BLAS_MODE\n');
    fprintf(1,'and try again (help makefile_mexcode_octave for details)\n\n');
    fprintf(1,'[PRESS ENTER to go ahead]\n');
    pause;
    status = system ( sprintf('bash %s/build_local_openblas.sh %s %s',path0,opts,suffix) );
    if(status==0)
        fprintf('\nSUCCEEDED!!!\n');
    else
        fprintf('\nFAILED\n');
    end
    fprintf(1,'[PRESS ENTER to go ahead]\n');
    pause;
end

end

function BLAS_MODE = get_BLAS_mode(path0)
lines = get_BLAS_config(path0);
if(length(lines)<1)
    BLAS_MODE = 3;
    return;
end
switch(lower(lines{1}))
    case 'netlib'
        BLAS_MODE = 1;
    case 'openblas'
        BLAS_MODE = 2;
    case 'openblas-local'
        BLAS_MODE = 3;
    case 'mkl'
        BLAS_MODE = 5;
    otherwise
        error(sprintf('Unable to parse option <%s> in config.octave. Choose one of [ netlib | openblas | openblas-local ] or delete the config file to use defaults',lines{1}));
end
end

function suffix = get_BLAS_suffix(path0)
lines = get_BLAS_config(path0);
if(length(lines)<2)
    suffix = 'local';
else
    suffix = lines{2};
end
end

function lines = get_BLAS_config(path0)
config = sprintf('%s/config.octave',path0);
if( exist(config,'file')~=2 )
    warning(sprintf('Unable to open config file ''%s'' for reading. Creating it now with default options',config));
    fid = fopen(config,'w');
    fprintf(fid,'openblas-local\n');
    fprintf(fid,'local\n');
    fclose(fid);
end
fid = fopen(config,'r');
str = fgetl(fid);
lns = 0;
while(str~=-1)
    lns = lns+1;
    lines{lns} = str;
    str = fgetl(fid);
end
fclose(fid);
end





















