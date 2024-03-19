% - Script to build mex functions
% - NOTE: if using windows, either MinGW or Cywin/gcc are strongly
%         recommended.

% ------------------------------------------------------------------------------
mtlroot = matlabroot;
if(isunix)
    libsdir = sprintf('%s/bin/glnxa64',mtlroot);
    mkllib  = sprintf('%s/mkl.so',libsdir);
    blaslinks = {'-lmwlapack','-lmwblas'};
    trlinks = {'-lpthread',mkllib};
    trflags = {'CXXFLAGS="$CXXFLAGS -D_USE_MKL_THREAD_CONTROL"'};
end

if(ispc)
    blaslinks = {'-lmwlapack','-lmwblas'};
    trlinks = {};
    trflags = {'CXXFLAGS="$CXXFLAGS -fopenmp -D_USE_OMP_THREAD_CONTROL"','LDFLAGS="$LDFLAGS -fopenmp"'};
end

if( ismac )
    blaslinks = {'-lmwlapack','-lmwblas'};
    % Check if mkl.dylib is present. If it is the case, this is an old
    % Mac machine with an Intel processor, hence we can rely on MKL as
    % we do for linux (but we must explicitly tell the computer through
    % the -D_HAS_MKL_BLAS flag).
    libsdir = sprintf('%s/bin/maci64',mtlroot);
    mkllib  = sprintf('%s/mkl.dylib',libsdir);
    if( exist(mkllib,'file') ~= 0 )
        % MAC with intel processor
        trlinks = {'-lpthread',mkllib};
        trflags = {'CXXFLAGS="$CXXFLAGS -D_USE_MKL_THREAD_CONTROL"'};
    else
        % MAC with ARM processor
        trlinks = {'-lpthread','-lmwopenblas'};
        trflags = {'CXXFLAGS="$CXXFLAGS -D_USE_OPENBLAS_THREAD_CONTROL"'};
    end
end
% ------------------------------------------------------------------------------
oldpath = pwd;
% ------------------------------------------------------------------------------
cd( [fileparts(which('setup__DMRIMatlab_toolbox')),'/mexcode'] );
warning('off','MATLAB:mex:GccVersion_link');
mid = 1;
% -----------------
modules(mid).name = 'mexGenerateSHMatrix';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'mexTripleSHProdReal';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'mexAllSquaredSHFactors';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'mexTestPosODFsGrads';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../mathsmex/posODFsMaths.cxx'};
modules(mid).links = blaslinks;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).flags = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'signal2sqrtshodf';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../mathsmex/posODFsMaths.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MiSFIT/utils'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'sh2squaredsh';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/sh'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'shodf2samples';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = trlinks;
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/sh'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'sh2hot_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/sh2hothardcodes.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'hot2sh_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/sh2hothardcodes.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'hot2signal_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'atti2hydidsi_';
modules(mid).src = './hydidsi';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../gcv/compute_gcv.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/HYDI-DSI'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'hydidsiQIntegrals_';
modules(mid).src = './hydidsi';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/HYDI-DSI'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'atti2dti_';
modules(mid).src = './tensor';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tensor'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'dti2spectrum_';
modules(mid).src = './tensor';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../threads/threadHelper.cpp'};
modules(mid).links  = [ blaslinks, trlinks ];
modules(mid).flags = trflags;
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tensor'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------

for m=1:mid-1
    clean_my_mex_module(modules(m));
end

for m=1:mid-1
    build_my_mex_module(modules(m));
end

cd(oldpath);

function clean_my_mex_module(module)
if(~module.ignore)
    todelete = fullfile(module.dest,[module.name,'.',mexext]);
    if(exist(todelete,'file'))
        fprintf(1,'Removing <%s> from <%s>\n',[module.name,'.',mexext],module.dest);
        delete(todelete);
    else
        fprintf(1,'Couldn''t find <%s> in <%s> [SKIP DELETION]\n',[module.name,'.',mexext],module.dest);
    end
end
end

function build_my_mex_module(module)
if(~module.ignore)
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
end
