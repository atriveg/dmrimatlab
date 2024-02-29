% - Script to build mex functions
% - NOTE: if using windows, either MinGW or Cywin/gcc are strongly
%         recommended.
% - NOTE: for Mac, we have disabled lpthreads until we find a safe way to
%         prevent BLAS/LAPACK from using their own inner threads. This is
%         actually implemented within "mexToMathsTypes.h"
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
modules(mid).others = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'mexTripleSHProdReal';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).others = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'mexAllSquaredSHFactors';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
modules(mid).links = {};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).others = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'mexTestPosODFsGrads';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../mathsmex/posODFsMaths.cxx'};
modules(mid).links = {'-lmwlapack','-lmwblas'};
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tests'];
modules(mid).others = {};
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'signal2sqrtshodf';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx','../mathsmex/posODFsMaths.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/MiSFIT/utils'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'sh2squaredsh';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/matrixCalculus.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/sh'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'shodf2samples';
modules(mid).src = './sh';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread'};
    modules(mid).others = {};
else
    modules(mid).links = {};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/sh'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'sh2hot_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/sh2hothardcodes.cxx','../mathsmex/matrixCalculus.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'hot2sh_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/sh2hothardcodes.cxx','../mathsmex/matrixCalculus.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'hot2signal_';
modules(mid).src = './hot';
modules(mid).depends = {'../mathsmex/sphericalHarmonics.cxx','../mathsmex/sh2hot.cxx','../mathsmex/matrixCalculus.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/hot'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'atti2hydidsi_';
modules(mid).src = './hydidsi';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx','../gcv/compute_gcv.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/HYDI-DSI'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'hydidsiQIntegrals_';
modules(mid).src = './hydidsi';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/HYDI-DSI'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'atti2dti_';
modules(mid).src = './tensor';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
modules(mid).dest = [fileparts(which('setup__DMRIMatlab_toolbox')),'/tensor'];
modules(mid).ignore = false;
mid = mid+1;
% -----------------
modules(mid).name = 'dti2spectrum_';
modules(mid).src = './tensor';
modules(mid).depends = {'../mathsmex/matrixCalculus.cxx'};
if( ~ispc && ~ismac )
    modules(mid).links = {'-lpthread','-lmwlapack','-lmwblas'};
    % This is necessary for calling OpenMP thread managing functions:
    modules(mid).others = {'CXXFLAGS="$CXXFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"'};
else
    modules(mid).links = {'-lmwlapack','-lmwblas'};
    modules(mid).others = {};
end
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
            module.others{:}, ...
            [module.name,'.cpp'], ...
            module.depends{:}, ...
            module.links{:} ...
            );
    else
        LPATH = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft');
        mex( '-R2018a', ...
            module.others{:}, ...
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
