function [atti,gi,bi] = dwi2atti(dwi,gi0,bi0,varargin)
% function [atti,gi,bi] = dwi2atti(dwi, gi0, bi0, 'option1', value1, ...)
%
%   Computes the attenuation signal atti from a set of diffusion weighted
%   images dwi:
%
%      dwi: a MxNxPx(G+B) double array comprising both the DWI channels and
%           the baseline or baselines (which will be averaged in case they
%           are more than one, B>1)
%      gi0: a (G+B)x3 matrix with the set of input gradients, including
%           those corresponding to the baselines. In case the module
%           of the n-th row is different from 1, the n-th b-value,
%           bi0(n), will be accordingly scaled.
%      bi0: a (G+B)x1 vector of b-values, including the zero-values
%           corresponding to the baselines. In case a scalar is passed, it
%           is assumed the same b-value is used for the whole data set.
%
%      atti: a MxNxPxG double array with the Si/S0 attenuation signal at
%           each voxel within the MxNxP image frame and for each i=1...G
%           gradient direction. In case B>1 (more than one baseline image
%           is present in the data set), all baselines are averaged.
%      gi:  a Gx3 matrix with the gradient directions exluding those
%           corresponding to baselines
%      bi:  a Gx1 vector with the b-values for each gradient direction,
%           excluding those corresponding to baselines
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      b0th: b-values below this threshold will be considered to be 
%           baseline images (default: 1). Units are s/mm^2.

% -------------------------------------------------------------------------
if(nargin<3)
    error('At least dwi, gi0, and bi0 must be supplied');
end
% -------------------------------------------------------------------------
[~,~,~,G] = size(dwi);
if(numel(bi0)==1)
    bi0 = bi0*ones(G,1);
end
if( size(gi0,1)~=G || size(bi0,1)~=G )
    error('The size of the gradients table and b-values must match the 4-th dimension of the DWI');
end
if( size(gi0,2)~=3 )
    error('gi0 should have size (G+B)x3');
end
if( size(bi0,2)~=1 )
    error('bi0 should have size (G+B)x1');
end
% -------------------------------------------------------------------------
opt.b0th = 1;      optchk.b0th = [true,true];    % 1x1 double
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------

modules = sqrt(sum(gi0.*gi0,2)); % Gx1
bsidx   = (abs(bi0)<opt.b0th);   % baselines
gridx   = ~bsidx;                % gradients

baseline = mean( dwi(:,:,:,bsidx), 4 ); % MxNxP
mask     = (baseline>1.0e-6);
baseline(~mask) = inf;

atti     = dwi(:,:,:,gridx);  % MxNxPxG
modules  = modules(gridx);    % Gx1
gi       = gi0(gridx,:);      % Gx3
bi       = bi0(gridx).*modules;

if(is_broadcast_available)
    atti = atti./baseline;
    gi   = gi./modules;
else
    atti = bsxfun( @(x,y)(x./y), atti, baseline );
    gi =   bsxfun( @(x,y)(x./y), gi, modules );
end
