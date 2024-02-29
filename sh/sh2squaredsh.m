% function shodf = sh2squaredsh( shsqrt[, mask] )
%
%   Computes the SH coefficients of a non-negative ODF, shodf, from the SH
%   coefficients of its squared root, shsqrt, using the analytical
%   relationships between them provided by Wigner's symbols.
%
%   NOTE: if the squared root of the ODF spans up to a maximum order L
%   within the SH basis, the ODF itself spans up to an order 2L.
%
%   NOTE: if the ODF needs to be unit-mass, then the SH coefficients of its
%   squared root must be a unit-norm vector.
%
%   INPUTS:
%
%   shsqrt: N1 x N2 x N3 x ... x ND x K, with K=(L+1)(L+2)/2 for some
%      non-negative, even integer L. The SH coefficients of the squared
%      root of the ODF
%   mask: N1 x N2 x N3 x ... x ND, a mask to avoid unnecesary computations.
%      Zero valued mask positions are filled with zeros in the output
%      array.
%
%   OUTPUTS:
%
%   shodf: N1 x N2 x N3 x ... x ND x Kp, with Kp=(2L+1)2(L+2)/2. The SH
%      coefficients of the ODF.
%
%   This is a mex function implemented in sh2squaredsh.c and the .h/.cxx
%   files included therein.
function varargout = sh2squaredsh(varargin) %#ok<STOUT>
error('Please, build the mex code for this function by using the script in the ''mexcode'' folder');
end