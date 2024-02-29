% shodf2samples(SH,N)
% shodf2samples(SH,N,c)
% shodf2samples(SH,N,c,seed)
% shodf2samples(SH,N,[],seed)
% c = shodf2samples(___)
% [theta,phi] = shodf2samples(___)
% [theta,phi,c] = shodf2samples(___)
%
% Generates N random samples (i.e. N random spatial orientations)
% distributed according to the Oriention Distribution Function (ODF) 
% whose SH coefficients are provided by the 1xK array SH, using the
% rejection method as described by:
%
%      Luc Devroye. "Non-uniform random variate generation". 
%      Springer Science+Business Media, LLC. 1986.
%
%   INPUTS:
%
%      - SH: A 1 X K vector of doubles with the coefficients of the ODF,
%            where K=(L+1)(L+2)/2 for some even integer L>=0.
%      - N: A 1 x 1 integer, the number of random samples to generate.
%      - c (optional): the constant used by the rejection method. If it is
%            ommitted or an empty array [] is passed, the function will
%            internally compute it (and can be returned as an output).
%      - seed (optional): a 1x3 vector used to initialize the state of the
%            random number generator (the C++ drand48_r() routine is used,
%            so that these 3 elements are at the end of the day combined
%            into a 48 bits integer). If ommitted, the function itself will
%            choose a proper initial state based on the CPU time.
%
%   OUTPUTS:
%
%      - c: A 1 x 1 double with the constant finally used by the rejection
%            method.
%      - theta, phi: Nx1 each, two vectors with doubles within the ranges
%            0 <= theta < pi and 0 <= phi < 2*pi that represent spatial
%            orientations in spherical coordinates:
%                  x = sin(theta)*cos(phi);
%                  y = sin(theta)*sin(phi);
%                  z = cos(theta);
%            These orientations are not evenly distributed within the unit
%            sphere, but they follow instead the probability law described
%            by the ODF.
%
%   NOTE: This is implemented as a multi-threaded mex function. For the
%   fastest performance, please consider calling the function just once
%   with a large enough N instead of repeteadly calling the function
%   passing the same ODF. Yet, externally providing a seed can improve the
%   computational performance. You can check the <test_shodf2samples.m>
%   script for a simple demonstration.
function varargout = shodf2samples(varargin) %#ok<STOUT>
error('Please, build the mex code for this function by using the script in the ''mexcode'' folder');
end