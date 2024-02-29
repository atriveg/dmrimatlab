% function [psi,nit,mu,Q,gradnorm,conderr] = ...
%                       signal2sqrtshodf(E,lambda,pshell,psi0,Ylm[,opts])
%
% Fits the squared root of a strictly non-negative ODF (i.e. an ODF that
% doesn't take negative values at any point of its continuous domain) with
% unit mass to the signal E in the SH basis according to a spherical 
% convolution model represented by convolution coefficients lambda and 
% shell indices pshell. The non-negative SH coefficients of the ODF itself 
% can be retrieved from the output psi of this function with a call to 
% sh2squaredsh.
%
% NOTE: if the squared root of the ODF expands up to an even order L, the 
% ODF itself will expand up to an even order 2L, hence L+1 convolution
% coefficients are required: 0, 2, 4, ..., 2L.
%
% NOTE: the unit mass constraint over the ODF translates in a unit norm
% constraint of the SH coefficients of its squared root.
%
% INPUTS:
%
%    E: N x FOV, each column is a vector of measurements within the FOV
%    lambda: (L+1) x P x FOV, the convolution factors. According to
%       Funk-Hecke's theorem, the SH coefficients of the ODF at each voxel,
%       phi(l,m), are related to the SH coeffcients of the signal E,
%       c(l,m), as:
%           c(l,m) = lambda(l,p)*phi(l,m),
%       where the convolution factor is assumed to be different at each of
%       the P subsets of dMRI measurements. In practice, this is often used
%       to describe P "shells", so that each of the N dMRI measurements
%       belongs to one and only one shell.
%    pshell: N x 1, a pointer to determine which of the P shells each dMRI
%       measurement comes from. NOTE: this pointer follows C's zero-based
%       indexing, meaning that the values of pshell MUST be in the range
%       [0,P-1].
%    psi0: K x FOV, with K = (L+1)(L+2)/2, this is an intial estimate of 
%       the output psi, since this is computed by means of numerical 
%       optimization using either Newton-Raphson's or Levenberg-Marquardt's
%       iterations.
%    Ylm: N x Kp, with Kp = (2L+1)(2L+2)/2, the SH matrix corresponding to
%       the evaluation of the SH basis (up to order 2L) at the orientations
%       of the N measurements.
%    opts: an optional structure with parameters to the algorithm:
%       opts.nu: the Laplace-Beltrami penalty over the non-negative ODF
%          s(default: 0.001);
%       opts.algorithm: either 'N' for Newton-Raphson's or 'L' for
%          Levenberg-Marquardt's . In the former case, the unit-mass 
%          constraint of the ODF is imposed via a Lagrange multiplier. In 
%          the former, the DC component of the squared root of the ODF is 
%          written as sqrt(1-||psi_nonDC||^2), and hence a pure least 
%          squares problem can be solved (default: 'L').
%       opts.T: the maximum number of iterations with either of the
%          previous methods (default: 100).
%       opts.thCost: with L-M's, the convergence threshold for the cost
%          function (default: 0.0001).
%       opts.thGrad: with L-M's or N-R's, the convergence threshold for the
%          gradient (default: 0.0001).
%       opts.thCond: with N-R's, the convergence threshold for the
%          unit-mass constraint of the ODF (default: 0.0001);
%       opts.rho0: the initial value of the adaptive damping factor for
%          both optimization methods (default: 1.0e-12).
%       opts.minrcn: minimum reciprocal condition number before a matrix is
%          considered singular and the damping factor is increased
%          (default: 1.0e-5).
%       opts.psi0eps: for L-M's, the minimum value allowed for psi(0,0),
%          the DC component of the squared root of the ODF (default:
%          1.0e-4).
%       opts.maxthreads: ONLY IN POSIX SYSTEMS the algorithm is run with
%          multiple threads. This is the maximum allowed number of threads,
%          which can indeed be reduced if it exceeds the number of logical
%          cores (default: the number of logical cores in the machine).
%
% OUTPUTS:
%
%    psi: K x FOV, the estimated SH coefficients.
%    nit: 1 x FOV, the number of iterations it took toconverge at each
%       voxel, or -1 if the algorithm failed to converge.
%    mu: 1 x FOV, for N-R's, the final value of the Lagrange multiplier
%       corresponding to the unit-mass constraint.
%    Q: 1 x FOV, the final value of the cost function (L-M's) or the
%       Lagrangian (N-R's).
%    gradnorm: 1 x FOV, the final value of the gradient of the cost
%       function (L-M's) or the Lagrangian (N-R's).
%    conderr: 1 x FOV, the final value of the error in the fulfillment of
%       the unit mass constraint (for N-R's).
%
% This is a mex function whose implementation you may find in file
% signal2sqrtshodf.c and the .h/.cxx files included therein.
%
function varargout = signal2sqrtshodf(varargin) %#ok<STOUT>
error('Please, build the mex code for this function by using the script in the ''mexcode'' folder');
end