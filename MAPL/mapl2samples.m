% mapl2samples(coeffs,dti,tau,N)
% mapl2samples(coeffs,dti,tau,N,c)
% mapl2samples(coeffs,dti,tau,N,c,seed)
% mapl2samples(coeffs,dti,tau,N,[],seed)
% mapl2samples(coeffs,dti,tau,N,*,seed,maxthreads)
% mapl2samples(coeffs,dti,tau,N,*,[],maxthreads)
% c = mapl2samples(___)
% [x,y,z] = mapl2samples(___)
% [x,y,z,c] = mapl2samples(___)
%
% Generates N random samples (i.e. N random 3-D displacements)
% distributed according to the Ensemble Average Propagator (EAP)
% whose MAP-MRI coefficients are provided by the 1xK array coeffs,
% using the rejection method as described by:
%
%      Luc Devroye. "Non-uniform random variate generation". 
%      Springer Science+Business Media, LLC. 1986.
%
%   INPUTS:
%
%      - coeffs: A 1 x K vector of doubles with the coefficients of the ODF,
%            where K=(Nmax+2)*(Nmax+4)*(2*Nmax+3)/24 for some even
%            integer Nmax>=0. These are the MAPL cofficients as estimated
%            with atti2mapl.
%      - dti: A 1 x 6 vector of doubles with the (unique) coefficients of
%            the diffusion tensor computed at the voxel to be simulated,
%            which is used to stretch and orient the space. This is also
%            returned by the atti2mapl function. NOTE: You can generate
%            DTI-based samples by simply passing 1 (i.e. a 1x1 matrix whose
%            only entry is 1) as the first input (coeffs).
%      - tau: the effective diffusion time you want to simulate in seconds.
%      - N: A 1 x 1 integer, the number of random samples to generate.
%      - c (optional): the constant used by the rejection method. If it is
%            ommitted or an empty array [] is passed, the function will
%            internally compute it (and can be returned as an output). It
%            is recommended that you let the function decide.
%      - seed (optional): a 1x3 vector used to initialize the state of the
%            random number generator (the C++ drand48_r() routine is used,
%            so that these 3 elements are at the end of the day combined
%            into a 48 bits integer). If ommitted (or empty), the function
%            itself will choose a proper initial state based on the CPU 
%            time.
%      - maxthreads (optional): a 1x1 integer with the maximum number of
%            threads used. If ommitted, as many threads as available CPUs
%            are used. If a number greater than the number of CPUs is
%            passed, it is cropped to the number of CPUs.
%
%   OUTPUTS:
%
%      - c: A 1 x 1 double with the constant finally used by the rejection
%            method.
%      - x, y, z: N x 1 each, three vectors with doubles so that the vector
%            [x(n),y(n),z(n)]^T, n = 1..N is the n-th random displacement
%            generated.
%
%   NOTE: This is implemented as a multi-threaded mex function. For the
%   fastest performance, please consider calling the function just once
%   with a large enough N instead of repeteadly calling the function
%   passing the same MAPL coeffs. Yet, externally providing a seed can
%   improve the computational performance. You can check the
%   <test_mapl2samples.m> script for a simple demonstration.
function varargout = mapl2samples(varargin) %#ok<STOUT>
error('Please, build the mex code for this function by using the script in the ''mexcode'' folder');
end
