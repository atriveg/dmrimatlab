% hydidsi2samples(eap,dti,lattice,Q,N,interp)
% hydidsi2samples(eap,dti,lattice,Q,N,interp,c)
% hydidsi2samples(eap,dti,lattice,Q,N,[],c)
% hydidsi2samples(eap,lattice,Q,N,*,c,seed)
% hydidsi2samples(eap,dti,lattice,Q,N,*,[],seed)
% hydidsi2samples(eap,dti,lattice,Q,N,*,*,seed,maxthreads)
% hydidsi2samples(eap,dti,lattice,Q,N,*,*,[],maxthreads)
% c = hydidsi2samples(___)
% [x,y,z] = hydidsi2samples(___)
% [x,y,z,c] = hydidsi2samples(___)
%
% Generates N random samples (i.e. N random 3-D displacements)
% distributed according to the Ensemble Average Propagator (EAP)
% whose HYDI-DSI coefficients are provided by the 1xK array eap
% (with grid parameters lattice, Qx, Qy and Qz, all of them
% returned by atti2hydidsi), using the rejection method as
% described by:
%
%      Luc Devroye. "Non-uniform random variate generation".
%      Springer Science+Business Media, LLC. 1986.
%
%   INPUTS:
%
%      - eap: A 1 x K vector of doubles with the coefficients of the EAP,
%            where K=[(2*Nx+1)*(2*Ny+1)*(2*Nz+1)+1]/2.
%      - dti: A 1 x 6 vector of doubles with the (unique) coefficients of
%            the diffusion tensor computed at the voxel to be simulated,
%            which is used to stretch and orient the space. This is also
%            returned by the atti2hydidsi function.
%      - lattice: a 1 x 3 vector, lattice = [Nx,Ny,Nz], describing the
%            sampling of the EAP in the cartesian grid. Note this is either
%            fed to, or returned by, atti2hydidsi, so it must be consistent
%            with the size of the eap array.
%      - Q: a 1 x 3 vector, Q = [Qx,Qy,Qz] with the bandwidths returned by
%            atti2hydidsi at each voxel.
%      - N: A 1 x 1 integer, the number of random samples to generate.
%      - interp: a 1 x 1 char, either 'n', 'l' or 's', telling the algorihtm
%            if the EAP at continuous points should be interpolated by
%            nearest neighbor, linear or cubic B-spline interpolation (leave
%            empty for dedault 'l').
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
%   passing the same HYDI-DSI coeffs. Yet, externally providing a seed can
%   improve the computational performance. You can check the
%   <test_hydidsi2samples.m> script for a simple demonstration.
function varargout = hydidsi2samples(eap,dti,lattice,Q,N,intrp,c,seed,maxthreads)
if(nargin<4)
    error('At least the eap, the dti, the lattice and the bandwidth factors must be supplied');
end
% ----------------
assert( isequal(size(dti),[1,6]), 'The dti coeffs must be 1x6' );
assert( isequal(size(lattice),[1,3]), 'The lattice must be 1x3' );
assert( isequal(lattice,round(lattice)), 'The lattice must contain only integer numbers' );
assert( isequal(size(Q),[1,3]), 'The bandwidths vector must be 1x3' );
K = (prod(2*lattice+1)+1)/2;
assert( isequal(size(eap),[1,K]), sprintf('The eap coeffs must be 1 x %d',K) );
% ----------------
switch(nargout)
    case 0
        % Nothing to do
        return;
    case 1
        % Just compute an estime of c
        varargout{1} = hydidsi2samples_compute_c(eap,lattice,Q);
        return;
    otherwise
        % Actually generate samples
        assert(nargin>=5,'At least the eap, the dti, the lattice, the bandwidth factors and N must be supplied');
        assert(isscalar(N),'The number of samples, N, must be a scalar');
        assert(N==round(N),'The number of samples, N, must be integer');
        if( nargin<6 || isempty(intrp) )
            intrp = 'l';
        else
            assert(isscalar(intrp),'intrp must be a single char');
            assert(ischar(intrp),'intrp must be a single char');
            assert( any(intrp=='nls'), 'intrp must be either ''n'', ''l'' or ''s''');
        end
        sigma = [];
        if( nargin<7 || isempty(c) ) % Need to estimate c
            [c,sigma] = hydidsi2samples_compute_c(eap,lattice,Q);
        else
            assert(isscalar(c),'c must be a scalar');
            assert(c>1,'c must be greater than 1');
        end
        if( nargin<8 || isempty(seed) ) % Need to create a seed
            seed = datetime('now').Second;
            seed = seed - floor(seed);
            seed = 1e6*seed; % microseconds
            seed = [ seed, seed*7907, seed*6991 ];
        else
            assert( isequal(size(seed),[1,3]), 'The seed must be a 1x3 vector' );
            seed = round(seed);
        end
        if( nargin<9 || isempty(maxthreads) )
            maxthreads = 1e6; % Internally adjusted by the mex function
        else
            assert(isscalar(maxthreads),'maxthreads must be a scalar');
            maxthreads = round(maxthreads);
        end
        % Get the variance of the Gaussian envelope, if necessary:
        if(isempty(sigma))
            [~,sigma] = hydidsi2samples_compute_c(eap,lattice,Q);
        end
        % Call the mex function:
        [x,y,z] = hydidsi2samples_( double(eap)/prod(Q), double(lattice), double(c), ...
            double(sigma), double(N), char(intrp), double(seed), double(maxthreads) );
        % The mex function creates normalized samples in [-Ni,Ni], but we
        % want them to be in [-Ni/Qi,Ni/Qi]. Hence:
        x = x/Q(1);
        y = y/Q(2);
        z = z/Q(3);
        % Yet, the lattice refers to a rotated space where z is aligned
        % with the main eigenvector, so that we need to rotate these
        % samples to undo this transform:
        [u1,u2,u3] = dti2spectrum( reshape(dti,[1,1,1,6]) ); % In descending order
        U  = [u3(:),u2(:),u1(:)]; % We need them in ascending order!!
        U(:,1) = cross(U(:,2),U(:,3)); % Make sure U is a rotation matrix
        samples = [x,y,z]*(U');
        % Return the samples:
        for n=1:min(nargout,3)
            varargout{n} = samples(:,n); %#ok<AGROW>
        end
        % If necessary, return c:
        if(nargout>3)
            varargout{4} = c;
        end
end

end

% -------------------------------------------------------------------------
function [c,sigma] = hydidsi2samples_compute_c(eap,lattice,Q)
% Compute the {x,y,z} normalized values where the eap has
% been evaluated:
x = -lattice(1):lattice(1);
y = -lattice(2):lattice(2);
z = -lattice(3):lattice(3);
[x,y,z] = meshgrid(x,y,z);
x = x( (numel(x)+1)/2 : end );
y = y( (numel(y)+1)/2 : end );
z = z( (numel(z)+1)/2 : end );
% Normalize the eap to get rid of the dependency with Q:
eap = eap./prod(Q);
% The value of the normalized Gaussian will be nearly:
sigma0  = sqrt(-0.5/log(0.05));
% But we will probe values near this one:
factors = exp(   linspace( log(0.9), log(2), 21 )   );
cs      = zeros(size(factors));
% For each of these factors:
Nx2    = lattice(1)*lattice(1);
Ny2    = lattice(2)*lattice(2);
Nz2    = lattice(3)*lattice(3);
Gconst = sqrt(8*pi*pi*pi*Nx2*Ny2*Nz2);
for f = 1:length(factors)
    % Compute the Gaussian kernel:
    sigma  = sigma0*factors(f);
    sigma2 = sigma*sigma;
    sigma3 = sigma2*sigma;
    G      = exp(-x.*x/(2*sigma2*Nx2)).*exp(-y.*y/(2*sigma2*Ny2)).*...
        exp(-z.*z/(2*sigma2*Nz2))/Gconst/sigma3;
    % Compute the value of c as the constant that guarantees
    % c*G >= eap for all the nodes in the lattice
    cs(f)  = max(eap./G);
end
% Keep the best combination, which means the minimum value
% of c:
[~,pos] = min(cs);
c       = cs(pos);
sigma   = sigma0*factors(pos);
end
