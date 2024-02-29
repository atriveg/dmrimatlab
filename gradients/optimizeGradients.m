function gi = optimizeGradients(gi,varargin)
% function gi = optimizeGradients(gi,'opt1',value1,'opt2',value2)
%
%    Optimizes a Gx3 set of gradients by maximizing the overall distance
%    from one gradient to the others. An energy function defined as:
%
%         Q = 1/2 * sum_{i\not=j} 1/sin^p(alpha_ij) [1]
%
%    is minimized to that end, being alpha_ij the angle between g(i,:) and
%    g(j,:).
%
%       gi: Gx3 matrix of normalized gradients (just in the north
%           hemisphere, since it is assumed that gi(i,:) implies -gi(i,:))
%
%    Options:
%
%       p: 1x1 double, it defines the energy term in [1] (default: 1)
%       iters: 1x1 double, maximum number of iterations for the
%          optimization (default: 100)
%       lambda: 1x1 double, Levenberg-Marquardt regularization term in the 
%          Newton-Raphson step. This is just the initial value, it is
%          adaptively updated as the iterations proceed (default: 1.0)
%       delta: maximum allowed change in the parameters before the
%          iterations stop (default: 1.0e-5)
%       exclude: a vector of unique indices within the range [1,G] that are
%          excluded from the optimization, i.e., they will remain unchanged
%          throughout the entire process. By default, exclude=1, so that
%          the first gradient is anchored and the others are moved
%          (default: [1])
%       verbose: 1x1 logical, print additional information (default: false)

assert(size(gi,2)==3,'gi must be Gx3');
G = size(gi,1);
% Parse the optional input arguments:
opt.p = 1;           optchk.p = [true,true];        % always 1x1 double
opt.iters = 100;     optchk.iters = [true,true];    % always 1x1 double
opt.lambda = 1.0;    optchk.lambda = [true,true];   % always 1x1 double
opt.delta = 1.0e-5;  optchk.delta = [true,true];    % always 1x1 double
opt.exclude = 1;     optchk.exclude = [true,false]; % vector of integers with free size
opt.verbose = false; optchk.verbose = [true,true];  % always 1x1 boolean
opt = custom_parse_inputs(opt,optchk,varargin{:});

xi = gi(:,1);
yi = gi(:,2);
zi = gi(:,3);
ti = atan2( sqrt(xi.*xi+yi.*yi), zi );
pi = atan2( yi, xi );

if(opt.verbose)
    txtlength = 0;
end
idxs = setdiff((1:G)',opt.exclude(:));
idxl = [idxs;idxs+G];
lmb  = opt.lambda;
for it=1:opt.iters
    % Compute tigonometric functions of theta and phi
    cti = cos(ti);
    sti = sin(ti);
    cpi = cos(pi);
    spi = sin(pi);
    % Compute the angles between each pair of gradients:
    [Caij,Saij] = compute_solid_angles_from_spherical(cti,sti,cpi,spi);
    % Compute the overall cost:
    Q  = compute_energy_from_angles(Saij,opt.p);
    % Compute the gradient of the cost:
    DQ  = compute_energy_gradient_from_angles(cti,sti,cpi,spi,Caij,Saij,opt.p,false);
    % DEBUG: DQt = compute_energy_gradient_from_angles(cti,sti,cpi,spi,Caij,Saij,opt.p,true);
    % Compute the Hessian of the cost:
    DDQ  = compute_energy_hessian_from_angles(cti,sti,cpi,spi,Caij,Saij,opt.p,false);
    % DEBUG: DDQt = compute_energy_hessian_from_angles(cti,sti,cpi,spi,Caij,Saij,opt.p,true);
    % (Regularized) Newton-Raphson step. One gradient direction will remain
    % fixed, and the remaining ones will update:
    LSM     = DDQ(idxl,idxl);
    LSMnorm = sqrt( trace(LSM*LSM')/(size(LSM,1)) );
    LSM     = LSM + lmb*LSMnorm*eye(size(LSM,1));
    LSV     = DQ(idxl);
    parc    = [ti(idxs);pi(idxs)];
    parn    = parc - LSM\LSV;
    % Check if the cost function actually decreases:
    tin         = ti;
    pin         = pi;
    tin(idxs)   = parn( 1:length(parn)/2 );
    pin(idxs)   = parn( (1:length(parn)/2) + length(parn)/2 );
    cti         = cos(tin);
    sti         = sin(tin);
    cpi         = cos(pin);
    spi         = sin(pin);
    [~,Saij] = compute_solid_angles_from_spherical(cti,sti,cpi,spi);
    Qn          = compute_energy_from_angles(Saij,opt.p);
    if(Qn<=Q)
        ti    = tin;
        pi    = pin;
        lmb   = lmb/2;
        delta = max(abs(parn-parc));
        if(opt.verbose)
            fprintf(1,repmat('\b',[1,txtlength]));
            msg = sprintf('Iter: %d; Q: %1.5f -> %1.5f; delta: %1.4f; lambda=%1.5f',it,Q,Qn,delta,lmb);
            fprintf(1,msg);
            txtlength = length(msg);
        end
        if(delta<opt.delta)
            break;
        end
    else
        lmb = lmb*2;
    end
end

if(opt.verbose)
    fprintf(1,repmat('\b',[1,txtlength]));
end
    
gi = [ sin(ti).*cos(pi), sin(ti).*sin(pi), cos(ti) ];

function Q = compute_energy_from_angles(Saij,E)
%
% Saij: GxG, Saij(i,j) = Saij(j,i) = sin(alpha_ij)
% E: 1x1, exponent defining the energy term
Q = 1./(Saij.^E);
Q(logical(eye(size(Saij,1)))) = 0;
Q = sum(Q(:))/2;

function DQ = compute_energy_gradient_from_angles(cti,sti,cpi,spi,Caij,Saij,E,numerical)
%
% cti: Gx1, cos(thetai)
% sti: Gx1, sin(thetai)
% cpi: Gx1, cos(phii)
% spi: Gx1, sin(phii)
% Caij: GxG, Caij(i,j) = Caij(j,i) = cos(alpha_ij)
% Saij: GxG, Saij(i,j) = Saij(j,i) = sin(alpha_ij)
% E: 1x1, exponent defining the energy term
G  = size(cti,1);
if(numerical)
    DQ = zeros(2*G,1);
    tc = atan2(sti,cti);
    pc = atan2(spi,cpi);
    dh = pi/100000;
    for n=1:G
        tn = tc;
        tp = tc;
        tn(n) = tn(n) + dh/2;
        tp(n) = tp(n) - dh/2;
        [~,Saijn] = compute_solid_angles_from_spherical(cos(tn),sin(tn),cpi,spi);
        [~,Saijp] = compute_solid_angles_from_spherical(cos(tp),sin(tp),cpi,spi);
        DQ(n,1) = (compute_energy_from_angles(Saijn,E)-compute_energy_from_angles(Saijp,E))/dh;
        pn = pc;
        pp = pc;
        pn(n) = pn(n) + dh/2;
        pp(n) = pp(n) - dh/2;
        [~,Saijn] = compute_solid_angles_from_spherical(cti,sti,cos(pn),sin(pn));
        [~,Saijp] = compute_solid_angles_from_spherical(cti,sti,cos(pp),sin(pp));
        DQ(n+G,1) = (compute_energy_from_angles(Saijn,E)-compute_energy_from_angles(Saijp,E))/dh;
    end
else
    Gfact = E*(Caij./(Saij.^(E+2)));
    Gfact(logical(eye(G))) = 0;
    DQ1   = (cti*sti').*(cpi*cpi'+spi*spi') - sti*cti';
    DQ2   = -(sti*sti').*(spi*cpi'-cpi*spi');
    DQ    = sum([Gfact.*DQ1;Gfact.*DQ2],2);
end

function DDQ = compute_energy_hessian_from_angles(cti,sti,cpi,spi,Caij,Saij,E,numerical)
%
% cti: Gx1, cos(thetai)
% sti: Gx1, sin(thetai)
% cpi: Gx1, cos(phii)
% spi: Gx1, sin(phii)
% Caij: GxG, Caij(i,j) = Caij(j,i) = cos(alpha_ij)
% Saij: GxG, Saij(i,j) = Saij(j,i) = sin(alpha_ij)
% E: 1x1, exponent defining the energy term
G  = size(cti,1);
if(numerical)
    DDQ = zeros(2*G,2*G);
    tc = atan2(sti,cti);
    pc = atan2(spi,cpi);
    dh = pi/100000;
    for n=1:G
        tn = tc;
        tp = tc;
        tn(n) = tn(n) + dh/2;
        tp(n) = tp(n) - dh/2;
        [Caijn,Saijn] = compute_solid_angles_from_spherical(cos(tn),sin(tn),cpi,spi);
        DQn = compute_energy_gradient_from_angles(cos(tn),sin(tn),cpi,spi,Caijn,Saijn,E,true);
        [Caijp,Saijp] = compute_solid_angles_from_spherical(cos(tp),sin(tp),cpi,spi);
        DQp = compute_energy_gradient_from_angles(cos(tp),sin(tp),cpi,spi,Caijp,Saijp,E,true);
        DDQ(:,n) = (DQn-DQp)/dh;
        pn = pc;
        pp = pc;
        pn(n) = pn(n) + dh/2;
        pp(n) = pp(n) - dh/2;
        [Caijn,Saijn] = compute_solid_angles_from_spherical(cti,sti,cos(pn),sin(pn));
        DQn = compute_energy_gradient_from_angles(cti,sti,cos(pn),sin(pn),Caijn,Saijn,E,true);
        [Caijp,Saijp] = compute_solid_angles_from_spherical(cti,sti,cos(pp),sin(pp));
        DQp = compute_energy_gradient_from_angles(cti,sti,cos(pp),sin(pp),Caijp,Saijp,E,true);
        DDQ(:,n+G) = (DQn-DQp)/dh;
    end
else
    Gfact1 = E*(Caij./(Saij.^(E+2)));
    Gfact1(logical(eye(G))) = 0;
    H11_1  = Gfact1.*( (cti*cti').*(cpi*cpi'+spi*spi') + sti*sti' );
    H22_1  = Gfact1.*( (sti*sti').*(cpi*cpi'+spi*spi') );
    H12_1  = Gfact1.*( (cti*sti').*(spi*cpi'-cpi*spi') );
    Gfact2 = E*(Saij.*Saij+(E+2)*Caij.*Caij)./(Saij.^(E+4));
    Gfact2(logical(eye(G))) = 0;
    %%% -----
    H11_2  = Gfact2.* ...
        ( (cti*sti').*(cpi*cpi'+spi*spi') - sti*cti' ).* ...
        ( (sti*cti').*(cpi*cpi'+spi*spi') - cti*sti' );
    H22_2  = Gfact2.* ...
        ( -(sti*sti').*(spi*cpi'-cpi*spi') ).* ...
        (  (sti*sti').*(spi*cpi'-cpi*spi') );
    H12_2  = Gfact2.* ...
        ( (cti*sti').*(cpi*cpi'+spi*spi') - sti*cti' ).* ...
        (  (sti*sti').*(spi*cpi'-cpi*spi') );
    %%% -----
    dH11_1  = Gfact1.*( -(sti*sti').*(cpi*cpi'+spi*spi') - cti*cti' );
    dH22_1  = Gfact1.*( -(sti*sti').*(cpi*cpi'+spi*spi') );
    dH12_1  = Gfact1.*( -(cti*sti').*(spi*cpi'-cpi*spi') );
    dH11_2  = Gfact2.* ...
        ( (cti*sti').*(cpi*cpi'+spi*spi') - sti*cti' ).* ...
        ( (cti*sti').*(cpi*cpi'+spi*spi') - sti*cti' );
    dH22_2  = Gfact2.* ...
        ( -(sti*sti').*(spi*cpi'-cpi*spi') ).* ...
        ( -(sti*sti').*(spi*cpi'-cpi*spi') );
    dH12_2  = Gfact2.* ...
        (  (cti*sti').*(cpi*cpi'+spi*spi') - sti*cti' ).* ...
        ( -(sti*sti').*(spi*cpi'-cpi*spi') );
    d11 = sum(dH11_1,2)+sum(dH11_2,2);
    d22 = sum(dH22_1,2)+sum(dH22_2,2);
    d12 = sum(dH12_1,2)+sum(dH12_2,2);
    %%% -----
    H11_1(logical(eye(G))) = d11;
    H22_1(logical(eye(G))) = d22;
    H12_1(logical(eye(G))) = d12;
    %%% -----
    DDQ = [H11_1,H12_1;H12_1',H22_1]+[H11_2,H12_2;H12_2',H22_2];
end

function [Caij,Saij] = compute_solid_angles_from_spherical(cti,sti,cpi,spi)
%
% cti: Gx1, cos(thetai)
% sti: Gx1, sin(thetai)
% cpi: Gx1, cos(phii)
% spi: Gx1, sin(phii)
Caij = (sti*sti').*(cpi*cpi'+spi*spi') + cti*cti'; % GxG, sym.
Saij = sqrt(1-Caij.*Caij); % GxG, sym






