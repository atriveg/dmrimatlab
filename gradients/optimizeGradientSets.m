function gi = optimizeGradientSets(gi,sets,varargin)
% function gi = optimizeGradientSets(gi,sets,'opt1',value1,'opt2',value2)
%
%    This function is used to optimally distribute N sets of gradient
%    directions so that each set is (roughly) uniformly distributed
%    and their combinations are uniformly distributed as well. Here,
%    "uniformly distributed" means the energy term:
%
%         Q[G] = 1/2 * sum_{i\not=j} 1/sin^p(alpha_ij) [1]
%
%    is as small as possible (alpha_ij is the angle between gradients g_i
%    and g_j in a given collection G made up of one or several sets of
%    gradients).
%
%     - Let G_1, G_2, ..., G_N be the N sets of gradient directions.
%     - Let G_{i1,i2,...,iK} = Union(G_i1,G_i2,...,G_iK).
%     - Let y_{n,j} be each of the (N over n) combinations of n elements 
%       taken from 1,2,...,N.
%
%    Then, from the definition in eq. [1], we minimize:
%
%      C = sum_{n=1..N} 
%               w(n) * sum_{j=1..(N over n)} Q[G_{y_{n,j}}]/sc_{n,j} [2]
%
%    i.e.:
%
%     - w(1) penalizes the individual costs of each set.
%     - w(2) penalizes pair-wise combinations.
%     - w(3) penalizes three-fold combinations.
%     - w(end) penalizes the cost of the N sets altogether.
%
%    The scales sc_{n,j} are internally computed so that all possible
%    combinations of sets have roughly the same impact in the overall
%    cost regardless of the total number of gradients they include.
%    
%    This function can be used for two purposes: (1) designing multi-shell
%    sampling schemes for which each shell fills the gaps in the
%    orientations space not covered by the others or (2) designing
%    single-shell sampling schemes that can be bootstrapped for validation
%    purposes.
%
%     EXAMPLES:
%
%       >> Ns     = [8,16,32,32,64,64];
%       >> [G,ps] = spiralPhyllotaxis(Ns);
%       >> G      = optimizeGradientSets( G, ps, ...
%                         'weights', [1000,100,0,0,0,1], ...
%                         'verbose', true );
%
%       (Taken from test_spiralPhyllotaxis.m, you may run the script to
%       visually check the results. See also test_spiralPhyllotaxis2.m 
%       and test_spiralPhyllotaxis3.m).
%
%
%     INPUTS:
%
%       gi: Gx3 matrix of normalized gradients (just in the north
%           hemisphere, since it is assuemd that gi(i,:) implies -gi(i,:))
%       sets: Gx1, vector of integers in the range 1..S, each entry
%           indicates the set each gradient belongs to. It is important
%           that the variable fulfills unique(sets) = [1;2;...;S], with
%           S>1.
%
%     OUTPUTS:
%
%       gi: Gx3, the optimized gradients.
%
%     Options:
%
%       p: 1x1 double, it defines the energy term in [1] (default: 1)
%       weights: Sx1, with S the number of different sets of gradients.
%          These are the weights applied to the different combinations of
%          gradients in eq. [2] (default: [1;0;0;...;1])
%       nocomb: if set true, then not all possible combinations of sets,
%          but only consecutive ones, are mixed together in eq. [2]. For
%          example, pair-wise costs are computed over gradient sets {1,2},
%          {2,3}, {3,4}, ..., and {N-1,N}, but not for {1,3}, {1,4}, {2,4},
%          etcetera. This option might be helpful to design progressive
%          samplings where a new set complements the previous ones
%          (default: false)
%
%     The optimization uses a Newton-Raphson/Levenberg-Marquardt
%     procedure to minimize the cost in [2], for which the following
%     parameters apply:
%
%       iters: 1x1 double, maximum number of iterations for the
%          optimization (default: 1000)
%       lambda: 1x1 double, Levenberg-Marquardt regularization term in the 
%          Newton-Raphson step. This is just the initial value, it is
%          adaptively updated as the iterations proceed (default: 10.0)
%       delta: maximum allowed change in the parameters before the
%          iterations stop (default: 1.0e-9)
%
%       verbose: 1x1 logical, print additional information (default: false)

% ----------------------------------------------------------------------------
if(nargin<2)
    error('At least gi and sets must be provided');
end
assert(size(gi,2)==3,'gi must be Gx3');
sets = sets(:);
assert(size(gi,1)==numel(sets),'The number of elements in sets must match the number of entries in gi');
[sidx,ngs,~] = unique(sort(sets),'last');
ngs  = diff([0;ngs]);
S    = max(sidx);
if(S==1)
    error('At least two different sets must be provided');
end
if(~isequal(sidx,(1:S)'))
    error('The entries in "sets" should be in the range {1, 2, ..., S}');
end
weights      = zeros(S,1);
weights(1)   = 1;
weights(end) = 1;
% ----------------------------------------------------------------------------
G = size(gi,1);
% Parse the optional input arguments:
opt.weights = weights; optchk.weights = [true,false];
opt.nocomb = false;    optchk.nocomb = [true,true];
opt.p = 1;             optchk.p = [true,true];          % always 1x1 double
opt.iters = 1000;      optchk.iters = [true,true];      % always 1x1 double
opt.lambda = 10.0;     optchk.lambda = [true,true];     % always 1x1 double
opt.delta = 1.0e-9;    optchk.delta = [true,true];      % always 1x1 double
opt.verbose = false;   optchk.verbose = [true,true];    % always 1x1 boolean
opt = custom_parse_inputs(opt,optchk,varargin{:});
% ----------------------------------------------------------------------------
if(numel(opt.weights)~=numel(weights))
    error('The weights provided do not match the number of unique elements in "sets"');
end
% ----------------------------------------------------------------------------
xi = gi(:,1);
yi = gi(:,2);
zi = gi(:,3);
ti = atan2( sqrt(xi.*xi+yi.*yi), zi );
pj = atan2( yi, xi );
% ----------------------------------------------------------------------------
combs = cell(S,1);
if(opt.nocomb)
    for s=1:S
        ccomb = zeros(S-s+1,s);
        for r=1:S-s+1
            ccomb(r,:) = r:r+s-1;
        end
        combs{s,1} = ccomb;
    end
else
    for s=1:S
        combs{s,1} = combnk(1:S,s);
    end
end
% ----------------------------------------------------------------------------
scales = cell(S,1);
for s=1:S
    ccomb = combs{s,1};
    scale = zeros( size(ccomb,1), 1 );
    for r=1:size(ccomb,1)
        gtmp = spiralPhyllotaxis(   sum( ngs(ccomb(r,:)) )   );
        xtmp = gtmp(:,1);
        ytmp = gtmp(:,2);
        ztmp = gtmp(:,3);
        ttmp = atan2( sqrt(xtmp.*xtmp+ytmp.*ytmp), ztmp );
        ptmp = atan2( ytmp, xtmp );
        cttmp = cos(ttmp);
        sttmp = sin(ttmp);
        cptmp = cos(ptmp);
        sptmp = sin(ptmp);
        [~,Stmp] = compute_solid_angles_from_spherical(cttmp,sttmp,cptmp,sptmp);
        % Compute the overall cost:
        Qtmp  = compute_energy_from_angles(Stmp,opt.p);
        scale(r,1) = Qtmp;
    end
    scales{s,1} = scale;
end
% ----------------------------------------------------------------------------
if(opt.verbose)
    txtlength = 0;
end
idxs = (2:G)';
idxl = [idxs;idxs+G];
lmb  = opt.lambda;
for it=1:opt.iters
    % ------------------------------------------------------
    % Initiallize the cost, the gradient and the Hessian
    Q   = 0;
    DQ  = zeros(2*G,1);
    DDQ = zeros(2*G,2*G);
    % ------------------------------------------------------
    % Compute the global cost, gradient and Hessian by summing up all
    % possible combinations
    for s=1:S
        % Check if the weight associated to this combination is greater
        % than 0, otherwise ignore
        if(opt.weights(s)>eps)
            % Need to proceed
            crrcomb = combs{s,1};   % Each row is a combination
            crscale = scales{s,1};  % Each row is a scaling factor for the cost
            for r=1:size(crrcomb,1) % So, for each combination
                gidx = false(G,1);
                for c=1:size(crrcomb,2)
                    gidx = ( gidx | (abs(sets-crrcomb(r,c))<0.1) );
                end
                gidxd = [gidx;gidx];
                % Now gidx is true for all the gradients involved in the
                % present combination
                % -------------
                % Compute the contribution of this combination to the
                % overall cost, gradient, and Hessian
                % Compute tigonometric functions of theta and phi
                cti = cos(ti(gidx));
                sti = sin(ti(gidx));
                cpi = cos(pj(gidx));
                spi = sin(pj(gidx));                
                % Compute the angles between each pair of gradients:
                [Caij,Saij] = compute_solid_angles_from_spherical(cti,sti,cpi,spi);
                % Compute the overall cost:
                Qc   = compute_energy_from_angles(Saij,opt.p)/crscale(r);
                % Compute the gradient of the cost:
                DQc  = compute_energy_gradient_from_angles(cti,sti,cpi,spi,Caij,Saij,opt.p,false)/crscale(r);
                % Compute the Hessian of the cost:
                DDQc = compute_energy_hessian_from_angles(cti,sti,cpi,spi,Caij,Saij,opt.p,false)/crscale(r);
                % Sum to the overall terms:
                Q = Q + (opt.weights(s))*Qc;
                DQ(gidxd,1)  = DQ(gidxd,1) + (opt.weights(s))*DQc;
                DDQ(gidxd,gidxd) = DDQ(gidxd,gidxd) + (opt.weights(s))*DDQc;
            end
        end
    end
    % ------------------------------------------------------
    % (Regularized) Newton-Raphson step. One gradient direction will remain
    % fixed, and the remaining ones will update:
    LSM     = DDQ(idxl,idxl);
    LSMnorm = sqrt( trace(LSM*LSM')/(size(LSM,1)) );
    LSM     = LSM + lmb*LSMnorm*eye(size(LSM,1));
    LSV     = DQ(idxl);
    parc    = [ti(idxs);pj(idxs)];
    parn    = parc - LSM\LSV;
    % ------------------------------------------------------
    % Check if the cost function actually decreases:
    tin         = ti;
    pin         = pj;
    tin(idxs)   = parn( 1:length(parn)/2 );
    pin(idxs)   = parn( (1:length(parn)/2) + length(parn)/2 );
    % Recompute the cost using all combinations:
    Qn          = 0;
    for s=1:S
        % Check if the weight associated to this combination is greater
        % than 0, otherwise ignore
        if(opt.weights(s)>eps)
            % Need to proceed
            crrcomb = combs{s,1};   % Each row is a combination
            crscale = scales{s,1};  % Each row is a scaling factor for the cost
            for r=1:size(crrcomb,1) % So, for each combination
                gidx = false(G,1);
                for c=1:size(crrcomb,2)
                    gidx = ( gidx | (abs(sets-crrcomb(r,c))<0.1) );
                end
                % -------------
                cti = cos(tin(gidx));
                sti = sin(tin(gidx));
                cpi = cos(pin(gidx));
                spi = sin(pin(gidx));
                % -------------
                [~,Saij] = compute_solid_angles_from_spherical(cti,sti,cpi,spi);
                Qc  = compute_energy_from_angles(Saij,opt.p)/crscale(r);
                % -------------
                Qn  = Qn + (opt.weights(s))*Qc;
            end
        end
    end
    % Check:    
    if(Qn<=Q)
        ti    = tin;
        pj    = pin;
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
    
gi = [ sin(ti).*cos(pj), sin(ti).*sin(pj), cos(ti) ];

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
Saij = sqrt(1-min(Caij.*Caij,1)); % GxG, sym






