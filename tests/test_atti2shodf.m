% test_dwi2freewater3
rng(13,'twister');

if(exist('quaternion','file')~=2)
    error('This test script uses the quaternion class to produce uniform random rotations: http://www.lpi.tel.uva.es/node/626');
end

N = 6; % Number of synthetic voxels to be tested
PSNR = 20; % Peak signal to noise ratio in the baseline
sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI (noiseless baseline assumed)
L = 8; % Order of the spherical harmonics
lambda = 0.006; % Tikhonov regularisation for Spherical Harmonics
type = 'copdt';
b0   = 3000;

q1 = [1,ii,jj,kk]*randn(4,N); % 1xN
q1 = q1./abs(q1); % Random rotations for the first compartment


[gi,~,~] = icosamplesSphere(5,'O1',true); % Gx3
G  = size(gi,1);
bi = b0*ones(G,1);

x1 = imag(q1.*ii.*conj(q1)); % 3xN
y1 = imag(q1.*jj.*conj(q1)); % 3xN
z1 = imag(q1.*kk.*conj(q1)); % 3xN
x2 = imag(q1.*jj.*conj(q1)); % 3xN
y2 = imag(q1.*kk.*conj(q1)); % 3xN
z2 = imag(q1.*ii.*conj(q1)); % 3xN
l1 = 1.7e-3;
l2 = 0.4e-3;
S  = zeros(G,N);
for n=1:N
    D1 = l1*(x1(:,n)*x1(:,n)')+l2*(y1(:,n)*y1(:,n)')+l2*(z1(:,n)*z1(:,n)');
    D2 = l1*(x2(:,n)*x2(:,n)')+l2*(y2(:,n)*y2(:,n)')+l2*(z2(:,n)*z2(:,n)');
    S1 = D1*(gi');
    S2 = D2*(gi');
    S1 = bi.*sum(gi.*S1',2);
    S2 = bi.*sum(gi.*S2',2);
    S1 = exp(-S1);
    S2 = exp(-S2);
    S(:,n) = (S1+S2)/2;
end

noise1 = sigma*(randn(G,N)+1i*randn(G,N));
noise2 = sigma*(randn(G,N)+1i*randn(G,N));
S      = abs(S+noise1)./abs(1+noise2);
S      = reshape(S',[N,1,1,G]);

%%% -----------------------------------------------------------------------
NV   = 1000;
nvts = [6,12,18,42,66,162,258,642,1026,2562];
src  = [1,2,1,2,1,2,1,2,1,2];
pos  = find(nvts>=NV,1,'first');
src  = src(pos);
pos  = ceil(pos/2)-1;
if(src==1) % From octahedron
    polyhedron = createOctahedron;
else % From icosahedron
    polyhedron = createIcosahedron;
end
for l=1:pos % Refine as much as needed
    polyhedron = refinePolyhedron(polyhedron);
    % Improve vertices distribution:
    %%% -------------------------------------------------------------------
    grads = polyhedron.vertices;
    moduli = sqrt(sum(grads.*grads,2));
    grads  = grads./repmat(moduli,[1,3]);
    % Find dot-products of each gradient with the remaining ones:
    dp   = grads(:,1)*(grads(:,1)') + grads(:,2)*(grads(:,2)') + grads(:,3)*(grads(:,3)');
    dp   = abs(dp); % The closer to 1, the more similar
    dp   = dp - diag(diag(dp)); % Avoid self-similarity
    idx1 = 1:size(grads,1);
    idx2 = 1:(size(grads,1)/2);
    for n=1:size(grads,1)/2
        pos = idx1(n);
        [~,mpos] = max(dp(pos,:)); % mpos is the most similar gradient to n, so...
        idx1 = setdiff(idx1,mpos); % ... remove it from the list
        idx2(n) = mpos;            % store the most similar gradient to idx(n)
    end
    grads = grads(idx1,:);
    vts   = grads;
    %%% -------------------------------------------------------------------
    vts2 = optimizeGradients(vts,'verbose',false);
    vts(idx1,:) =  vts2;
    vts(idx2,:) = -vts2;
    polyhedron.vertices = vts;
end
%%% -----------------------------------------------------------------------

odf  = atti2odf( S, gi, polyhedron.vertices, 'type', type, 'L', L, 'lambda', lambda, 'O1', 0 );
odf(odf<0) = 0;
vts   = polyhedron.vertices; % original vertices over the unit sphere

for n=1:N
    subplot(2,3,n);
    hold('on');
    odfn = squeeze(odf(n,1,1,:));
    DD   = max(odfn(:))*1.2;
    vts2   = repmat(odfn,[1,3]).*vts; % Distances are weighted depending on the ODF
    patch( ...
        'Vertices', vts2, ...
        'Faces',    polyhedron.facets, ...
        'EdgeColor', [0,0,0], ...
        'EdgeAlpha', 0, ...
        'FaceVertexCData', abs(vts), ...
        'FaceColor', 'flat', ...
        'FaceLighting', 'phong' ...
        );
    light;
    plot3([-x1(1,n),x1(1,n)]*DD,[-x1(2,n),x1(2,n)]*DD,[-x1(3,n),x1(3,n)]*DD);
    plot3([-x2(1,n),x2(1,n)]*DD,[-x2(2,n),x2(2,n)]*DD,[-x2(3,n),x2(3,n)]*DD);
    axis('equal');
    rotate3d('on');
    grid('on');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end
