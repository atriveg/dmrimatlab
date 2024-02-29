% test_dwi2freewater3
rng(23,'twister');

PSNR = 13.3; % Peak signal to noise ratio in the baseline
sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI (noiseless baseline assumed)
L = 6; % Order of the spherical harmonics
lambda = 0.006; % Tikhonov regularisation for Spherical Harmonics
type = 'opdt';
b0   = 1200;

RR = [...
    -0.791198578757938, -0.258568926543925, -0.554208371641335; ...
    -0.610518180022741,  0.386808203471646,  0.691120080440987; ...
     0.035670167270954,  0.885167511803816, -0.463903131282709];

[gi,~,~] = icosamplesSphere(5,'O1',true); % Gx3
G  = size(gi,1);
bi = b0*ones(G,1);

x1  = [1;0;0];
y1  = [0;1;0];
z1  = [0;0;1];
x2  = [0;1;0];
y2  = [0;0;1];
z2  = [1;0;0];
alp = pi*(70/180);
x2p = cos(alp)*x2 - sin(alp)*y2;
y2p = sin(alp)*x2 + cos(alp)*y2;
x2  = x2p;
y2  = y2p;
%%% ------------------------
x1 = RR*x1; y1 = RR*y1; z1 = RR*z1;
x2 = RR*x2; y2 = RR*y2; z2 = RR*z2;
%%% ------------------------
l1 = 1.7e-3;
l2 = 0.3e-3;
D1 = l1*(x1*x1') + l2*(y1*y1')+l2*(z1*z1');
D2 = l1*(x2*x2') + l2*(y2*y2')+l2*(z2*z2');
S1 = D1*(gi');
S2 = D2*(gi');
S1 = bi.*sum(gi.*S1',2);
S2 = bi.*sum(gi.*S2',2);
S1 = exp(-S1);
S2 = exp(-S2);
S  = (S1+S2)/2;

RR = [...
    -0.791198578757938, -0.258568926543925, -0.554208371641335; ...
    -0.610518180022741,  0.386808203471646,  0.691120080440987; ...
     0.035670167270954,  0.885167511803816, -0.463903131282709];
 
noise1r = sigma*randn(G,1);
noise1i = sigma*randn(G,1);
noise2r = sigma*randn(G,1);
noise2i = sigma*randn(G,1);
S       = abs(S+noise1r+1i*noise1i)./abs(1+noise2r+1i*noise2i);

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

vts  = polyhedron.vertices; % original vertices over the unit sphere
S    = reshape(S,[1,1,1,G]);

[odf,sh1]  = atti2odf( S, gi, polyhedron.vertices, 'type', type, 'L', L, 'lambda', lambda, 'O1', 0 );
odf  = squeeze(odf);
odf(odf<0) = 0;
subplot(1,2,1);
hold('on');
DD   = max(odf)*1.2;
vts2 = repmat(odf,[1,3]).*vts; % Distances are weighted depending on the ODF
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
plot3([-x1(1),x1(1)]*DD,[-x1(2),x1(2)]*DD,[-x1(3),x1(3)]*DD);
plot3([-x2(1),x2(1)]*DD,[-x2(2),x2(2)]*DD,[-x2(3),x2(3)]*DD);
axis('equal');
rotate3d('on');
grid('on');
xlabel('x');
ylabel('y');
zlabel('z');

[odf,sh2]  = atti2odf( S, ((RR')*(gi'))', ((RR')*(polyhedron.vertices'))', 'type', type, 'L', L, 'lambda', lambda, 'O1', 0 );
odf  = squeeze(odf);
odf(odf<0) = 0;
subplot(1,2,2);
hold('on');
DD   = max(odf)*1.2;
vts2 = repmat(odf,[1,3]).*vts; % Distances are weighted depending on the ODF
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
plot3([-x1(1),x1(1)]*DD,[-x1(2),x1(2)]*DD,[-x1(3),x1(3)]*DD);
plot3([-x2(1),x2(1)]*DD,[-x2(2),x2(2)]*DD,[-x2(3),x2(3)]*DD);
axis('equal');
rotate3d('on');
grid('on');
xlabel('x');
ylabel('y');
zlabel('z');
