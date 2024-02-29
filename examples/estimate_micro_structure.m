% estimate_micro_structure
clear;

% This script illustrates a simple pipeline to estimate the percentage of
% confined water (i.e. intracellular) in the white matter together with a
% simple model for micro-structure diffusion based on an impulse response
% with a parallel diffusivity and a perpendicular diffusivity.

% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose your favorite to load a multi-shell diffusion volume. At least 3
% different shells must be available to compute the 3 parameters of the
% model, i.e.: the fraction of intra-axon water, the parallel diffusivity
% and the perpendiculas diffusivity:
PATH = '/media/atriveg/DATOS/Documentos/DMRIDATA/HCP_MGH_1007/';
filename1 = [PATH,'HCPMGH_1007_B1000.mat'];
filename2 = [PATH,'HCPMGH_1007_B3000.mat'];
filename3 = [PATH,'HCPMGH_1007_B5000.mat'];
filename4 = [PATH,'HCPMGH_1007_baseline.mat'];
% -----
file1 = load(filename1);
fprintf(1,'Loaded first shell\n');
file2 = load(filename2);
fprintf(1,'Loaded second shell\n');
file3 = load(filename3);
fprintf(1,'Loaded third shell\n');
file4 = load(filename4);
fprintf(1,'Loaded baseline\n');
% -----
bi    = [file1.bvals1000'; file2.bvals3000'; file3.bvals5000']; % Gx1
gi    = [file1.bvecs1000;  file2.bvecs3000;  file3.bvecs5000];  % Gx3
% -----
[M,N,P,G1] = size(file1.dwi_b1000);
[~,~,~,G2] = size(file2.dwi_b3000);
[~,~,~,G3] = size(file3.dwi_b5000);
% -----
dwi      = zeros(M,N,P,G1+G2+G3);
baseline = mean(file4.baseline,4);
dwi(:,:,:,1:G1)             = file1.dwi_b1000;
dwi(:,:,:,G1+1:G1+G2)       = file2.dwi_b3000;
dwi(:,:,:,G1+G2+1:G1+G2+G3) = file3.dwi_b5000;
fprintf(1,'Arranged the whole DWI. Computing attenuations...\n');
% -----
if(is_broadcast_available)
    atti = dwi./baseline;
else
    atti = bsxfun( @(x,y)(x./y), dwi, baseline );
end
atti(isinf(atti)) = 0;
atti(isnan(atti)) = 0;
fprintf(1,'Computed attenuations\n');
% -----
dwi = dwi(:,:,:,1:G1);
dwi(:,:,:,2:G1+1) = dwi;
dwi(:,:,:,1)      = baseline;
fprintf(1,'Stripped first shell\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% At this point we have volume atti with size MxNxPxG with the whole
% (attenuation) data set, bi with size Gx1 with the corresponding b-values
% (all shells), and gi with size Gx3 with the corresponding gradient
% directions. You shouldn't need to change the code below this point.

% Mask the volume to avoid unnecessary computations:
mask       = (mean(dwi,4)>50);
clear('dwi');
[xx,yy,zz] = ndgrid(-5:5);
nhood      = (sqrt(xx.^2 + yy.^2 + zz.^2)<=5.0);
mask       = imdilate(mask,nhood);
mask       = imerode(mask,nhood);

fprintf(1,'Computed mask. Processing...\n');

[lpar,lperp,f] = atti2micro( atti, gi, bi, ...
    'tl', 1.0e-6, 'tu', 1-1.0e-6, 'ADC0', 3.0e-3, ...
    'usef', true, 'mask', mask, 'bth', 100, ...
    'verbose', true ...
    );
clear('atti');
fprintf(1,'All done\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lpar(lpar>3.0e-3) = 0;
lpar(lpar<0)      = 0;
lperp(lperp>lpar) = 0;
lperp(lperp<0)    = 0;
tau   = 20.0e-3;
uRTAP = (1/(4*pi*tau))./sqrt(lperp.*lperp);
uRTPP = (1/sqrt(4*pi*tau))./sqrt(lpar);
uRTOP = uRTAP.*uRTPP;
mul   = (lpar+2*lperp)/3;
uFA   = sqrt(1/2)*sqrt(2*(lpar-lperp).*(lpar-lperp))./sqrt(lpar.*lpar+2*lperp.*lperp);
%%%%%
uRTAP(~mask) = 0;
uRTPP(~mask) = 0;
uRTOP(~mask) = 0;
uFA(~mask)   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 2;
C = 3;
N = R*C;
SL  = size(f,3);
SL0 = 15;
dsl = round((SL-SL0)/(N+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(R,C,1);
for n=1:N
    subplot(R,C,n);
    sl = SL - dsl*n;
    imagesc(f(:,:,sl)');
    colorbar;
    title(['f at slice ',num2str(sl)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
subplot(R,C,1);
for n=1:N
    subplot(R,C,n);
    sl = SL - dsl*n;
    imagesc(lpar(:,:,sl)');
    colorbar;
    title(['l_{||} at slice ',num2str(sl)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
subplot(R,C,1);
for n=1:N
    subplot(R,C,n);
    sl = SL - dsl*n;
    imagesc(lperp(:,:,sl)');
    colorbar;
    title(['l_{\perp} at slice ',num2str(sl)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
subplot(R,C,1);
for n=1:N
    subplot(R,C,n);
    sl = SL - dsl*n;
    imagesc(uFA(:,:,sl)',[0.6,1]);
    colorbar;
    title(['uFA at slice ',num2str(sl)]);
end

    
    
    
    
    