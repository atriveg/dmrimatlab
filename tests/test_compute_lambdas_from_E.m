% test_compute_lambdas_from_E

b       = [500;1000;2500];
lperp0  = linspace(0,3,100).*1.0e-3;
deltal0 = linspace(0,2,100).*1.0e-3;
[lperp0,deltal0] = meshgrid(lperp0,deltal0);
lpar0   = lperp0 + deltal0;

val1    = lperp0.*permute(b,[2,3,1]);
val2    = sqrt(deltal0.*permute(b,[2,3,1]));
val1    = exp(-val1);
val2    = (sqrt(pi)/2)*erf(val2)./val2;
val2(isnan(val2)) = 1;
val2(isinf(val2)) = 1;
E       = val1.*val2;

ADC0MAX = 3.0e-3;
mlperp  = 0.05e-3;

mu      = 5.0e-5;
NMAX    = 400;
[lpar,lperp,nit] = compute_lambdas_from_E(E,b,mu,mlperp,NMAX,1.0e-9,1.0e-9,ADC0MAX);

errpar  = abs(lpar-lpar0);
errperp = abs(lperp-lperp0);

mask          = ( (lpar0>ADC0MAX) | (lperp0>ADC0MAX) | (lperp0<mlperp) );
errpar(mask)  = 0;
errperp(mask) = 0;

figure(1001);
subplot(2,3,1); imagesc(lperp0(1,:),deltal0(:,1)',lperp0);  colorbar; xlabel('\lambda_{\perp}'); ylabel('\lambda_{||}-\lambda_{\perp}'); title('True \lambda_{\perp}');
subplot(2,3,2); imagesc(lperp0(1,:),deltal0(:,1)',lperp);   colorbar; xlabel('\lambda_{\perp}'); ylabel('\lambda_{||}-\lambda_{\perp}'); title('Estimated \lambda_{\perp}');
subplot(2,3,3); imagesc(lperp0(1,:),deltal0(:,1)',errperp); colorbar; xlabel('\lambda_{\perp}'); ylabel('\lambda_{||}-\lambda_{\perp}'); title('Error in \lambda_{\perp}');
subplot(2,3,4); imagesc(lperp0(1,:),deltal0(:,1)',lpar0);   colorbar; xlabel('\lambda_{\perp}'); ylabel('\lambda_{||}-\lambda_{\perp}'); title('True \lambda_{||}');
subplot(2,3,5); imagesc(lperp0(1,:),deltal0(:,1)',lpar);    colorbar; xlabel('\lambda_{\perp}'); ylabel('\lambda_{||}-\lambda_{\perp}'); title('Estimated \lambda_{||}');
subplot(2,3,6); imagesc(lperp0(1,:),deltal0(:,1)',errpar);  colorbar; xlabel('\lambda_{\perp}'); ylabel('\lambda_{||}-\lambda_{\perp}'); title('Error in \lambda_{||}');

figure(2002);
subplot(2,3,1)
imagesc(nit,[0,NMAX]); colormap; colorbar;
subplot(2,3,2)
imagesc(nit,[0,NMAX/2]); colormap; colorbar;
subplot(2,3,3)
imagesc(nit,[0,NMAX/4]); colormap; colorbar;
subplot(2,3,4)
imagesc(nit,[0,NMAX/8]); colormap; colorbar;
subplot(2,3,5)
imagesc(nit,[0,NMAX/16]); colormap; colorbar;
subplot(2,3,6)
imagesc(nit,[0,NMAX/32]); colormap; colorbar;






