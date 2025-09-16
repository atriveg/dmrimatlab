function plotVolBW_demoDTI(V,scale)
if(nargin<2)
    scale = [];
end
if(isempty(scale))
    scale = quantile(V(:),[0.025,0.975]);
end
z   = [3,9,14];
IMG = [V(:,:,z(1))',V(:,:,z(2))',V(:,:,z(3))'];
imshow(IMG,scale);
colormap(parula);
colorbar;
end
