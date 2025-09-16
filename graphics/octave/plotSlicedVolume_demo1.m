function plotSlicedVolume_demo1(hf,volume,mask,slices,flipY)
figure(hf);

if(nargin<4)
    slices = [3,5,8,10,13,15]; % Slices to show
end

if(nargin<5)
    flipY = false;
end

if(nargin<3)
    mask = true(size(volume));
end
volume = volume.*mask;

IMG = [];
for k=1:length(slices)
    img = squeeze(volume(:,:,slices(k),:));
    % Transpose the image (swap xy -> ji) to meet Matlab
    % convention:
    img = permute(img,[2,1,3]);
    if(flipY)
        img = img(end:-1:1,:,:);
    end
    IMG = cat(2,IMG,img);
end

if(size(IMG,3)>1)
    % If the images is RGB, represent as it is:
    imshow(IMG);
else
    % Otherwise, find a proper dynamic range:
    window = quantile( volume(mask), [0.05,0.95] );
    if(window(1)==window(2))
        window(1) = 0;
    end
    % And plot:
    imshow(IMG,window);
    colorbar;
end
drawnow;
end
