function plotOverlaidLabelmap_demo1(labels,fa)

slices = [ ...
    3, 5, 7, 9, 11, ...
    30, 31, 32, 34, 35, ...
    36, 37, 51, 53, 55 ];
colors = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;.8,.4,.2];
IMG = [];
RGB = zeros( size(fa,2), size(fa,1), 3 );

fa(isnan(fa)) = 0;
fa(isinf(fa)) = 0;

pos = 1;
for r = 1:3
    ROW = [];
    for c = 1:5
        % -------------------------
        R  = fa(:,:,slices(pos));
        G  = fa(:,:,slices(pos));
        B  = fa(:,:,slices(pos));
        LB = labels(:,:,slices(pos));
        % -------------------------
        for l=1:7 % For each label
            R(LB==l) = colors(l,1);
            G(LB==l) = colors(l,2);
            B(LB==l) = colors(l,3);
        end
        % -------------------------
        RGB(end:-1:1,:,1) = R';
        RGB(end:-1:1,:,2) = G';
        RGB(end:-1:1,:,3) = B';
        % -------------------------
        ROW = cat(2,ROW,RGB);
        % -------------------------
        pos = pos+1;
    end
    IMG = cat(1,IMG,ROW);
end
imshow(IMG);
drawnow;
end
