function plotOverlaidLabelmap(labels,fa)
%slices = 7:4:54;
slices = [ ...
    3, 5, 7, 9, 11, ...
    30, 31, 32, 34, 35, ...
    36, 37, 51, 53, 55 ];
colors = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1;.8,.4,.2];
tiledlayout( 3, 5, 'TileSpacing', 'none' );
RGB = zeros( size(fa,2), size(fa,1), 3 );
for k=1:length(slices)
    nexttile;
    R  = fa(:,:,slices(k));
    G  = fa(:,:,slices(k));
    B  = fa(:,:,slices(k));
    LB = labels(:,:,slices(k));
    for l=1:7 % For each label
        R(LB==l) = colors(l,1);
        G(LB==l) = colors(l,2);
        B(LB==l) = colors(l,3);
    end
    RGB(end:-1:1,:,1) = R';
    RGB(end:-1:1,:,2) = G';
    RGB(end:-1:1,:,3) = B';
    imshow(RGB);
end
drawnow;
end
