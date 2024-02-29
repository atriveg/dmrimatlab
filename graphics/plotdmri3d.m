function [ha,ho] = plotdmri3d(SH,i1,i2,i3,varargin)
% function [ha,ho] = plotdmri3d(SH,i1,i2,i3,'opt1',val1,'opt2',val2,...)
%
%    This is a rather complex function to represent diffusion MRI data in
%    3D, either scalars (FA, MD, ...), color-coded information, and/or 
%    glyphs (DTI-like ellipsoids or ODF-like glyphs). The function is
%    mainly intended to work with data read from NRRD files, so that we use
%    the same names and conventions as in the NRRD standard. Yor may check
%    some of the case uses in the examples folder for some further
%    information.
%
%        SH: This is a MxNxPxK array of doubles with SH coefficients
%            describing some symmetric orientation function, so that
%            K=(L+1)(L+2)/2 for L=0,2,4,... This volume should be obtained
%            with some of the functions in the sh folder (probably
%            atti2odf.m or atti2shadc.m; the latter, with L=2, will mimick
%            a DTI estimation).
%               NOTE: this volume is used as a reference and cannot be left
%               empty. In case you don't want glyphs to be represented, use
%               the option 'glyphs' described below.
%
%        i1, i2, i3: 1x2 arrays of integers each describing the piece of
%            data (slice) to plot with each call to this function. One of
%            these arguments must describe a singleton dimension, i.e.,
%            either i1(1)=i1(2), or i2(1)=i2(2), o i3(1)=i3(2). A slice of
%            the SH volume taken as SH(i1(1):i1(2),i2(1):i2(2),i3(1):i3(2))
%            will be represented in the 3-D view.
%
%        ha: 1x1 double with the axes handle where the graphics will be
%            plotted (can be passed in subsequent calls to the function to
%            stack several slices/representations in one figure).
%
%        ho: a 1x1 struct of handles to the graphical objects created by
%            the function when it is called. Each field may or may not be
%            present depending on the options passed (i.e. depending if
%            each of the following items are actually created, see below):
%               ho.bbox: handles to the bounding box
%               ho.rastxt: handles to the R/L, A/P, I/S flags
%               ho.bgimage: handle to the surface of the backgound image
%               ho.glyphs: handles to the glyphs objects
%
%
%   ADDITIONAL ARGUMENTS PASSED AS 'OPTION',VALUE PAIRS
%
%   + Anatomy-related:
%
%        origin: 3x1 vector of doubles (default: [0;0;0])
%        direction: 3x3 matrix of doubles (default: eye(3))
%            These two parameters together relate the index space within
%            the (3+1)-D array SH and the real world, physical coordinates
%            of the imaged volume; if m\in[i1(1),i1(2)], n\in[i2(1),i2(2)],
%            and p\in[i3(1),i3(2)], then the voxel (m,n,p) corresponds to
%            the physical location:
%               [x;y;z] = [ direction, origin ]*[m;n;p;1]
%        space: a string of the form: RAS, LPS, RAI or
%            right-anterior-superior, left-posterior-superior,
%            right-anterior-inferior (8 combinations of right/left,
%            anterior/posterior, superior/inferior) describing the
%            anatomical convention. R means the increasing direction of the
%            x axis goes from left to right; A means the increasing
%            direction of the y axis goes from back to front; S means the
%            increasing direction o the z axis goes upwards (default: RAS,
%            which is the same as right-anterior-superior)
%        mframe: the measurment frame, a 3x3 rotation matrix that relates
%            the diffusion measurements to the physical space (default: 
%            eye(3))
%
%   + General plot options:
%
%        ha: axes where to plot the current slice (may be combined with the
%            output of this function to stack several slices in one
%            figure). Leave empty to create new axes (default: [])
%        bbox: whether (true) or not (false) plotting a 3-D bounding box
%            with the limits of the dMRI volume in both its own i-j-k space
%            and the RAS space together with 6 flags indicating the
%            left/right, anterior/posterior, and inferior/superior
%            directions (default: true)
%        rastxt: if the bounding box is drawn, this is the font size of the
%            R/L, A/P, I/S flags used to indicate anatomical directions. If 
%            set to 0, just remove the flags (default: 20)
%        mask: a MxNxP logical used to mask the plot. Voxels where the mask
%            is false are removed from the plot (default: true(M,N,P))
%
%   + Plot options related to the scalar/color-coded background image:
%
%        bgsh: this is a MxNxPxK2 image of SH that may be used to plot a
%            background image (MD, FA, color-coded, others) along with the
%            desired glyphs. It has to be either an empty matrix, [], or its
%            first three dimensions should match those of the first input
%            to the function, SH (see below) (default: [])
%        bgimage: selects the background (scalar or color-coded) image.
%            Its value may be:
%              - 'none' (default) to avoid plotting any background image.
%                In this case bgsh should be left empty ('bgsh',[])
%              - a string representing some sort of scalar diffusion
%                measurement, such as: 'md', 'fa', 'gfa', or 'color', which
%                will be computed from bgsh (then, it cannot be left empty)
%              - a whole MxNxP volume of doubles with the actual scalar
%                volume to represent, whose first three sizes M, N, P, must
%                match those of the first input to th function, SH. In this
%                case bgsh should be left empty ('bgsh', []).
%            Use cases of these two arguments together are:
%              plotdmri3d(..., 'bgsh', SH, 'bgimage', 'fa', ... )
%              plotdmri3d(..., 'bgsh', SH, 'bgimage', 'color', ... )
%              plotdmri3d(..., 'bgsh', [], 'bgimage', 'none', ... )
%              plotdmri3d(..., 'bgsh', [], 'bgimage', VOL, ... )
%            where size(SH)=[M,N,P,K2] (K2=6,15,28,45...) and
%            size(VOL)=[M,N,P].
%        bgalpha: 1x1 scalar within the range [0,1] with the alpha value to
%            apply to the background image (default: 1)
%
%   + Plot options related to the dMRI glyphs:
%
%        glyphs: whether (true) or not (false) actually plot the
%            ellipsoids/ODF-like glyphs for the SH volume (default: true)
%        clip: whether (true) or not (false) clip negative values of the
%            ODF to 0 (default: true).
%        contrast: this is a function handle that takes one input and gives
%            one output. It is applied to the signal reconstructed from the
%            SH coefficients passed. This option may be used, for example,
%            if the squared root of the ODF has been computed using 
%            signal2sqrtshodf within MiSFIT. In that case, @(x)(x.*x)
%            should be passed. Pass an empty array, [], to avoid using this
%            contrast function (default: []).
%        angleres: the angular resolution used to plot the glyphs. This is 
%            just an integer telling the minimum number of vertices to use.
%            In practice, either a icosahedron or an octahedron is
%            successively refined to obtain a number of vertices greater 
%            than or equal to angleres. In any case, angleres over 642 will
%            be cropped to angleres = 642 (default: 162)
%        glyphspc: 3x1 vector of integers. This is a subsampling factor
%            used to avoid too many glyphs to be represented. If glyphspc =
%            [s1,s2,s3], only one glyph each s1 voxels in the 1st diection,
%            s2 voxels in the 2nd direction, and s3 voxels in the 3rd
%            direction is actually plotted (default: [1;1;1])
%        glyphsc: 1x1 double, a constant scaling factor applied to all the
%            glyphs reprsented (default: 1)
%        glyphscvar: describes an additional scaling applied to the glyphs
%            depending on some scalar measurement derived from the 
%            diffusion data. Its value may be:
%              - 'none' (default), no aditional scaling is applied
%              - 'ga', the glyphs are scaled depending on a generalized
%                 anisotropy computed from the first input, SH, as:
%                             sqrt(sum(SH(:,:,:,2:end).^2,4)./sum(SH.^2,4))
%              - 'bg', the glyphs are scaled by the background image
%                 described by the 'bgimage' optional argument. NOTE: you
%                 can set the 'bgalpha' to 0 to just apply this scaling to
%                 the glyphs without actually plotting the image itself.
%                 NOTE2: you shouldn't use the 'bg' option if the selected
%                 background image was 'color', since in this case the
%                 scalar values are more or less random indices to a
%                 colormap.
%
%    NOTE: Each time the function is called, the colormap of the axes it
%    uses is modified (this is necessary in order to simultaneously
%    represent OFS and FAs or color-coded info, for example). In case you
%    add your own plots to these axes, be prepared to struggle with some
%    weird coloring effects for this reason (they can be worked around by
%    adding new rows to the color map and using 'CDataMapping', 'direct'
%    when plotting colormap-based graphics). 

[M,N,P,L,OK,msg] = plotdmri3d_check_mandatory_arguments(SH,i1,i2,i3);
if(~OK), error(msg); end

% Parse the optional input arguments:
%   Anatomy:
opt.origin = [0;0;0];      optchk.origin = [true,true];     % always 3x1 double
opt.direction = eye(3);    optchk.direction = [true,true];  % always 3x3 double
opt.space = 'RAS';         optchk.space = [true,false];     % string, variable length
opt.mframe = eye(3);       optchk.mframe = [true,true];     % always 3x3 double
% General plot:
opt.ha = [];               optchk.ha = [false,false];
opt.bbox = true;           optchk.bbox = [true,true];        % always 1x1 boolean
opt.rastxt = 20;           optchk.rastxt = [true,true];      % always 1x1 double
opt.mask = true(M,N,P);    optchk.mask = [true,true];        % always MxNxP boolean
% Background plot:
opt.bgsh = [];             optchk.bgsh = [true,false];       % double with variable size
opt.bgimage = 'none';      optchk.bgimage = [false,false];
opt.bgalpha = 1;           optchk.bgalpha = [true,true];     % always 1x1 double
% Glyphs plot:
opt.glyphs = true;         optchk.glyphs = [true,true];      % always 1x1 boolean
opt.clip = true;           optchk.clip = [true,true];        % always 1x1 boolean
opt.contrast = [];         optchk.contrast = [false,false];
opt.angleres = 162;        optchk.angleres = [true,true];    % always 1x1 double
opt.glyphspc = [1;1;1];    optchk.glyphspc = [true,true];    % always 3x1 double
opt.glyphsc = 1;           optchk.glyphsc = [true,true];     % always 1x1 double
opt.glyphscvar = 'none';   optchk.glyphscvar = [true,false]; % variable size string
% Automatically parse:
opt = custom_parse_inputs(opt,optchk,varargin{:});

space = plotdmri3d_space_string2matrix(opt.space);

if(isempty(opt.ha))
    hf = figure;
    ha = axes('Parent',hf);
else
    ha = opt.ha;
end
set(ha,'NextPlot','add');

if(~isempty(opt.contrast))
    if(~isa(opt.contrast,'function_handle'))
        error('''contrast'' can be either emptyor a function handle');
    end
end

if(opt.bbox)
    corners = [ 1,1,1; 2,1,1; 2,2,1; 1,2,1; 1,1,2; 2,1,2; 2,2,2; 1,2,2 ]'; %3x8
    idxbbox = [1,2,3,4,1,5,6,7,8,5,8,4,3,7,6,2]; %1x16
    % ------------
    limits  = [0,M-1;0,N-1;0,P-1];
    xyz     = [ limits(1,corners(1,idxbbox)); limits(2,corners(2,idxbbox)); limits(3,corners(3,idxbbox)) ]; %3x14
    bboxxyz = [opt.direction,opt.origin]*[xyz;ones(1,16)]; %3x16
    ho.bbox(1) = line(bboxxyz(1,:),bboxxyz(2,:),bboxxyz(3,:),'Color',[.0,.0,.0],'LineWidth',1,'LineStyle',':','Parent',ha);
    % ------------
    % ------------
    limits  = [ min(bboxxyz(1,:)),max(bboxxyz(1,:)); min(bboxxyz(2,:)),max(bboxxyz(2,:)); min(bboxxyz(3,:)),max(bboxxyz(3,:)) ];
    bboxras = [ limits(1,corners(1,idxbbox)); limits(2,corners(2,idxbbox)); limits(3,corners(3,idxbbox)) ]; %3x14
    ho.bbox(2) = line(bboxras(1,:),bboxras(2,:),bboxras(3,:),'Color',[.0,.0,.5],'LineWidth',1,'LineStyle',':','Parent',ha);
    % ------------
    if(opt.rastxt>0)
        ho.rastxt(1) = text(limits(1,space(1)+1),mean(limits(2,:)),mean(limits(3,:)),'R','FontSize',opt.rastxt,'Color',[.0,.0,.8],'FontWeight','bold','Parent',ha);
        ho.rastxt(2) = text(limits(1,2-space(1)),mean(limits(2,:)),mean(limits(3,:)),'L','FontSize',opt.rastxt,'Color',[.0,.0,.8],'FontWeight','bold','Parent',ha);
        ho.rastxt(3) = text(mean(limits(1,:)),limits(2,space(2)+1),mean(limits(3,:)),'A','FontSize',opt.rastxt,'Color',[.0,.0,.8],'FontWeight','bold','Parent',ha);
        ho.rastxt(4) = text(mean(limits(1,:)),limits(2,2-space(2)),mean(limits(3,:)),'P','FontSize',opt.rastxt,'Color',[.0,.0,.8],'FontWeight','bold','Parent',ha);
        ho.rastxt(5) = text(mean(limits(1,:)),mean(limits(2,:)),limits(3,space(3)+1),'S','FontSize',opt.rastxt,'Color',[.0,.0,.8],'FontWeight','bold','Parent',ha);
        ho.rastxt(6) = text(mean(limits(1,:)),mean(limits(2,:)),limits(3,2-space(3)),'I','FontSize',opt.rastxt,'Color',[.0,.0,.8],'FontWeight','bold','Parent',ha);
    end
end

% Check the consistency of the background image:
if(ischar(opt.bgimage))
    if(~strcmp(opt.bgimage,'none')) % opt.bgsh cannot be empty!
        [OK,msg] = plotdmri3d_check_bgsh_size(opt.bgsh,SH);
        assert(OK,msg);
        
    end
else % opt.bgimage must have the proper size
    assert(isequal(size(opt.bgimage),[M,N,P]),'The background image (bgimage) must be a string or a 3-D volume with size matching that of SH');
end
    
bgmap = gray(256);
if(ischar(opt.bgimage))
    if(strcmp(opt.bgimage,'none'))
        bgimage = [];
    elseif(strcmp(opt.bgimage,'fa'))
        bgimage = sqrt(sum(opt.bgsh(:,:,:,2:6).*opt.bgsh(:,:,:,2:6),4)./sum(opt.bgsh(:,:,:,1:6).*opt.bgsh(:,:,:,1:6),4));
    elseif(strcmp(opt.bgimage,'gfa'))
        bgimage = sqrt(sum(opt.bgsh(:,:,:,2:end).*opt.bgsh(:,:,:,2:end),4)./sum(opt.bgsh(:,:,:,1:end).*opt.bgsh(:,:,:,1:end),4));
    elseif(strcmp(opt.bgimage,'md'))
        bgimage = opt.bgsh(:,:,:,1);
    elseif(strcmp(opt.bgimage,'color'))
        bgimage = sqrt(sum(opt.bgsh(:,:,:,2:6).*opt.bgsh(:,:,:,2:6),4)./sum(opt.bgsh(:,:,:,1:6).*opt.bgsh(:,:,:,1:6),4));
        [grads,~,~] = icosamplesSphere(6,'O1',true,'verbose',false);
        signal = sh2signal( opt.bgsh(:,:,:,1:6), grads, 'mask', opt.mask );
        grads = ((opt.mframe)*(grads'))';
        [~,signal] = max(signal,[],4);
        rc = abs(grads(:,1));
        gc = abs(grads(:,2));
        bc = abs(grads(:,3));
        R  = rc(signal);
        G  = gc(signal);
        B  = bc(signal);
    end
else
    bgimage  = opt.bgimage;
end

if(~isempty(bgimage))
    bgimage(~opt.mask) = nan;
    idx = i1(1):i1(2);
    idy = i2(1):i2(2);
    idz = i3(1):i3(2);
    CS  = squeeze(bgimage(idx,idy,idz));
    if(strcmp(opt.bgimage,'color'))
        rgbmask = isnan(CS);
        RGB = zeros(size(CS,1),size(CS,2),3);
        R = squeeze(R(idx,idy,idz)).*CS; R(rgbmask) = 0;
        G = squeeze(G(idx,idy,idz)).*CS; G(rgbmask) = 0;
        B = squeeze(B(idx,idy,idz)).*CS; B(rgbmask) = 0;
        RGB(:,:,1) = R;
        RGB(:,:,2) = G;
        RGB(:,:,3) = B;
        [CS,bgmap] = rgb2ind(RGB,256);
        CS = double(CS);
        CS(rgbmask) = nan;
    end
    [idx,idy,idz] = meshgrid(idx,idy,idz);
    [MS,NS,PS] = size(idx);
    xyz = [idx(:)'-1;idy(:)'-1;idz(:)'-1];
    xyz = [opt.direction,opt.origin]*[xyz;ones(1,MS*NS*PS)];
    XS  = reshape(xyz(1,:),[MS,NS,PS]);
    YS  = reshape(xyz(2,:),[MS,NS,PS]);
    ZS  = reshape(xyz(3,:),[MS,NS,PS]);
    XS2 = squeeze(XS);
    YS2 = squeeze(YS);
    ZS2 = squeeze(ZS);
    if(PS==1)
        XS2 = XS';
        YS2 = YS';
        ZS2 = ZS';
    end
    cmap  = colormap(ha);
    nclrs = size(cmap,1);
    cmap  = [cmap;bgmap];
    alpha = ones(size(CS))*opt.bgalpha;
    alpha(isnan(CS)) = 0;
    if(strcmp(opt.bgimage,'color'))
        CS = double(CS) + nclrs + 1;
        ho.bgimage = surf(XS2,YS2,ZS2,CS,'LineStyle','none','FaceColor','texturemap','FaceAlpha','texturemap','AlphaData',alpha,'Parent',ha,'CDataMapping','direct');        
    else
        CS = 255*( CS - min(CS(~isnan(CS))) )/( max(CS(~isnan(CS))) - min(CS(~isnan(CS))) ) + nclrs + 1;
        ho.bgimage = surf(XS2,YS2,ZS2,CS,'LineStyle','none','FaceColor','interp','FaceAlpha','interp','AlphaData',alpha,'Parent',ha,'CDataMapping','direct');
    end
    colormap(ha,cmap);
end

if(opt.glyphs)
    % SCALING -------------------------------------------------------------
    scales = ones(M,N,P)*nthroot(prod(diag(opt.direction)),3)*(opt.glyphsc);
    switch(opt.glyphscvar)
        case 'none'
        case 'ga'
            scales = scales.*sqrt(sum(SH(:,:,:,2:end).*SH(:,:,:,2:end),4)./sum(SH.*SH,4));
        case 'bg'
            if(~isempty(bgimage))
                bgimage = (bgimage-min(bgimage(:)))./(max(bgimage(:))-min(bgimage(:)));
                scales  = scales.*bgimage;
            end
        otherwise
            error(['Wrong variable scaling option ''glyphscvar'': ',opt.glyphscvar]);
    end
    scales(~opt.mask) = 0;
    % END SCALING ---------------------------------------------------------
    pol   = plotdmri3d_sphere_sampling(opt.angleres); % Will be used to draw patches
    vts   = ((opt.mframe)*(pol.vertices'))'; % original vertices over the unit sphere
    NV    = size(vts,1);
    SHmat = GenerateSHMatrix( L, pol.vertices );
    x     = i1(1):opt.glyphspc(1):i1(2);
    y     = i2(1):opt.glyphspc(2):i2(2);
    z     = i3(1):opt.glyphspc(3):i3(2);
    mask  = opt.mask(x,y,z);
    [y,x,z] = meshgrid(y,x,z);
    x = x(mask);
    y = y(mask);
    z = z(mask);
    NGL = numel(x);
    xyz = [opt.direction,opt.origin]*[x'-1;y'-1;z'-1;ones(1,NGL)];
    ho.glyphs = zeros(1,NGL);
    for n=1:NGL % Draw glyphs one by one
        SHcoef = SH(x(n),y(n),z(n),:);
        SHcoef = SHcoef(:);
        SHcoef = SHmat*SHcoef;
        if(~isempty(opt.contrast))
            SHcoef = opt.contrast(SHcoef);
        end
        if(opt.clip)
            SHcoef(SHcoef<0) = 0;
        end
        SHcoef = scales(x(n),y(n),z(n))*SHcoef/max(SHcoef(:));
        vts2   = repmat(SHcoef,[1,3]).*vts; % Distances are weighted depending on the ODF
        vts2   = vts2 + repmat(xyz(:,n)',[NV,1]);
        ho.glyphs(n) = patch( ...
            'Vertices', vts2, ...
            'Faces',    pol.facets, ...
            'EdgeColor', [0,0,0], ...
            'EdgeAlpha', 0, ...
            'FaceVertexCData', abs(vts), ...
            'FaceColor', 'interp', ...
            'FaceLighting', 'phong' ...
            );
    end
end

% ---------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------
function [M,N,P,L,OK,msg] = plotdmri3d_check_mandatory_arguments(SH,i1,i2,i3)
[M,N,P,K] = size(SH);
OK  = true;
msg = '';
L   = (sqrt(8*(K-1)+9)-3)/2; % SH order
if( abs(L-round(L))>1.0e-9 || L<2 )
    msg = 'Weird size of the SH volume. Its fourth dimension should have size 6, 15, 28, 45, ..., (L+1)(L+2)/2, with L=2,4,6,...';
    OK  = false;
    return;
end
if(~isequal(size(i1),[1,2]))
    msg = 'i1 must have size 1x2';
    OK  = false;
    return;
end
if(~isequal(size(i2),[1,2]))
    msg = 'i2 must have size 1x2';
    OK  = false;
    return;
end
if(~isequal(size(i3),[1,2]))
    msg = 'i3 must have size 1x2';
    OK  = false;
    return;
end
i1 = round(i1);
i2 = round(i2);
i3 = round(i3);
d1 = i1(2)-i1(1);
d2 = i2(2)-i2(1);
d3 = i3(2)-i3(1);
if( (d1~=0) && (d2~=0) && (d3~=0) )
    msg = 'One of i1, i2, or i3 must represent a singleton dimension';
    OK  = false;
    return;
end
if( any(i1<1) || any(i1>M) || (i1(1)>i1(2)) )
    msg = 'Values of i1 out of range or non-increasing';
    OK  = false;
    return;
end
if( any(i2<1) || any(i2>N) || (i2(1)>i2(2)) )
    msg = 'Values of i2 out of range or non-increasing';
    OK  = false;
    return;
end
if( any(i3<1) || any(i3>P) || (i3(1)>i3(2)) )
    msg = 'Values of i3 out of range or non-increasing';
    OK  = false;
    return;
end

% ---------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------
function [OK,msg] = plotdmri3d_check_bgsh_size(bgsh,SH)
[M,N,P,~] = size(SH);
[M2,N2,P2,K] = size(bgsh);
msg = '';
OK  = true;
if(~isequal([M,N,P],[M2,N2,P2]))
    OK = false;
    msg = 'The size of the background image doesn''t match that of the SH volume';
    return;
end
L   = (sqrt(8*(K-1)+9)-3)/2; % SH order
if( abs(L-round(L))>1.0e-9 || L<2 )
    msg = 'Weird size of the backgound SH volume. Its fourth dimension should have size 6, 15, 28, 45, ..., (L+1)(L+2)/2, with L=2,4,6,...';
    OK  = false;
    return;
end

% ---------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------
function space = plotdmri3d_space_string2matrix(spacestr)
space = ones(3,1);
if(length(spacestr)>3)
    rest = spacestr;
    [current,rest] = strtok(rest,'-');
    switch(current)
        case 'right'
            spacestr = 'r';
        case 'left'
            spacestr = 'l';
        otherwise
            error('Wrong space string');
    end
    [current,rest] = strtok(rest,'-');
    switch(current)
        case 'anterior'
            spacestr(2) = 'a';
        case 'posterior'
            spacestr(2) = 'p';
        otherwise
            error('Wrong space string');
    end
    [current,~] = strtok(rest,'-');
    switch(current)
        case 'superior'
            spacestr(3) = 's';
        case 'inferior'
            spacestr(3) = 'i';
        otherwise
            error('Wrong space string');
    end
elseif(length(spacestr)<3)
    error('Wrong space string');
end
spacestr = lower(spacestr);
if(spacestr(1)=='r'),space(1)=1; elseif(spacestr(1)=='l'), space(1)=0; else error('Wrong space string'); end %#ok<SEPEX>
if(spacestr(2)=='a'),space(2)=1; elseif(spacestr(2)=='p'), space(2)=0; else error('Wrong space string'); end %#ok<SEPEX>
if(spacestr(3)=='s'),space(3)=1; elseif(spacestr(3)=='i'), space(3)=0; else error('Wrong space string'); end %#ok<SEPEX>

% ---------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------
function polyhedron = plotdmri3d_sphere_sampling(N)

if(N>642), N=642; end
nvts = [6,12,18,42,66,162,258,642,1026,2562];
src  = [1,2,1,2,1,2,1,2,1,2];
pos  = find(nvts>=N,1,'first');
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
    [vts,idx1,idx2] = plotdmri3d_remove_duplicated_grads(polyhedron.vertices);
    vts2 = optimizeGradients(vts,'verbose',false);
    vts(idx1,:) =  vts2;
    vts(idx2,:) = -vts2;
    polyhedron.vertices = vts;
end

% ---------------------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------------------
function [grads,idx,idx2] = plotdmri3d_remove_duplicated_grads(grads)
moduli = sqrt(sum(grads.*grads,2));
grads  = grads./repmat(moduli,[1,3]);
% Find dot-products of each gradient with the remaining ones:
dp   = grads(:,1)*(grads(:,1)') + grads(:,2)*(grads(:,2)') + grads(:,3)*(grads(:,3)');
dp   = abs(dp); % The closer to 1, the more similar
dp   = dp - diag(diag(dp)); % Avoid self-similarity
idx  = 1:size(grads,1);
idx2 = 1:(size(grads,1)/2);
for n=1:size(grads,1)/2
    pos = idx(n);
    [~,mpos] = max(dp(pos,:)); % mpos is the most similar gradient to n, so...
    idx = setdiff(idx,mpos);   % ... remove it from the list
    idx2(n) = mpos;            % store the most similar gradient to idx(n)
end
grads = grads(idx,:);
