%function histogram(rsz,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',[.5,.0,.0],'LineWidth',2);
function hf = histogram(data,nbins,varargin);
if(nargin>1)
    if(ischar(nbins))
        varargin = { nbins, varargin{:} };
        nbins = max( round(numel(data)/1000), 10 );
    end
else
    nbins = max( round(numel(data)/1000), 10 );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.Normalization = 'pdf';   optchk.Normalization = [true,false];
opt.DisplayStyle = 'stairs'; optchk.DisplayStyle = [true,false];
opt.FaceColor = [.0,.0,.5];  optchk.EdgeColor = [true,true];
opt.FaceAlpha = 0.5;         optchk.FaceAlpha = [true,true];
opt.EdgeColor = [.0,.0,.0];  optchk.EdgeColor = [true,true];
opt.EdgeAlpha = 1.0;         optchk.EdgeAlpha = [true,true];
opt.LineWidth = 2;           optchk.LineWidth = [true,true];
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nn,xx] = hist( data(:), nbins );
nn = nn/sum(nn)/(xx(2)-xx(1));
hf = bar(xx,nn,'hist','edgealpha',opt.EdgeAlpha,'edgecolor',opt.EdgeColor,'facecolor',opt.FaceColor,'facealpha',opt.FaceAlpha);
end
