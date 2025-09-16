%function histogram(rsz,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',[.5,.0,.0],'LineWidth',2);
function hf = histogram2(data1,data2,nbins,varargin);
if(nargin>2)
    if(ischar(nbins))
        varargin = { nbins, varargin{:} };
        nbins = max( round(sqrt(numel(data1))/100), 10 );
        nbins = [nbins,nbins];
    end
else
    nbins = max( round(sqrt(numel(data1))/100), 10 );
    nbins = [nbins,nbins];
end
if(length(nbins)==1)
    nbins = [nbins,nbins];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.Normalization = 'pdf';   optchk.Normalization = [true,false];
opt.DisplayStyle = 'stairs'; optchk.DisplayStyle = [true,false];
opt.EdgeColor = [.5,.0,.0];  optchk.EdgeColor = [true,true];
opt.LineWidth = 2;           optchk.LineWidth = [true,true];
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmri_hist3( [data1(:),data2(:)], nbins );
end
