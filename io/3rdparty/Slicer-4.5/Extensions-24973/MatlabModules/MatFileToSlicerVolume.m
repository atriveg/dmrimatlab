function MatFileToSlicerVolume(inputParams)

IMG   = load(inputParams.slicermatfile);
names = fieldnames(IMG);
img   = IMG.(names{1});
if(isfield(inputParams,'outputvolume'))
    cli_imagewrite(inputParams.outputvolume, img);
elseif(isfield(inputParams,'outputvolumeDWI'))
    cli_imagewrite(inputParams.outputvolumeDWI, img);
elseif(isfield(inputParams,'outputvolumeDTI'))
    cli_imagewrite(inputParams.outputvolumeDTI, img);
elseif(isfield(inputParams,'outputvolumeLabl'))
    cli_imagewrite(inputParams.outputvolumeLabel, img);
end
