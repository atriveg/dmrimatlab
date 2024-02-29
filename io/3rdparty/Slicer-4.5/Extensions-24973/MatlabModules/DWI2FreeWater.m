function outputParams=DWI2FreeWater(inputParams)


imgDWI = cli_imageread(inputParams.inputvolumeDWI);

tl = inputParams.tl;
tu = inputParams.tu;
L = inputParams.L;
lambda = inputParams.lambda;
fmin = inputParams.fmin;
fmax = inputParams.fmax;
fbins = inputParams.fbins;
fdepth = inputParams.fdepth;
O2 = logical(inputParams.O2);
fiters = inputParams.fiters;
fth = inputParams.fth;
chunksz = inputParams.chunksz;

addpath('/home/atriveg/Dropbox/Code/FreeWaterMatlab');
[dwi,gi,bi] = SlicerVol2DWI(imgDWI);

useMask = false;
if(isfield(inputParams,'inputvolumeMask'))
	if(~isempty(inputParams.inputvolumeMask))
		useMask = true;
	end
end

if(useMask)
	mask = cli_imageread(inputParams.inputvolumeMask);
	mask = logical(mask.pixelData);
	[f,SH] = dwi2freewater( dwi, gi, bi, ...
   		'tl', tl, 'tu', tu, ...
   		'L', L, 'lambda', lambda, ...
		'flogerr', false, ...
   		'fmin', fmin, 'fmax', fmax, 'fbins', fbins, 'fdepth', fdepth, ...
                'O2', O2, 'fiters', fiters, 'fth', fth, ...
   		'chunksz', chunksz, ...
		'verbose', true, 'mask', mask );
	DTI = SHADC2DTI( SH, 'chunksz', chunksz, 'mask', mask );
else
	[f,SH] = dwi2freewater( dwi, gi, bi, ...
   		'tl', tl, 'tu', tu, ...
   		'L', L, 'lambda', lambda, ...
		'flogerr', false, ...
		'fmin', fmin, 'fmax', fmax, 'fbins', fbins, 'fdepth', fdepth, ...
   		'O2', O2, 'fiters', fiters, 'fth', fth, ...
   		'chunksz', chunksz, ...
		'verbose', true );
	DTI = SHADC2DTI( SH, 'chunksz', chunksz );
end

fIMG.metaData.type = 'float';
fIMG.metaData.dimension = 3;
fIMG.metaData.space = imgDWI.metaData.space;
fIMG.metaData.sizes = sprintf('%d %d %d',size(f,1),size(f,2),size(f,3));
fIMG.metaData.space_directions = strrep(imgDWI.metaData.space_directions,'none ','');
fIMG.metaData.kinds = strrep(imgDWI.metaData.kinds,'list ','');
fIMG.metaData.endian = imgDWI.metaData.endian;
fIMG.metaData.encoding = 'raw';
fIMG.metaData.space_origin = imgDWI.metaData.space_origin;

fIMG.metaDataFieldNames.space_directions = imgDWI.metaDataFieldNames.space_directions;
fIMG.metaDataFieldNames.space_origin = imgDWI.metaDataFieldNames.space_origin;

fIMG.pixelData = single(f);

fIMG.ijkToLpsTransform = imgDWI.ijkToLpsTransform;

cli_imagewrite(inputParams.outputvolumeF, fIMG );

DTIIMG.metaData.type = 'float';
DTIIMG.metaData.dimension = 4;
DTIIMG.metaData.space = fIMG.metaData.space;
DTIIMG.metaData.sizes = sprintf('%d %d %d %d',9, size(DTI,1),size(DTI,2),size(DTI,3));
DTIIMG.metaData.space_directions = ['none ',fIMG.metaData.space_directions];
DTIIMG.metaData.kinds = ['3D-matrix ',fIMG.metaData.kinds];
DTIIMG.metaData.endian = fIMG.metaData.endian;
DTIIMG.metaData.encoding = 'raw';
DTIIMG.metaData.space_origin = fIMG.metaData.space_origin;
DTIIMG.metaData.measurement_frame = imgDWI.metaData.measurement_frame;

DTIIMG.metaDataFieldNames.space_directions = imgDWI.metaDataFieldNames.space_directions;
DTIIMG.metaDataFieldNames.space_origin = imgDWI.metaDataFieldNames.space_origin;
DTIIMG.metaDataFieldNames.measurement_frame = imgDWI.metaDataFieldNames.measurement_frame;

DTI = reshape(DTI,[size(DTI,1),size(DTI,2),size(DTI,3),9]);
DTI = permute(DTI,[4,1,2,3]);
DTIIMG.pixelData = single(DTI);

DTIIMG.ijkToLpsTransform = imgDWI.ijkToLpsTransform;

cli_imagewrite(inputParams.outputvolumeDTI, DTIIMG );

