function compressedPixelData = nrrd_gzip_compress(pixelData,DataType)

cn = find(strcmp(DataType,{'double','single','logical','char','int8','uint8',...
    'int16','uint16','int32','uint32','int64','uint64'}));

if cn == 3 || cn == 4
    pixelData = uint8(pixelData);
    DataType  = 'uint8';
end

pixelData = cast(pixelData,DataType);

% Create a temporary folder to work in:
tmpDir = tempname;
try
    mkdir(tmpDir);
catch
    error('Unable to create tmp folder');
end

% Create a tmp file to write raw data:
rawFileName = sprintf('%s%snrrdRawData.bin',tmpDir,filesep);
fid = fopen(rawFileName,'w');
if(fid<0)
    error('Unable to write raw data to tmp file');
end

% Save raw data:
fwrite( fid, pixelData, DataType );
fclose(fid);

% gzip the raw data to a tmp file:
try
    gzipFileName = gzip(rawFileName);
catch
    error('Unable to gzip raw data from the tmp file');
end

% Read the compressed file to a byte stream:
fid = fopen(gzipFileName{1},'r');
if(fid<0)
    error('Unable to read compressed data from tmp file');
end
compressedPixelData = uint8(fread( fid, Inf, 'uint8' ));
fclose(fid);

% Clean-up:
sf = check_software_platform;
if(sf==2) % Octave
    confirm_recursive_rmdir( false, 'local' );
end
rmdir(tmpDir,'s');

end

