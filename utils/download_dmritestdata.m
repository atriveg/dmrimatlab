% This script automatically downloads all test data used throughout the
% examples and demos of the toolbox. It only needs to be run once. The data
% pieces are considerably large, so the download might take long.
%
%    IMPORTANT NOTES:
%
% Test data in this folder proceed from publicly available databases, and 
% should be used under the terms and conditions of their owners.
%
% * test_data.mat is an excerpt from volume HCP MGH-1007 of the Human 
%   Connectome Project (HCP) (http://www.humanconnectomeproject.org, 
%   https://db.humanconnectome.org). HCP: Principal Investigators: Bruce 
%   Rosen, M.D., Ph.D., Martinos Center at Massachusetts General Hospital, 
%   Arthur W. Toga, Ph.D., University of Southern California, Van J. 
%   Weeden, MD, Martinos Center at Massachusetts General Hospital). HCP 
%   funding was provided by the National Institute of Dental and 
%   Craniofacial Research (NIDCR), the National Institute of Mental Health 
%   (NIMH), and the National Institute of Neurological Disorders and Stroke
%   (NINDS). HCP is the result of efforts of co-investigators from the 
%   University of Southern California, Martinos Center for Biomedical 
%   Imaging at Massachusetts General Hospital (MGH), Washington University,
%   and the University of Minnesota. HCP data are disseminated by the 
%   Laboratory of Neuro Imaging at the University of Southern California.
% 
%   In case you use dataset test_data.mat for your purposes, please
%   consider citing:
% 
%       Fan, Q., Witzel, T., Nummenmaa, A., Van Dijk, K. R., Van Horn, J. 
%       D., Drews, M. K., Somerville, L. H., Sheridan, M. A., Santillana,
%       R. M., Snyder, J., Hedden, T., Shaw, E. E., Hollinshead, M. O., 
%       Renvall, V., Zanzonico, R., Keil, B., Cauley, S., Polimeni, J. R.,
%       Tisdall, D., Buckner, R. L., Wedeen, V. J., Wald, L. L., Toga, A.
%       W., Rosen, B. R., 2016. MGH–USC Human Connectome Project datasets 
%       with ultra-high b-value diffusion MRI. NeuroImage 124, 1108–1114.
% 
% * test_data2.mat is an excerpt from a dataset of the Cross-scanner and 
%   cross-protocol diffusion MRI data harmonisation database by the Cardiff
%   University - Brain  Research  Imaging  Centre  (CUBRIC): 
%   https://www.cardiff.ac.uk/cardiff-university-brain-research-imaging-
%         centre/research/projects/cross-scanner-and-cross-protocol-
%         diffusion-MRI-data-harmonisation.
%   The data were acquired at the UK National Facility for In Vivo MR
%   Imaging of Human Tissue Microstructure funded by the EPSRC (grant 
%   EP/M029778/1), and The Wolfson Foundation.
% 
%   In case you use dataset test_data.mat for your purposes, please
%   consider citing:
% 
%       Tax, C. M., Grussu, F., Kaden, E., Ning, L., Rudrapatna, U., Evans,
%       C. J., St-Jean, S., Leemans, A., Koppers,S., Merhof, D., et al., 
%       2019. Cross-scanner and cross-protocol diffusion MRI data 
%       harmonisation: A benchmark database and evaluation of algorithms. 
%       NeuroImage 195, 285–299.
% 
% * test_data3.mat is an excerpt from volume HCP WuMinn-139839 of the Human
%   Connectome Project (HCP) (http://www.humanconnectomeproject.org,
%   https://db.humanconnectome.org). HCP: Principal Investigators: Bruce 
%   Rosen, M.D., Ph.D., Martinos Center at Massachusetts General Hospital, 
%   Arthur W. Toga, Ph.D., University of Southern California, Van J. 
%   Weeden, MD, Martinos Center at Massachusetts General Hospital). HCP 
%   funding was provided by the National Institute of Dental and 
%   Craniofacial Research (NIDCR), the National Institute of Mental Health 
%   (NIMH), and the National Institute of Neurological Disorders and Stroke
%   (NINDS). HCP is the result of efforts of co-investigators from the
%   University of Southern California, Martinos Center for Biomedical
%   Imaging at Massachusetts General Hospital (MGH), Washington University,
%   and the University of Minnesota. HCP data are disseminated by the
%   Laboratory of Neuro Imaging at the University of Southern California.
% 
%   In case you use dataset test_data.mat for your purposes, please
%   consider citing:
% 
%       Van Essen, D. C., Smith, S. M., Barch, D. M., Behrens, T. E.,
%       Yacoub, E., Ugurbil, K., 2013. The WU-Minn Human Connectome
%       Project: An overview. NeuroImage 80, 62–79.
%
% * test_data4.mat is an already processed volume (subejct 1, session 1)
%   taken from the Cross-scanner and cross-protocol diffusion MRI data
%   harmonisation database by the  Cardiff  University - Brain  Research  
%   Imaging  Centre  (CUBRIC): https://www.cardiff.ac.uk/cardiff-
%     university-brain-research-imaging-centre/research/projects/
%     cross-scanner-and-cross-protocol-diffusion-MRI-data-harmonisation.
%   This means you won't find the raw DWI data set, but signal models 
%   already estimated from it, namely: the diffusion tensor volume and the
%   lamda parallel/perpendicular and SH ODF coefficientes computed with 
%   MiSFIT. The original, raw data were acquired at the UK National 
%   Facility for In Vivo MR Imaging of Human Tissue Microstructure funded 
%   by the EPSRC (grant EP/M029778/1), and The Wolfson Foundation.
%
%   In case you use dataset test_data.mat for your purposes, please
%   consider citing:
% 
%       Tax, C. M., Grussu, F., Kaden, E., Ning, L., Rudrapatna, U., Evans,
%       C. J., St-Jean, S., Leemans, A., Koppers,S., Merhof, D., et al., 
%       2019. Cross-scanner and cross-protocol diffusion MRI data 
%       harmonisation: A benchmark database and evaluation of algorithms. 
%       NeuroImage 195, 285–299.

help('download_dmritestdata');

oldpath = pwd;

cd(fileparts(which('setup__DMRIMatlab_toolbox')));
if(isfolder('./data'))
    fprintf(1,'Data folder already exists. Skip creation\n');
    cd('./data');
else
    [success,message,~] = mkdir('./data');
    if(~success)
        fprintf(1,'Could not create data folder: \n');
        fprintf(1,'   %s: \n',message)
    else
        cd('./data');
    end
end

BASEURL = 'http://www.lpi.tel.uva.es/~atriveg/';

DFILE = 'README.txt';
try
    fprintf(1,'[DOWNLOADING...] %s\n',[BASEURL,DFILE]);
    fname = websave(DFILE,[BASEURL,DFILE]);
    fname = strip_fname(fname);
    if(~strcmp(fname,DFILE))
        try
            delete(fname);
        catch
        end
        cd(oldpath);
        error('Something went wrong. Couldn''t download test data [FATAL]');
    end
catch
    try
        delete([DFILE,'.html']);
    catch
    end
    cd(oldpath);
    error('Something went wrong. Couldn''t download test data [FATAL]');
end

DFILE = 'test_data.mat';
try
    fprintf(1,'[DOWNLOADING...] %s\n',[BASEURL,DFILE]);
    fname = websave(DFILE,[BASEURL,DFILE]);
    fname = strip_fname(fname);
    if(~strcmp(fname,DFILE))
        try
            delete(fname);
        catch
        end
        fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
    end
catch
    try
        delete([DFILE,'.html']);
    catch
    end
    fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
end

DFILE = 'test_data2.mat';
try
    fprintf(1,'[DOWNLOADING...] %s\n',[BASEURL,DFILE]);
    fname = websave(DFILE,[BASEURL,DFILE]);
    fname = strip_fname(fname);
    if(~strcmp(fname,DFILE))
        try
            delete(fname);
        catch
        end
        fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
    end
catch
    try
        delete([DFILE,'.html']);
    catch
    end
    fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
end

DFILE = 'test_data3.mat';
try
    fprintf(1,'[DOWNLOADING...] %s\n',[BASEURL,DFILE]);
    fname = websave(DFILE,[BASEURL,DFILE]);
    fname = strip_fname(fname);
    if(~strcmp(fname,DFILE))
        try
            delete(fname);
        catch
        end
        fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
    end
catch
    try
        delete([DFILE,'.html']);
    catch
    end
    fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
end

DFILE = 'test_data4.mat';
try
    fprintf(1,'[DOWNLOADING...] %s\n',[BASEURL,DFILE]);
    fname = websave(DFILE,[BASEURL,DFILE]);
    fname = strip_fname(fname);
    if(~strcmp(fname,DFILE))
        try
            delete(fname);
        catch
        end
        fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
    end
catch
    try
        delete([DFILE,'.html']);
    catch
    end
    fprintf(1,'Something went wrong. Couldn''t download: %s [SKIP]\n',DFILE);
end

addpath(pwd);
cd(oldpath);

function fname = strip_fname(fname)
[~,ffile,fext] = fileparts(fname);
fname = [ffile,fext];
end