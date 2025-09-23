function dmrimatlab_demo

sf = check_software_platform;

if(sf==1)
    hf      = figure( 'Position', [600, 600, 325, 375] );
    btnsize = [10 0 300 50];
    FSZ     = 12;
    FSZ2    = FSZ+4;
else
    check_jupyter_available;
    hf      = figure( 'Position', [600, 600, 325, 400] );
    btnsize = [10 0 300 50];
    FSZ     = 18;
    FSZ2    = 22;
end
hp = uipanel ('title', 'Launch demo', 'position', [0 0 1 1], ...
                'parent', hf, 'FontSize', FSZ2, 'FontAngle', 'italic' );

gp = uibuttongroup (hp, 'Position', [ 0 0 1 1]);
% Create a buttons in the group
if(sf==2)
    b01 = uicontrol(gp, 'style', 'text', ...
                    'string', 'EXE:', ...
                    'Position', [ 10 300 50 50 ], ...
                    'FontSize', FSZ, 'enable', 'on' );
    b02 = uicontrol(gp, 'style', 'edit', ...
                    'string', 'NONE', ...
                    'Position', [ 60 300 250 50 ], ...
                    'FontSize', FSZ-2, 'enable', 'on', ...
                    'HorizontalAlignment', 'left', ...
                    'callback', {@check_demo_config_file, false} );
    check_demo_config_file(b02,[],true);
else
    b02 = [];
end
b1 = uicontrol(gp, 'style', 'pushbutton', ...
                'string', 'A tour on DMRIMatlab', ...
                'Position', btnsize + [0 250 0 0], ...
                'FontSize', FSZ, 'enable', 'on', ...
                'callback', {@dmrimatlab_demo_callback, 1, b02} );
b2 = uicontrol(gp, 'style', 'pushbutton', ...
                'string', 'Demo DT-MRI', ...
                'Position', btnsize + [0 200 0 0], ...
                'FontSize', FSZ , ...
                'callback', {@dmrimatlab_demo_callback, 2, b02} );
b3 = uicontrol(gp, 'style', 'pushbutton', ...
                'string', 'Demo AMURA', ...
                'Position', btnsize + [0 150 0 0], ...
                'FontSize', FSZ , ...
                'callback', {@dmrimatlab_demo_callback, 3, b02} );
b4 = uicontrol(gp, 'style', 'pushbutton', ...
                'string', 'Demo MiSFIT', ...
                'Position', btnsize + [0 100 0 0], ...
                'FontSize', FSZ , ...
                'callback', {@dmrimatlab_demo_callback, 4, b02} );
b5 = uicontrol(gp, 'style', 'pushbutton', ...
                'string', 'Demo MAP-MRI/MAPL', ...
                'Position', btnsize + [0 50 0 0], ...
                'FontSize', FSZ , ...
                'callback', {@dmrimatlab_demo_callback, 5, b02} );
b6 = uicontrol(gp, 'style', 'pushbutton', ...
                'string', 'Demo HYDI-DSI', ...
                'Position', btnsize, ...
                'FontSize', FSZ , ...
                'callback', {@dmrimatlab_demo_callback, 6, b02} );
end

% ----------------------------------------------------------------------
function check_jupyter_available
% Check if jupyter-lab is available
[status, output] = system ('jupyter-lab --version');
if(status==0)
    [status,~] = system('jupyter kernelspec list | grep ''^\s*octave''');
    if(status==0)
        fprintf(1,'Demos will run over jupyter-lab version %s\n',output);
    else
        fprintf(1,'It seems that Jupyter''s kernel for Octave is not available, Please install it\n');
        error('Jupyter''s kernel for Octave is not available');
    end
else
    fprintf(1,'It seems that jupyter-lab is not available\n');
    fprintf(1,'Please install jupyter-lab and Jupyter''s kernel for Octave\n');
    error('jupyter-lab is not available');
end
end

% ----------------------------------------------------------------------
function retval = check_demo_config_file(h,evt,flag)
% -------------------------------------
retval      = 0;
path0       = mfilename('fullpath');
[path0,~,~] = fileparts(path0);
config      = sprintf('%s/config.octave',path0);
% -------------------------------------
if(exist(config,'file')~=2)
    if(flag)
        octave_exe = '/usr/bin/octave';
    else
        octave_exe = get(h,'string');
    end
    fprintf(1,'Config file not found. Creating: <%s>\n',config);
    fid = fopen(config,'w');
    fprintf(fid,'OCTAVE_EXE=%s\n',octave_exe);
    fclose(fid);
end
% -------------------------------------
if(flag)
    octave_exe = '/usr/bin/octave';
    fid = fopen(config,'r');
    str = fgetl(fid);
    while(str~=-1)
        str = strtrim(str);
        % --------------
        idx = strfind(str,'=');
        if(length(idx)~=1)
            fclose(fid);
            error( sprintf('Bad line in config file: <%s>. Must be PROPERTY=value',str) );
        end
        if(idx==1)
            fclose(fid);
            error( sprintf('Bad line in config file: <%s>. Must be PROPERTY=value',str) );
        end
        if(idx==length(str))
            fclose(fid);
            error( sprintf('Bad line in config file: <%s>. Must be PROPERTY=value',str) );
        end
        % --------------
        property   = strtrim(str(1:idx-1));
        octave_exe = strtrim(str(idx+1:end));
        % --------------
        if(~strcmp(property,'OCTAVE_EXE'))
            fclose(fid);
            error( sprintf('Bad property in config file: %s',property) );
        end
        % --------------
        str = fgetl(fid);
    end
    fclose(fid);
else
    octave_exe = get(h,'string');
    fid = fopen(config,'w');
    fprintf(fid,'OCTAVE_EXE=%s\n',octave_exe);
    fclose(fid);
end
% -------------------------------------
set(h,'string',octave_exe);
end

% ----------------------------------------------------------------------
function retval = dmrimatlab_demo_callback(h,evt,idx,b0)
switch(idx)
    case 1,
        filename = 'tour_00index';
    case 2,
        filename = 'demo_DTMRI';
    case 3,
        filename = 'demo_AMURA';
    case 4,
        filename = 'demo_MiSFIT';
    case 5,
        filename = 'demo_MAPL';
    case 6,
        filename = 'demo_HYDI_DSI';
    otherwise
        error('Wrong callback value');
end
if(~isempty(b0))
    octave_exe = get(b0,'string');
else
    octave_exe =[];
end
retval = dmrimatlab_run_demo(octave_exe,filename);
end

% ----------------------------------------------------------------------
function retval = dmrimatlab_run_demo(octave_exe,filename)
path0       = mfilename('fullpath');
[path0,~,~] = fileparts(path0);
sf          = check_software_platform;
if(sf==1)
    open(sprintf('%s.mlx',filename));
    retval = 0;
else
    fprintf(1,'Running demo file: %s/%s.ipynb\n',path0,filename)
    comm = sprintf('export OCTAVE_EXECUTABLE=%s; jupyter-lab %s/%s.ipynb > /dev/null 2> /dev/null &',octave_exe,path0,filename);
    [retval, output] = system (comm);
    if(retval~=0)
        warning(sprintf('System call returned an error code %d:\n',retval));
        warning(sprintf('%s',output));
    end
end
end

