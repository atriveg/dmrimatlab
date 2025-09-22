function dmrimatlab_demo

sf = check_software_platform;
if(sf==2)
    % Check if jupyter-lab is available
    [status, output] = system ('jupyter-lab --version');
    if(status==0)
        fprintf(1,'Demos will run over jupyter-lab version %s\n',output);
    else
        fprintf(1,'It seems that jupyter-lab is not available\n');
        fprintf(1,'Please install jupyter-lab and Jupyter''s kernel for Octave\n');
        error('jupyter-lab is not available');
    end
end

hf = figure( 'Position', [600, 600, 275, 350] );
hp = uipanel ('title', 'Launch demo', 'position', [0 0 1 1], ...
                'parent', hf, 'FontSize', 22, 'FontAngle', 'italic' );

gp = uibuttongroup (hp, 'Position', [ 0 0 1 1]);
% Create a buttons in the group
b1 = uicontrol (gp, 'style', 'pushbutton', ...
                'string', 'A tour on dmrimatlab', ...
                'Position', [ 10 250 250 50 ], ...
                'FontSize', 16, 'enable', 'on', ...
                'callback', {@dmrimatlab_demo_callback, 1} );
b2 = uicontrol (gp, 'style', 'pushbutton', ...
                'string', 'Demo DT-MRI', ...
                'Position', [ 10 200 250 50 ], ...
                'FontSize', 16 , ...
                'callback', {@dmrimatlab_demo_callback, 2} );
b3 = uicontrol (gp, 'style', 'pushbutton', ...
                'string', 'Demo AMURA', ...
                'Position', [ 10 150 250 50 ], ...
                'FontSize', 16 , ...
                'callback', {@dmrimatlab_demo_callback, 3} );
b4 = uicontrol (gp, 'style', 'pushbutton', ...
                'string', 'Demo MiSFIT', ...
                'Position', [ 10 100 250 50 ], ...
                'FontSize', 16 , ...
                'callback', {@dmrimatlab_demo_callback, 4} );
b5 = uicontrol (gp, 'style', 'pushbutton', ...
                'string', 'Demo MAP-MRI/MAPL', ...
                'Position', [ 10 050 250 50 ], ...
                'FontSize', 16 , ...
                'callback', {@dmrimatlab_demo_callback, 5} );
b6 = uicontrol (gp, 'style', 'pushbutton', ...
                'string', 'Demo HYDI-DSI', ...
                'Position', [ 10 000 250 50 ], ...
                'FontSize', 16 , ...
                'callback', {@dmrimatlab_demo_callback, 6} );
end

% ----------------------------------------------------------------------
function retval = dmrimatlab_demo_callback(h,evt,idx)
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
retval = dmrimatlab_run_demo(filename);
end

% ----------------------------------------------------------------------
function retval = dmrimatlab_run_demo(filename)
path0       = mfilename('fullpath');
[path0,~,~] = fileparts(path0);
sf          = check_software_platform;
if(sf==1)
    open(sprintf('%s.mlx',filename));
    retval = 0;
else
    fprintf(1,'Running demo file: %s/%s.ipynb\n',path0,filename)
    comm = sprintf('export OCTAVE_EXECUTABLE=/usr/bin/octave; jupyter-lab %s/%s.ipynb > /dev/null 2> /dev/null &',path0,filename);
    [retval, output] = system (comm);
    if(retval~=0)
        warning(sprintf('System call returned an error code %d:\n',retval));
        warning(sprintf('%s',output));
    end
end
end

