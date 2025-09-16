function brdcst = is_broadcast_available_test
% function available = is_broadcast_available_test
%
%   Checks if the current matlab version suppots boadcast as in numpy
%
%   Suggestion: in the very beginning of your function, use:
%      >> is_broadcast_available = is_broadcast_available_test
%   within your function to avoid repeated calls to this function.
global is_broadcast_available_test_var; %#ok<GVMIS> 
if(isempty(is_broadcast_available_test_var))
    [sf,~,vs] = check_software_platform;
    if(sf==1) % This is Matlab
        % Broadcast is available only from release 2016b aka 9.1
        if(vs>9.05)
            is_broadcast_available_test_var = true;
        else
            is_broadcast_available_test_var = false;
        end
    elseif(sf==2) % This is Octave
        % Broadcast is available only from release 3.6.0
        if(vs>3.59999)
            is_broadcast_available_test_var = true;
        else
            is_broadcast_available_test_var = false;
        end
    else % This is unknown
        error('You''re running dmrimatlab from an unknown software platform');
    end
end
brdcst = is_broadcast_available_test_var;
