function brdcst = is_broadcast_available_test
% function available = is_broadcast_available_test
%
%   Checks if the current matlab version suppots boadcast as in numpy
%
%   Suggestion: in the very beginning of your function, use:
%      >> is_broadcast_available = is_broadcast_available_test
%   within your function to avoid repeated calls to this function.
global is_broadcast_available_test_var;
if(isempty(is_broadcast_available_test_var))
    % Broadcast is available only from release 2016b
    vers   = version('-release');
    versn  = round(str2double(vers(1:4)));
    versl  = vers(5);
    is_broadcast_available_test_var = true;
    if(versn<2016)
        is_broadcast_available_test_var = false;
    elseif(versn==2016)
        if(versl~='b')
            is_broadcast_available_test_var = false;
        end
    end
end
brdcst = is_broadcast_available_test_var;
