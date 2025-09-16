function use = use_parallel_test
% function use = use_parallel_test
%
%   Checks if the library will use parallel computing
%
%   Suggestion: in the very beginning of your function, use:
%      >> use_parallel = use_parallel_test
%   within your function to avoid repeated calls to this function
global use_parallel_test_var; %#ok<GVMIS> 
if(isempty(use_parallel_test_var))
    use = false;
else
    use = use_parallel_test_var;
end
