function opt = custom_parse_inputs(opt,optchk,varargin)
% function opt = custom_parse_inputs(opt,optchk,varargin)
%
%    Function to systematically parse option/value-like optional parameters
%    within an arbitrary function.
%
%    USAGE:
%
%    - Say your function has 2 mandatory arguments marg1 and marg2
%      and an undetermined number of optional arguments arg1,
%      arg2, ..., argN with values val1, val2, ... valN. Then your function
%      header should look like:
%        function [...] = myFun( marg1, marg2, varargin ); 
%      to be called something like:
%        [...] = myFun(mval1,mval2,'arg1',val1,'arg2',val2,...);
%    - In your function, create a structure opt whose fields are the
%      optional parameters you want to be parsed (with their default
%      values):
%        opt.arg1 = val1;
%        opt.arg2 = val2;
%         ...
%        opt.argN = valN;
%      Each field may have arbitrary data types and sizes.
%    - The function can check the data type and/or the size of the optional
%      values passed compared to those of the default values in the opt
%      structure. The structure optchk (second argument to the function) is
%      used for that (pass just an empty structure [] to avoid any
%      cheking):
%          + For each argi whose type and/or size you want to check, just
%            add a corresponding field to optchk and assign a two-component
%            vector of booleans:
%                optchk.argi = [true,false]
%            will check if the value passed is the same class as the 
%            default value originally placed in opt.argi, otherwise it will
%            throw a warning message.
%          +     optchk.argi = [false,true]
%            will check if the value passed is the same size as the
%            original default value placed in opt.argi, otherwise it will
%            throw a warning message.
%          +     optchk.argi = [true,true]
%            will check for the data type of the value passed. ONLY IN THE
%            CASE it is the same type as the default value, the size of the
%            value is further checked.
%          +     optchk.argi = [false,false]
%            is exactly the same as omitting the field argi from the
%            structure optchk.
%    - Finally, just call this function to parse the optional arguments:
%        opt = custom_parse_inputs(opt,optsizes,varargin{:});


st = dbstack; % Check what function is calling this one
if(length(st)>1)
    calling_fcn = upper(st(2).name);
else % Called from command window
    calling_fcn = [];
end

N = length(varargin);
if(rem(N,2)>0.5)
    error([calling_fcn,': Optional arguments must be passed as ''name'',value pairs']);
end

% Parse arguments pair by pair:
for n=1:N/2
    argument = varargin{2*n-1};
    value    = varargin{2*n};
    if(~ischar(argument))
        error([calling_fcn,': The name of optional arguments must be a character string']);
    end
    if(~isfield(opt,argument))
        error([calling_fcn,': argument <',argument,'> is not recognized as a valid optional argument. Optional arguments are case-sensitive']);
    end
    % Optional consistency checking
    if(~isempty(optchk)) % Checking required for some arguments
        if(isfield(optchk,argument)) % Checking required for this particular argument
            chcks = logical(optchk.(argument));
            % Make sure the user has introduced a proper structure field:
            if(numel(chcks)<1)
                chcks = [false,false];
            elseif(numel(chcks)<2)
                chcks(2) = false;
            elseif((numel(chcks)>2))
                chcks = chcks(1:2);
            end
            if(chcks(1)) % Need to check the data type
                old_class = class(opt.(argument));
                if(~isa(value,old_class)) % Wrong data type; size will not be checked
                    warning([calling_fcn,': The default data type for argument <',argument,'> is <',old_class,'> but a <',class(value),'> was passed']);
                elseif(chcks(2)) % Appropriate data type; size may need to be checked
                    old_size = size(opt.(argument));
                    new_size = size(value);
                    if(~isequal(old_size,new_size))
                        warning([calling_fcn,': The default array size for argument <',argument,'> is [',num2str(old_size),'] but a value with size [',num2str(new_size),'] was passed']);
                    end
                end
            else % No need to check the data type, but perhaps we need to check the size
                if(chcks(2))
                    old_size = size(opt.(argument));
                    new_size = size(value);
                    if(~isequal(old_size,new_size))
                        warning([calling_fcn,': The default array size for argument <',argument,'> is [',num2str(old_size),'] but a value with size [',num2str(new_size),'] was passed']);
                    end
                end
            end
        end
    end
    opt.(argument) = value;
end