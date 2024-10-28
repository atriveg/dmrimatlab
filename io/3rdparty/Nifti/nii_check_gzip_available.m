function success = nii_check_gzip_available(file)

if(nargin<1)
    file = '';
end

success = true;

v = version;
p = strfind(v,'.');
if(~isempty(p))
    v = v(1:p(1)+1);
end

% For some reason, usejava('jvm') randomly fails in deployed
% executables, so we have to work around by avoiding the 
% check on the Java Virtual Machine:
if(isdeployed)
    withjava = true;
else
    withjava = usejava('jvm');
end

if(   (str2double(v)<7.1) || ~withjava   )
    fprintf(1,'   >> I won''t be able to gunzip the file %s:\n',file);
    if( str2double(v)<7.1 )
        fprintf(1,'       >> Required version is at least 7.1 but I found version %s\n',v);
    end
    if( ~withjava )
        fprintf(1,'       >> I require JVM to be running but usejava(''jvm'') returned %d\n\n',usejava('jvm'));
    end
    success = false;
else
    if(~isempty(file))
        fprintf(1,'I will gunzip file %s\n',file);
    end
end

end