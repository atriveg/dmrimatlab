function out = nii_findstr(in1,in2)
out = ~isempty(strfind(in2,in1));
end
