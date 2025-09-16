function out = phong
sf = check_software_platform;
if(sf==1)
    out = 'phong';
else
    out = 'gouraud';
end
end
