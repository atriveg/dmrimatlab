function nii_rmdir(in,in2)
  [sf,sfname,vs] = check_software_platform;
  if(sf==2) % Octave
    confirm_recursive_rmdir( false, 'local' );
    rmdir(in,in2);
  else % Matlab
    rmdir(in,in2);
  end
end
