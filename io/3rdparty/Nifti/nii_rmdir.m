function nii_rmdir(in,in2)
  if(nargin<2)
    in2 = 's';
  end
  [sf,sfname,vs] = check_software_platform;
  if(sf==2) % Octave
    confirm_recursive_rmdir( false, 'local' );
  end
  rmdir(in,in2);
end
