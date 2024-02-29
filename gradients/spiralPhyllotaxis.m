function [G,ps] = spiralPhyllotaxis(Ns)
% function [G,ps] = spiralPhyllotaxis(Ns)
%
%   This function designs a set of N gradient directions evenly spaced in
%   (half) the unit sphere according to the botanics-inspired spiral
%   phyllotaxis described in:
%
%      D. Piccini, A. Littmann, S. Nielles-Vallespin, and M. O. Zenge
%      "Spiral Phyllotaxis: The Natural Way to Construct a 3D Radial
%         Trajectory in MRI"
%      Magnetic Resonance in Medicine 66: 1049-1056 (2011)
%
%   This design allows to create interleaved multi-shell samplings so that
%   each shell (roughly) uniformly covers the space of orientations without
%   repeating gradient directions from shell to shell. Contrarily to
%   icosahedron refinement-based schemes, Spiral Phyllotaxis allows
%   completely arbitrary numbers of gradient directions at each shell, so
%   that a flexible design is possible. On the other side of the coin, the
%   final stack of all gradients for all shells are not as uniformly
%   distributed as they become with the icosahedron strategy.
%
%   INPUTS:
%
%     Ns: A 1xS vector, with S the number of shells to design, each entry
%        being the number of gradients to cast into each shell, i.e.,
%        Ns = [n_1,n_2,...,n_S].
%
%   OUTPUTS:
%
%     G: An Nx3 matrix, with N = sum(Ns) = n_1+n2+...+n_S, each row being a
%        unit vector (i.e. a gradient direction).
%     ps: An Nx1 vector with values in the range 1...S, each position
%        indicating which shell the corresponding gradient belongs to.

%--------------------------------------------------------------------------
assert(ismatrix(Ns),'Ns must be a 1xS vector, with S the number of shells to design');
assert(size(Ns,1)<=1,'Ns must be a 1xS vector, with S the number of shells to design');
S = size(Ns,2); % The number of shells
N = sum(Ns);    % The total number of gradients
%--------------------------------------------------------------------------
gangle_deg = 180*(3-sqrt(5));    % The "golden angle" in degrees
gangle_rad = pi*gangle_deg/180;  % The "golden angle" in radians
phi       = ((1:N)')*gangle_rad; % The azimuth of the gradient directions
%--------------------------------------------------------------------------
% Design shell-by-shell
ps    = zeros(N,1);
G     = zeros(N,3);
pinit = 1;
for s=1:S
    % ------------------------------------
    % Put things in place:
    pend = pinit + Ns(s) - 1;
    ps(pinit:pend,1) = s;
    % ------------------------------------
    % Define the polar angle for this shell:
    theta = (pi/2)*sqrt( (1:Ns(s))'/Ns(s) );
    % ------------------------------------
    % Compute the gradients from the spherical coordinates:
    G(pinit:pend,1:3) = [ ...
        sin(theta).*cos(phi(pinit:pend,1)), ...
        sin(theta).*sin(phi(pinit:pend,1)), ...
        cos(theta)   ];
    % ------------------------------------
    % Prepare for the next shell:
    pinit = pend + 1;
    % ------------------------------------
end
%--------------------------------------------------------------------------
