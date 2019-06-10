function NormPlane = NormaltoPlaneAtInfinity(P, K)
%% NORMALTOPLANEATINFINITY computes the Normal to the Plane at Infinity

%   Input
%       P    - Projection Matrix of the Camra
%       K    - Approximate Values of Intrinsics (3 x 3)
%
%   Output
%       NormPlane  - Normal to the Plane at Infinity

%% Function starts here

% Compute the Quadric W^-1 (Image of the Absolute Quadric)
W_inv = K * K';

% Construct Symbolic Variables to Solve for Plane at Infinity
X = sym('X', 'real');
Y = sym('Y', 'real');
Z = sym('Z', 'real');
XX = sym('XX', 'real');

% Define Normal to Plane at Infinity
N = [X; Y; Z]; 

% Equation of Dual Absolute Quadric (DAQ)
Q = [W_inv, (W_inv * N); (N' * W_inv), (N' * W_inv * N)]; 

% Select Any One Projection Matrix
M = P(:, :, 5);

% DAQ AutoCalibration Equation
Calib = M * Q * M';

% Solve for [X, Y, Z] considering the System of Linear Equations
S = solve(Calib(1, 1) == (XX * W_inv(1, 1)), ...
          Calib(2, 2) == (XX * W_inv(2, 2)), ...
          Calib(1, 3) == (XX * W_inv(1, 3)), ...
          Calib(2, 3) == (XX * W_inv(2, 3)));


% Take only N = [X; Y; Z], Normal to Plane at Infinity
NormPlane = [double(S.X(1)); double(S.Y(1)); double(S.Z(1))];

end