function E = costFunctionDualAbsoluteQuadric(H, X)
%% COSTFUNCTIONDUALABSOLUTEQUADRIC computes DAQ Cost function to find the Optimal Intrinsic Parameters

%   Input
%       H    - Homography Matrix of Plane at Infinity
%       X    - Approximate Values of Intrinsics (1 x 5)
%
%   Output
%       E    - Computed Cost

%% Function starts here

% Transform Intrinsics to Matrix Form
K = [X(1) X(2) X(3); X(4) X(5) X(6); X(7) X(8) X(9)];

% Compute the Quadric W^-1 (Image of the Absolute Quadric)
W_inv = K*K';

% Initialize Cost
E = [];

% Compute the Cost using DAQ's Equation
for i = 1:length(H)
    
    % Compute DAQ
    DAQ =   H{i} * W_inv * H{i}';
    
    % Compute Cost
    E1 = DAQ - W_inv;
    
    % Append Computed Cost
    E = [E E1(1,:) E1(2,:) E1(3,:)] ;
    
end

end
