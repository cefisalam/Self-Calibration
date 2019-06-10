function M = skewSymmetricMatrix(V)
%% SKEWSYMMETRICMATRIX computes Skew Symmetric Matrix of the given Vector (To do Matrix Multiplication instead of Cross Product)

%   Input
%       V    - Input Vector (3 x 1)
%
%   Output
%       M    - Skew Symmetric Matrix (3 x 3)

%% Function starts here

% Transform Vector into Skew Symmetric Matrix 
M = [ 0   -V(3)   V(2)
     V(3)   0    -V(1) 
    -V(2)  V(1)    0  ];

end

