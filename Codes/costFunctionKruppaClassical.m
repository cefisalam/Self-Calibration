function E = costFunctionKruppaClassical(F, X)
%% COSTFUNCTIONKRUPPACLASSICAL computes Classical Kruppa Cost function to find the Optimal Intrinsic Parameters

%   Input
%       F    - Fundamental Matrix between given two Images
%       X    - Approximate Values of Intrinsics (1 x 5)
%
%   Output
%       E    - Computed Cost

%% Function starts here

% Transform Intrinsics to Matrix Form
K = [X(1) X(2) X(3); 0 X(4) X(5); 0 0 1];

% Compute the Conic W^-1 (Image of the Absolute Conic)
W_inv = K*K';

% Initialize Cost
E = [];

% Compute the Cost using Classical Kruppa's Equation
for i = 1:size(F,3)-1
    for j = i+1:size(F,3)
        
        % 1st term of Kruppa's Equation
        A = F(:,:,i,j) * W_inv * F(:,:,i,j)';
        
        A = A/norm(A,'fro'); % Compute Frobenius Norm
        
        % 2nd term of Kruppa's Equation
        [~,~,V] = svd(F(:,:,i,j)'); % Compute SVD of Fundamental Matrix
        
        V = V(:,end); % Last Colum Vector is the Epipole
        
        Epi = skewSymmetricMatrix(V);
        
        B = Epi * W_inv * Epi';
        
        B = B/norm(B,'fro'); % Compute Frobenius Norm
        
        % Compute Cost
        E1 = A - B;
        
        % Append Computed Cost
        E = [E E1(1,1:3) E1(2,2:3)];
        
    end
end

end
