function E = costFunctionKruppaSimplified(F, X)
%% COSTFUNCTIONKRUPPASIMPLIFIED computes Simplified Kruppa Cost function to find the Optimal Intrinsic Parameters

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

% Compute the Cost using Simplified Kruppa's Equation
for i = 1:size(F,3)-1
    for j = i+1:size(F,3)
        
        % Compute SVD of Fundamental Matrix
        [U,D,V] = svd(F(:,:,i,j)');
        
        % Seperate the Colummns of U & V
        u1 = U(:,1);
        u2 = U(:,2);
        u3 = U(:,3);
        
        v1 = V(:,1);
        v2 = V(:,2);
        v3 = V(:,3);
        
        % Singular Values r and s (3rd value, D(3,3) is 0)
        r = D(1,1);
        s = D(2,2);
        
        % Use the Simplified Kruppa's Equation
        A = (r^2 * v1' * W_inv * v1) * pinv(u2' * W_inv * u2);
        B = (r * s * v1' * W_inv * v2) * pinv(-u1' * W_inv * u2);
        C = (s^2 * v2' * W_inv * v2) * pinv(u1' * W_inv * u1);
        
        % Compute Cost
        E1 = A - B;
        E2 = B - C;
        E3 = C - A;
        
        % Append Computed Cost
        E = [E E1 E2 E3];
        
    end
end

end
