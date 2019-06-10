function H = homographyPlaneAtInfinity(F, NormPlane)
%% HOMOGRAPHYPLANEATINFINITY computes the Inter-Image Homography (Hij) of Plane at Infinity

%   Input
%       F         - Fundamental Matrix
%       NormPlane - Normal to the Plane at Infinity
%
%   Output
%       H    - Homography Matrix of Plane at Infinity

%% Function starts here

count = 1;

for i = 1:size(F,3)-1
    for j = i+1:size(F,3)
        
        % Compute SVD of Fundamental Matrix
        [~,~,V] = svd(F(:,:,i,j)');
        
        % Cosider the Last Column Vector as Epipole
        V = V(:,end);
        
        % Compute Homography {H = [e21] * F12 + e21 * n'}
        H{count} =  skewSymmetricMatrix(V) * F(:,:,i,j) + V * NormPlane';
        count = count+1;
    end
end

end

