clear, clc;

format long;

% Load the Given Data
load('data.mat');

% Approximate Intrinsics (Given)
K_approx = A;
disp('Approximate Intrinsics (Given)')
disp(K_approx)

% Intrinsics in Vector Form
X0 = [A(1,:) A(2,2:3)];

%% Self Calibration using Mendonca & Cipolla's Method to Find Intrinsics

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolX', 1e-10);

% Compute Intrinsics by Optimization(Method 1)
K_MC1 = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fs, X, '1'),X0,[],[],Options);

% Intrinsics in Matrix Form
K_MC1 = [K_MC1(1) K_MC1(2) K_MC1(3); 0 K_MC1(4) K_MC1(5); 0 0 1];
disp('Intrinsics computed by Mendonca & Cipolla`s Self Calibration(Method 1)')
disp(K_MC1)

% Compute Intrinsics by Optimization(Method 2)
K_MC2 = lsqnonlin(@(X) costFunctionMendoncaCipolla(Fs, X, '2'),X0,[],[],Options);

% Intrinsics in Matrix Form
K_MC2 = [K_MC2(1) K_MC2(2) K_MC2(3); 0 K_MC2(4) K_MC2(5); 0 0 1];
disp('Intrinsics computed by Mendonca & Cipolla`s Self Calibration(Method 2)')
disp(K_MC2)

%% Self Calibration using Classical Kruppa's Method to Find Intrinsics

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolFun', 1e-10,'TolX',1e-10);

% Compute Intrinsics by Optimization
K_CK = lsqnonlin(@(X) costFunctionKruppaClassical(Fs, X),X0,[],[],Options);

% Intrinsics in Matrix Form
K_CK = [K_CK(1) K_CK(2) K_CK(3); 0 K_CK(4) K_CK(5); 0 0 1];
disp('Intrinsics computed by Classical Kruppa`s Self Calibration')
disp(K_CK)

%% Self Calibration using Simplified Kruppa's Method to Find Intrinsics

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt','TolFun', 1e-20,'TolX',1e-20);

% Compute Intrinsics by Optimization
K_SK = lsqnonlin(@(X) costFunctionKruppaSimplified(Fs, X),X0,[],[],Options);

% Intrinsics in Matrix Form
K_SK = [K_SK(1) K_SK(2) K_SK(3); 0 K_SK(4) K_SK(5); 0 0 1];
disp('Intrinsics computed by Simplified Kruppa`s Self Calibration')
disp(K_SK)

%% Self Calibration using Dual Absolute Quadric Method to Find Intrinsics

X0 = [A(1,:) A(2,:) A(3,:)];

% Compute Normal to Plane at Infinity
Normal = NormaltoPlaneAtInfinity(PPM, A);

% Compute Homography of Plane at Infinity
H = homographyPlaneAtInfinity(Fs, Normal);

% Set Optimizer Options for Non-linear Least-square Technique
Options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','off','TolFun', 1e-16,'TolX',1e-16);

% Compute Intrinsics by Optimization
K_DAQ = lsqnonlin(@(X) costFunctionDualAbsoluteQuadric(H, X), X0, [], [], Options);

% Intrinsics in Matrix Form
K_DAQ = [K_DAQ(1) K_DAQ(2) K_DAQ(3); K_DAQ(4) K_DAQ(5) K_DAQ(6); K_DAQ(7) K_DAQ(8) K_DAQ(9)];

% Regain the Scale by setting the Last element to 1
K_DAQ = (1/K_DAQ(3,3))*K_DAQ;
disp('Intrinsics computed by DAQ Self Calibration')
disp(K_DAQ)
