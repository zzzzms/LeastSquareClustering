% ================================================================= %
% This demo shows the performance of LSC on political blog network.

% ========================= Acknowledgement =============================
% This code is written by Zhaiming Shen in March 2022 under Dr. Ming-Jun Lai's supervision.
% It is modified based on Dr. Daniel Mckenzie's original code. 
% =======================================================================

clear, clc, close all, format compact

addpath(genpath('../LSC'))
addpath(genpath('../Utilities'))
addpath(genpath('../Datasets'))

load('polblogs.mat');

% ================ Randomly permute the adjacency matrix =============== %
A=polblogsAdjMat;
[~,n]=size(A);
perm = randperm(n);
A = A(perm,perm);
[~,permInv] = sort(perm);
TrueCluster = permInv(1:586);   % the ground truth first cluster, after permutation.

% ========================== Draw seed vertices =================== %
Gamma = datasample(TrueCluster,1,'Replace',false); % seed vertices

% =============== Parameters ============== %
epsilon = 0.8;    
reject = 0.1;    

% ====================== Run LSC ================================= %
tic
Cluster = LeastSquareClustering(A,Gamma,586,epsilon,3,reject);
toc


% ================= Assess Accuracy ===================== %
accuracy = 100*length(intersect(Cluster,TrueCluster))/length(TrueCluster);
accuracy = 100*(1-(length(Cluster)-586*accuracy/100)*2/1224);

disp(['Found Cluster 1 by LSC with an accuracy of ',num2str(accuracy),'%'])

misclassified = length(Cluster)-length(intersect(Cluster,TrueCluster))

