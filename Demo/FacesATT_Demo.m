% ================================================================= %
% This demo shows the performance of LSC on ATT Human Face Images.

% ========================= Acknowledgement =============================
% This code is written by Zhaiming Shen in March 2022 under Dr. Ming-Jun Lai's supervision.
% It is modified based on Dr. Daniel Mckenzie's original code. 
% =======================================================================

clear, clc, close all, format compact

addpath(genpath('../LSC'))
addpath(genpath('../Utilities'))
addpath(genpath('../Datasets'))

load('FacesATT.mat');

% ============== Define all vectors of interest =========== %

k = 10;  %number of clusters
Iter = 10; %number of iterations

Accuracy_LSQ_mat = zeros(k,Iter);
Precision_LSQ_mat = zeros(k,Iter);
Recall_LSQ_mat = zeros(k,Iter);

% =============== Parameters ============== %
epsilon_lsq = 0.6;
reject_lsq = 0.7;

sample_frac = 0.1;  %change here to modify the percentage of seed vertices

%rand=randperm(100);
%L=L(rand,:);
%A=A(rand,rand);
%y=y(rand);

% =========== Find the ground truth clusters ======== %
TrueClusters = cell(k,1);
n0vec = zeros(k,1);
for a = 1:k
    Ctemp = find(y == a);
    TrueClusters{a} = Ctemp;
    n0vec(a) = length(Ctemp);   
end

for j=1:Iter
    for i=1:k
        TrueCluster = TrueClusters{i};
        n0 = length(TrueCluster); 
        n0_equal = 10;
        

        % ================ Draw Seed vertices =============== %
        Gamma = datasample(TrueCluster,ceil(sample_frac*n0_equal),'Replace',false);

        % ================= Run LSC ================= %
        Cluster_LSQ = LeastSquareClustering(A,Gamma,n0_equal,epsilon_lsq,3,reject_lsq);
        
        Precision_LSQ_mat(i,j) = length(intersect(Cluster_LSQ,TrueCluster))/length(Cluster_LSQ);
        Recall_LSQ_mat(i,j) = length(intersect(Cluster_LSQ,TrueCluster))/n0;
        Accuracy_LSQ_mat(i,j) = length(intersect(Cluster_LSQ,TrueCluster))/length(TrueCluster);
        
    end
end

% ==================== Determine Error ============================== %

Precision_LSQ = mean(Precision_LSQ_mat,'all');
Recall_LSQ = mean(Recall_LSQ_mat,'all');
F1_LSQ = 2*Precision_LSQ.*Recall_LSQ./(Precision_LSQ+Recall_LSQ);
Accuracy_LSQ = mean(Accuracy_LSQ_mat,'all');

LSQ = [Precision_LSQ,Recall_LSQ, F1_LSQ, Accuracy_LSQ]
Std = [std(mean(Accuracy_LSQ_mat,1))]