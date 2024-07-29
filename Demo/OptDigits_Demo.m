% ================================================================= %
% This demo shows the performance of LSC on OptDigits Data.

% ========================= Acknowledgement =============================
% This code is written by Zhaiming Shen in March 2022 under Dr. Ming-Jun Lai's supervision.
% It is modified based on Dr. Daniel Mckenzie's original code. 
% =======================================================================


clear, clc, close all, format compact

addpath(genpath('../LSC'))
addpath(genpath('../Utilities'))
addpath(genpath('../Datasets'))

load('OptDigits_KNN_Mult_K=15.mat')

% =============== Parameters ============== %
epsilon_lsq = 0.6;
reject_lsq = 0.8;

num_trials = 10;
num_sizes = 5;
k = 10;  %number of clusters
n = size(A,1);

% =========== Find the ground truth clusters ======== %
TrueClusters = cell(k,1);
n0vec = zeros(k,1);
for a = 1:k
    Ctemp = find(y== a-1);
    TrueClusters{a} = Ctemp;
    n0vec(a) = length(Ctemp);
    
end


% ============== Define all vectors of interest =========== %
time_LSC_mat = zeros(k,num_sizes); Jaccard_LSC_mat = zeros(k,num_sizes);


for m = 1:num_trials
for j = 1:num_sizes
    sample_frac = 0.005*j;
    for i=1:k
        TrueCluster = TrueClusters{i};
        n0 = length(TrueCluster);
        n0_test = 560;

        % ================ Draw Seed set =============== %
        Gamma = datasample(TrueCluster,ceil(sample_frac*n0),'Replace',false);
        
        tic
        Cluster_LSC = LeastSquareClustering(A,Gamma,n0,epsilon_lsq,3,reject_lsq);
        time_LSC_mat(i,j) = time_LSC_mat(i,j) + toc;
        Jaccard_LSC_mat(i,j) = Jaccard_LSC_mat(i,j) + Jaccard_Score(TrueCluster,Cluster_LSC)

    end
end
end

time_LSC_mat = time_LSC_mat/num_trials;
Jaccard_LSC_mat = Jaccard_LSC_mat/num_trials;

figure,
plot(0.5*[1:5],mean(Jaccard_LSC_mat),'LineWidth',3)
xlabel('Percentage of vertices used as seeds','FontSize',14)
ylabel('Average Jaccard Score','FontSize',14)
legend({'LSC'},'FontSize',14)

figure,
plot(0.5*[1:5],log(mean(time_LSC_mat)),'LineWidth',3)
xlabel('Percentage of vertices used as seeds','FontSize',14)
ylabel('Logarithm of Average Time in seconds','FontSize',14)
legend({'LSC'},'FontSize',14)
