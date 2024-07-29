% ================================================================= %
% This demo shows the performance of LSC on USPS Data.

% ========================= Acknowledgement =============================
% This code is written by Zhaiming Shen in March 2022 under Dr. Ming-Jun Lai's supervision.
% It is modified based on Dr. Daniel Mckenzie's original code. 
% =======================================================================

clear, clc, close all, format compact

addpath(genpath('../LSC'))
addpath(genpath('../Utilities'))
addpath(genpath('../Datasets'))

load('USPS_KNN_Mult_K=15_r=10.mat')

% ============ Parameters ============= %
epsilon_lsq = 0.8;
reject_lsq = 0.85;   


num_trials = 10;
num_sizes = 5;
k = 10;  %number of clusters

% =========== Find the ground truth clusters ======== %
TrueClusters = cell(k,1);
n0vec = zeros(k,1);
for a = 1:k
    Ctemp = find(y== a-1);
    TrueClusters{a} = Ctemp;
    n0vec(a) = length(Ctemp);
    
end


% ============== Define all vectors of interest =========== %
time_LSQ_mat = zeros(k,num_sizes,num_trials);
Cluster_LSQ_mat = cell(k,num_sizes,num_trials);
InterLength_LSQ_mat = zeros(k,num_sizes,num_trials);
Precision_LSQ_mat=zeros(k,num_sizes,num_trials);
Recall_LSQ_mat=zeros(k,num_sizes,num_trials);
F1_LSQ_mat=zeros(k,num_sizes,num_trials);

sample_frac = 0.001;

for m = 1:num_trials
    for j = 1:num_sizes
        for i=1:k
            TrueCluster = TrueClusters{i};
            n0 = length(TrueCluster); 
            %%n0_equal = 7000;

            % ================ Draw Seed set =============== %
            Gamma = datasample(TrueClusters{i},ceil(sample_frac*j*n0),'Replace',false);
        
            tic
            Cluster_LSQ_mat{i,j,m} = LeastSquareClustering(A,Gamma,n0,epsilon_lsq,3,reject_lsq);
            time_LSQ_mat(i,j,m) = time_LSQ_mat(i,j,m) + toc;
       
            Precision_LSQ_mat(i,j,m) = length(intersect(Cluster_LSQ_mat{i,j,m},TrueClusters{i}))/length(Cluster_LSQ_mat{i,j,m})
            Recall_LSQ_mat(i,j,m) = length(intersect(Cluster_LSQ_mat{i,j,m},TrueClusters{i}))/length(TrueClusters{i})
            F1_LSQ_mat(i,j,m) = 2*Precision_LSQ_mat(i,j,m)*Recall_LSQ_mat(i,j,m)/(Precision_LSQ_mat(i,j,m)+Recall_LSQ_mat(i,j,m));
          
            InterLength_LSQ_mat(i,j,m) = length(intersect(Cluster_LSQ_mat{i,j,m},TrueClusters{i}));

            length(Cluster_LSQ_mat{i,j,m})
        
        end
    end
end

% =================== Determine Error ===================== %
Precision_LSQ = mean(Precision_LSQ_mat,[1,3])
Recall_LSQ = mean(Recall_LSQ_mat,[1,3])
F1_LSQ = mean(F1_LSQ_mat,[1,3]) 

Accuracy_LSQ = sum(InterLength_LSQ_mat,[1,3])/(num_trials*11000)
Std_LSQ = std(mean(Recall_LSQ_mat,1),0,3)

Time_LSQ = mean(time_LSQ_mat,[1,3])

