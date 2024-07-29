% ================================================================= %
% This demo shows the performance of LSC on SBM with different sized clusters.

% ========================= Acknowledgement =============================
% This code is written by Zhaiming Shen in March 2022 under Dr. Ming-Jun Lai's supervision.
% It is modified based on Dr. Daniel Mckenzie's original code. 
% =======================================================================

clear, clc, close all, format compact

addpath(genpath('../LSC'))
addpath(genpath('../Utilities'))
% ============== Parameters ================= %
num_sizes = 5;                      % Number of different cluster sizes
num_trials = 10;                     % Number of trials to run for each size
Cluster_sizes = 200*[1:num_sizes];  % Vector of cluster sizes


% ============ Parameters ========== %
epsilon_lsq = 0.6;
reject_lsq = 0.8;

% ============== Define all matrices of interest =========== %
time_LSQ_mat = zeros(num_trials,num_sizes); Jaccard_LSQ_mat = zeros(num_trials,num_sizes);

for j = 1:num_sizes
    n1 = Cluster_sizes(j);
    n0vec = n1*[1,2,5];
    n = sum(n0vec);
    p_prime = ((log(n))^2)/2;
    q = 5*log(n)/n;
    P_diag = p_prime./n0vec - q;
    P = q*ones(3,3) + diag(P_diag);

%     P = [a*log(n)/n,b*log(n)/n,b*log(n)/n,b*log(n)/n;
%         b*log(n)/n,a*log(n)/n,b*log(n)/n,b*log(n)/n;
%         b*log(n)/n,b*log(n)/n,a*log(n)/n,b*log(n)/n;
%         b*log(n)/n,b*log(n)/n,b*log(n)/n,a*log(n)/n;];
    
%     P = [a*log(n)/n,b*log(n)/n,b*log(n)/n;
%         b*log(n)/n,a*log(n)/n,b*log(n)/n;
%         b*log(n)/n,b*log(n)/n,a*log(n)/n;];
    
    for i = 1:num_trials
        A = generateA2(n0vec,P);
        Im1 = mat2gray(full(A));
        %perm = 1:n;
        perm = randperm(n);
        A = A(perm,perm);
        
        % =============== Find ground truth Cluster ================ %
        [~,permInv] = sort(perm);
        TrueCluster = permInv(1:n1);
        
        % ============== ExtractSeed vertices ================ %
        Gamma1 = datasample(TrueCluster,1,'Replace',false);
        Gamma = datasample(TrueCluster,5,'Replace',false);
        Gam_NBD = find(A(Gamma1,:));
        
        % ========== Find Cluster with LSC =========== %        
        tic
        Cluster_LSQ = LeastSquareClustering(A,Gamma,n1,epsilon_lsq,3,reject_lsq);
        time_LSQ_mat(i,j) = toc;
        Jaccard_LSQ_mat(i,j) = Jaccard_Score(TrueCluster,Cluster_LSQ)
      
    end
end

% ======= Plot all for comparison ======== %
figure,
plot(200*[1:num_sizes],mean(Jaccard_LSQ_mat,1),'LineWidth',3)
legend({'LSC'},'FontSize',14)
ylabel('Jaccard Index')
xlabel('Size of Target Cluster')
set(gca, 'FontSize',14)

% ======= Plot all for times comparison ======== %
figure,
plot(200*[1:num_sizes],log(mean(time_LSQ_mat,1)),'LineWidth',3)
legend({'LSC'},'FontSize',14)
ylabel('logarithm of run time')
xlabel('Size of Target Cluster')
set(gca, 'FontSize',14)
