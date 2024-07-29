% This demo shows the performance of LSC on a small SSBM graph.

% ========================= Acknowledgement =============================
% This code is written by Zhaiming Shen in March 2022 under Dr. Ming-Jun Lai's supervision.
% It is modified based on Dr. Daniel Mckenzie's original code. 
% =======================================================================

clear, close all, clc, warning off, format compact
addpath(genpath('../LSC'))
addpath(genpath('../Utilities'))

% ============== Generating adjacency matrix of graph ============== %
n = 400;          % Set the number of vertices of the graph
k = 4;            % number of clusters
n0 = ceil(n/k);   % size of each cluster (equally sized)
p = 0.1;          % in-cluster connection probability
q = 0.01;         % between cluster connection probability
A = generateA(n,n0,p,q);

% uncomment below if you wish to visualize the adjacency matrix
%Im1 = mat2gray(full(A));
%imshow(imcomplement(Im1));
%title('The ground truth adjacency matrix')



% ================= Parameters ===================
epsilon_lsq = 0.6;   
reject_lsq = 0.8;    

Iter = 1; % number of iteration

Cluster_LSQ = cell(Iter,1); Accuracy_LSQ = zeros(1,Iter);

Time_LSQ = 0;


% ================ Randomly permute the adjacency matrix =============== %
perm = sort(randperm(n));
%perm = randperm(n);
A = A(perm,perm);
[~,permInv] = sort(perm);
TrueCluster = permInv(1:n0);   % the ground truth first cluster, after permutation.

for i = 1:Iter
    % ========================== Draw seed vertices =================== %
    Gamma = datasample(TrueCluster,1,'Replace',false); % seed vertices
    
    % ========= Visualize the graph with seed vertices highlighted ======== %
    G = graph(A);
    figure
    H = plot(G,'Layout','force','MarkerSize',4);
    highlight(H,Gamma,'NodeColor','r','MarkerSize',8);
    title('Graph with seed vertices highlighted','FontSize',14)

    % ====================== Run LSC ================================= %
    tic
    Cluster_LSQ{i} = LeastSquareClustering(A,Gamma,n0,epsilon_lsq,3,reject_lsq);
    Time_LSQ = Time_LSQ + toc;
    
    length(Cluster_LSQ{i})
    
    % ================= Assess Accuracy ===================== %
   
    Accuracy_LSQ(i) = 100*length(intersect(Cluster_LSQ{i},TrueCluster))/length(TrueCluster);
   
end

Time_LSQ = Time_LSQ/Iter;

disp(['Found Cluster 1 by LSC with an accuracy of ',num2str(mean(Accuracy_LSQ)),'%', ' with average time ', num2str(Time_LSQ), ' seconds'])

% ===============  Replot graph =========================== %
figure
H2 = plot(G,'Layout','force','MarkerSize',4);
highlight(H2,Cluster_LSQ{Iter},'NodeColor','r','MarkerSize',4);
highlight(H2,Gamma,'NodeColor','r','MarkerSize',8);
title('Cluster Found by LSC.')
