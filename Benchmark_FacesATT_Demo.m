clear, close all, clc

%This code is to use the least squares clustering method to sort a set of human faces from ATT. 
%It is written by Zhaiming Shen under Dr. Ming-Jun Lai's supervision in 2021. 
%Please refer to our paper "A Least Square Approach to Semi-supervised Local Cluster Extraction" for more details.

% =================== Parameters and load the data ===================== %
load('adj_matrix.mat');

rand=randperm(100);
L=L(rand,:);
A=A(rand,rand);
y=y(rand);

% ===== split images and save files ===== %
image=cell(100,1);
for i=1:100
    Ltemp=L(i,:);
    image{i}=reshape(Ltemp,[],46); 
    
    %save each individual image as separate file
    imwrite(image{i}, ['p',num2str(i),'.png']);
end

% ===== show the combined images before clustering ===== %
M=zeros(560,700);
J=0;
for i=1:10
    for j=1:10
       J=J+1; 
    fname=['p',num2str(J),'.png'];
    %A=imread(fname);
    M(56*(i-1)+1:56*i,46*(j-1)+1:46*j)=imread(fname);
    end
end

% show the combined image before clustering
figure, imshow(uint8(M))


% ===== for LsqSingleClusterPursuit ===== %
epsilon_lsq = 0.6;
reject_lsq = 0.7;

M=zeros(560,700);  %initialize canvas
k = 10;  %number of clusters
sample_frac = 0.2;  %change here to modify the percentage of seed vertices

% =========== Find the ground truth clusters ======== %
TrueClusters = cell(k,1);
n0vec = zeros(k,1);
for a = 1:k
    Ctemp = find(y == a);
    TrueClusters{a} = Ctemp;
    n0vec(a) = length(Ctemp);   
end


% ============== Define all vectors of interest =========== %
Precision_LSQ_mat = zeros(k,1);
Recall_LSQ_mat = zeros(k,1);


for i=1:k
    TrueCluster = TrueClusters{i};
    n0 = length(TrueCluster); 
    n0_equal = 10;

    % ================ Draw Seed vertices =============== %
    Gamma = datasample(TrueCluster,ceil(sample_frac*n0_equal),'Replace',false);

    % ================= LsqSingleClusterPursuit ================= %
    Cluster_LSQ = LeastSquareClustering(A,Gamma,n0_equal,epsilon_lsq,3,reject_lsq);
    
    Precision_LSQ_mat(i,1) = length(intersect(Cluster_LSQ,TrueCluster))/length(Cluster_LSQ);
    Recall_LSQ_mat(i,1) = length(intersect(Cluster_LSQ,TrueCluster))/n0;
  
    % ================ Draw the image after Clustering =============== %
    for j=1:length(Cluster_LSQ)
        fname=['p',num2str(Cluster_LSQ(j)),'.png'];
        M(56*(i-1)+1:56*i,46*(j-1)+1:46*j)=imread(fname);
     
    end
end

%show the clustered image in each row
figure, imshow((uint8(M)))

%Precison, Recall, F1 Score, and Accuracy
Precision_LSQ = mean(Precision_LSQ_mat)
Recall_LSQ = mean(Recall_LSQ_mat)
F1_LSQ = 2*Precision_LSQ.*Recall_LSQ./(Precision_LSQ+Recall_LSQ)



