function [C,v] = LeastSquareClusterPursuit(L,Gamma_a,Omega_a,n_a,reject)
% Subroutine for LsqSingleClusterPursuit that performs the pursuit step.

% ========================= Acknowledgement ==============================
% This code is based on the code of ClusterPursuit algorithm by Dr. Daniel
% Mckenzie, with the subroutine subspacepursuit is replaced by lsqr, and  
% also the parameter 'reject' in Mckenzie's ClusterPursuit is removed. 
% ========================================================================

% INPUT
% ==========================================
% L .................... Laplacian matrix
% Gamma_a .............. Labeled data for C_a
% Omega_a .............. Superset containing C_a
% n_a .................. (Estimated) size of C_a
%
% OUTPUT
% =========================================
% C ................... Estimate of C_a (including Gamma_a)
% v ................... vector of probabilities
% of not being in C
%

Phi = L(:,Omega_a);
%n = size(L,1);
yg = sum(Phi,2);
[~,I]=mink(abs(Phi')*abs(yg),floor(n_a/5)); % change 3 to change the portion of vertices to remove
Phi(:,I)=0;
%g = length(Gamma_a);
%sparsity = ceil(1.1*(length(Omega_a) - (n_a - g)));
%if sparsity <= 0
%    C = union(Omega_a,Gamma_a);
%    v = zeros(n,1);
%else
    %v = lai2012(Phi,yg,1);
    %v = lsqminnorm(Phi,yg);
    v = lsqr(Phi,yg,[],10);
   % v=pinv(Phi)*yg;
   % v = lsqminnorm(Phi,yg)
   % size(v)
   % [~, Iv] = maxk(v,floor(n_a*epsilon));
  %  Lambda_a = Iv;
    Lambda_a = v > reject;
    C = setdiff(Omega_a,Omega_a(Lambda_a));
    C = union(C,Gamma_a);
%end

end
