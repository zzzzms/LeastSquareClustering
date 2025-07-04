function score = Jaccard_Score(C1,C2)
%        score = Jaccard_Score(C1,C2)
% This function computes the Jaccard index/similarity between two subsets
% C1 and C2
% =================== Inputs ================= %
% C1 ................... A subset of [n]
% C2 ................... Another Subset of [n]
%
% =================== Output ================ %
% score .............. The JAccard score

top = length(intersect(C1,C2));
bottom = length(union(C1,C2));

score = top/bottom;
end