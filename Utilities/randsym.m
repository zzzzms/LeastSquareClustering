function A = randsym(n,p)
% This function generates a random symmetric binary matrix with expected
% row sums np and zeros on the diagonal. 
% Use instead of sprandsym when p is close to 1.

AA = rand(n);
T = +(triu(AA,1)> 1-p);
A = T + T';
end
