function [A_inv] = inv2(A)

[L,U,P] = lu(A); % LU decomposition with permutation matrix P
A_inv = U \ (L \ P); % inverse of A
end
