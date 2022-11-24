clc
A = [2,-1,1;3,3,9;3,3,5];
n = size(A, 1);
L = eye(n);
for k = 1 : n
    L(k + 1 : n, k) = A(k + 1 : n, k) / A(k, k);
    for l = k + 1 : n
        A(l, :) = A(l, :) - L(l, k) * A(k, :);
    end
end
U = A;
U
L