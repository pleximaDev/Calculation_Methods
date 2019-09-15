clc
clear all variables;

% Init
A = randn(4, 4);
b = randn(4, 1);

A = A * A'; % Making A positive-defined


%{
format long % short
format short
pi
c = pi - 5e-5
%}



%____________________________________________%
% Eigendecomposition of a matrix (spectral)
[V, D] = eig(A); %V - eigenvectors, D - eigenvalues
% [V2, D2, flag] = eigs(A) % differencies in positions of columns
Q = V;
L = D;
Q_inv = Q^(-1);

Ans = Q*Q_inv;
(Ans < 0.8);
Ans(Ans < 0.8) = 0;

disp(A)
B = Q * L * Q^(-1)
%____________________________________________%


%____________________________________________%
% Singular value decomposition
[U,S,V] = svd(A) 
% U - Unitary matrix with left-singular vectors;
% S - rectangular diagonal matrix with singular values;
% V - Unitary matrix with right-singular vectors in columns;
disp(A)
B = U * S * conj(V)'
%____________________________________________%


%____________________________________________%
% LU decomposition
[L,U] = lu(A)
disp(A)
B = L * U
%____________________________________________%


%____________________________________________%
% LUP decomposition
[L,U,P] = lu(A)
disp(A)
B = P' * L * U
%____________________________________________%


%____________________________________________%
% LL decomposition
L = chol(A,'lower')
U = chol(A, 'upper')
disp(A)
B = L * U
B = L * conj(L)'
B = conj(U)' * U
%____________________________________________%


%____________________________________________%
% LDL decomposition
[L,D] = ldl(A)
disp(A)
B = L * D * conj(L)'
%____________________________________________%


%____________________________________________%
% QR decomposition
[Q,R] = qr(A)
disp(A)
B = Q * R
%____________________________________________%

[L,U] = my_lu(A)


function [L,U] = my_lu(x)
[m, n] = size(x);
L = zeros(m, n);
U = zeros(m, n);
sumL = 0;
sumU = 0;
for i = 1 : 1 : n
    for j = i : 1 : n
        for k = 1 : 1 : i-1
            sumU = sumU + L(i, k) * U(k, j);
            sumL = sumL + L(j, k) * U(k, i);
        end
        U(i,j) = x(i,j) - sumU;
        L(i,i) = 1;
        if j>n
            break
        end
        L(j+1,i) = (1/U(i,i)) * (x(j+1, i) - sumL)
    end    
end
end




