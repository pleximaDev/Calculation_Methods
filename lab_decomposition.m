clc
clear variables;

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



[L,U] = lu(A)
fprintf("_____________________")
[L,U] = my_lu(A, "Doolittle")
A
L*U

fprintf("_____________________\n")
[L,U] = my_lu(A, "Crout")
A
L*U
fprintf("_____________________\n")




[my, p] = my_chol(A, "Cholesky_Banachiewicz")
[standard] = chol(A, 'lower')




[Q,R] = qr(A);
[Q, R] = my_qr(A);


function [L,U] = my_lu(x, method)
switch method
    case "Doolittle"
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
                sumU = 0;
                L(i,i) = 1;
                if j ~= n
                    j = j + 1;
                end
                L(j,i) = (1/U(i,i)) * (x(j, i) - sumL);
                sumL = 0;
                j = j - 1;
            end    
        end
        
    case "Crout"
        [m, n] = size(x);
        L = zeros(m, n);
        U = zeros(m, n);
        sumL = 0;
        sumU = 0;
        for i = 1 : 1 : n
            for j = i : 1 : n
                for k = 1 : 1 : i-1
                    sumL = sumL + L(j, k) * U(k, i);
                    sumU = sumU + L(i, k) * U(k, j);
                end
                L(j,i) = x(j,i) - sumU;
                sumU = 0;
                U(i,i) = 1;
                if j ~= n
                    j = j + 1;
                end
                U(i,j) = (1/L(i,i)) * (x(i, j) - sumL);
                sumL = 0;
                j = j - 1;
            end    
        end
        
    otherwise 
        fprintf("Error occured while entering method's name.");
end
end

function [L, p] = my_chol(x, method)
switch method
    case "Cholesky_Banachiewicz"
        p = true;
        if all(eig(x) <= 1e-9)
            p = false;
            L = zeros(m, n);
            return
        end
        [m, n] = size(x);
        L = zeros(m, n);
        sumij = 0;
        sumii = 0;
        j = 1;
        k = 1;
        L(1, 1) = (x(1, 1))^(1/2);
        for i = 2 : 1 : m
            j = 1;
            while j <= i
               while k < j
                   sumij = sumij + L(i, k) * L(j, k);
                   k = k + 1;
               end
               k = 1;
               while k < j
                   sumii = sumii + L(i,k) * L(i,k);
                   k = k + 1;
               end
               k = 1;
               
               L(i,j) = (1/L(j,j))*(x(i,j) - sumij);
               L(i,i) = (x(i,i) - sumii)^(1/2);
               sumij = 0;
               sumii = 0;
               j = j + 1;
            end
        end
    otherwise
        fprintf("Error occured while entering method's name.")
end
end

function [Q,R] = my_qr(x)
[m,n] = size(x);
Q = zeros(m,n);
R = zeros(m,n);
sum_proj = zeros(m,1);

for k = 1 : 1 : m
    for j = 1 : 1 : k-1
        sum_proj = sum_proj + ((x(:,k)' * Q(:,j))/(Q(:,j)'*Q(:,j))).*Q(:,j);
    end
    Q(:,k) = x(:, k) - sum_proj;
    Q(:,k) = -Q(:,k)/norm(Q(:,k));
    sum_proj = 0;
end
R = Q' * x;
%___________Unnecessary___________%
zero_logic = triu(ones(m,n));
logic = (zero_logic == 0);
R(logic) = 0;
%___________Unnecessary___________%
end





