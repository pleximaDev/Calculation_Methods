clc
clear variables;




A = [1, -0.2589, -0.3093; -0.2589, 1, -0.2705; -0.3093, -0.2705, 1];
b = ones(3, 1);
x = [2.2873; 2.2162; 2.3068];

% 
% A = rand(4)
% b = rand(size(A, 1), 1)
x0 = zeros(size(A, 1), 1);
eps = 1e-10;
Kmax = 1000;




[x, ok, k] = lab_slau_jacobi(A, b, x0, eps, Kmax)


[x, ok, k] = lab_slau_gauss_seidel(A, b, x0, eps, Kmax)

[x, ok, k] = lab_slau_sor(A, b, x0, eps, Kmax)

[x, ok, k] = lab_slau_mcg(A, b, x0, eps, Kmax)

[x, ok, k] = lab_slau_mbcg(A, b, x0, eps, Kmax)

    scalar = @(q, d) dot(conj(q)', d);
    q = rand(4, 1)
    d = rand(4, )
    scalar(3, 6)
dot(3, 6)

function [x, ok, k] = lab_slau_jacobi(A, b, x0, eps, Kmax)
[m, n] = size(A);
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
k = 0;
x_k = x0;
x_kplus1 = b;

x = x0;
B = -D^(-1) * (L + U);
g = D^(-1) * b;

% full() for sparse matrix
% if (eig(full(B)) < 1); ok = true; else; ok = false; end
ok = (all(eig(full(B)) < 1));


if ok 
    while(k < Kmax) && (all(abs(x_kplus1 - x_k)) > eps)
        x_k = x_kplus1;
        for i = 1 : 1 : m
            sum = 0;
            for j = 1 : 1 : n
                if j ~= i
                    sum = sum + A(i, j) * x_k(j);
                end
            end
            x_kplus1(i) = 1/A(i, i) * (b(i) - sum);
        end
        k = k + 1;
    end
end

x = x_k;

end

function [x, ok, k] = lab_slau_gauss_seidel(A, b, x0, eps, Kmax)
[m, n] = size(A);
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
x = x0;
k = 0;

P = -(L + D)^(-1) * U;
ok = norm(full(P)) < 1;

if ok
    while (k < Kmax) && (norm(A * x - b) > eps)
        for i = 1 : 1 : m
            d(i) = b(i)/A(i, i);
            sum1 = 0;
            for j = 1 : 1 : i - 1
                c(i, j) = (-A(i, j)/A(i, i)) * (j ~= i) + 0 * (j == i);
                sum1 = sum1 + c(i, j) * x(j);
            end
            sum2 = 0;
            for j = i + 1 : 1 : n
                c(i, j) = (-A(i, j)/A(i, i)) * (j ~= i) + 0 * (j == i);
                sum2 = sum2 + c(i, j) * x(j);
            end
            x(i) = sum1 + sum2 + d(i);
        end
        k = k + 1;
    end
end
end


function [x, ok, k] = lab_slau_sor(A, b, x0, eps, Kmax)
%Successive over-relaxation
[m, n] = size(A);
D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);
x = x0;
k = 0;
T = eye(m, n) - D^(-1) * A;
ro = max(abs(eigs(T)));
omega = 1 - (ro/(1 + sqrt(1 - ro^2)))^2;

P = -(L + D)^(-1) * U;
ok = norm(full(P)) < 1;

if ok
    while (k < Kmax) && (norm(A * x - b) > eps)
        for i = 1 : 1 : m
            sum1 = 0;
            for j = 1 : 1 : i - 1
                sum1 = sum1 + A(i, j) * x(j);
            end
            sum2 = 0;
            for j = i + 1 : 1 : n
                sum2 = sum2 + A(i, j) * x(j);
            end
            x(i) = (1 - omega) * x(i) + omega/A(i, i) * (b(i) - sum1 - sum2);
        end
        k = k + 1;
    end
end
end

function [x, ok, k] = lab_slau_mcg(A, b, x0, eps, Kmax)
    x = x0;
    r = b - A * x;
    z = r;
    k = 0;
    
    ok = isreal(A) * (all(eig(A) > 1));
    
    while(k < Kmax) &&  (norm(A * x - b) > eps)
        alpha = dot(r, r)/dot(A*z, z);
        x = x + alpha * z;
        temp = dot(r, r);
        r = r - alpha * A * z;
        betta = dot(r, r)/(temp);
        z = r + betta * z;
        k = k + 1;
    end
end

function [x, ok, k] = lab_slau_mbcg(A, b, x0, eps, Kmax)

    x = x0;
    r = b - A * x0;
    p = r;
    z = r;
    s = r;
    
    k = 0;
    ok = isreal(A) * (all(eig(A) > 1));
    
    while(k < Kmax) &&  (norm(A * x - b) > eps)
        alpha = dot(p, r)/dot(s, A * z);
        x = x + alpha * z;
        temp = dot(p, r);
        r = r - alpha * A * z;
        p = p - alpha * A' * s;
        betta = dot(p, r)/temp;
        z = r + betta * z;
        s = p + betta * s;
        
        k = k + 1;
    end
end
    

function [x, ok, k] = lab_slau_smbcg(A, b, x0, eps, Kmax)
    
    r = b - A * x;
    r_tilda = r;
    ro = 1;
    alpha = 1; 
    omega = 1;
    v = 0;
    p = 0;
    
    scalar = @(q, d) dot(conj(q)', d);

    while(k < Kmax) &&  (norm(A * x - b) > eps)
        ro_prev = ro;
        ro = dot(r_tilda, r);
        betta = (ro/ro_prev)*(alpha/omega);
        p = r + betta * (p - omega * v);
        v = A * p;
        alpha = ro/dot(r_tilda, v);
        s = r - alpa * v;
        t = A * s;
        
    end
    
    
    
end
