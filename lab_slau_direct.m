clc
clear variables;
close all force;

A = [1, -0.2589, -0.3093; -0.2589, 1, -0.2705; -0.3093, -0.2705, 1];
b = ones(3, 1);
x = [2.2873; 2.2162; 2.3068]


[x, ok] = lab_slau_gauss(A, b);


[x, ok] = lab_slau_gauss_jordan(A, b)

[x, ok] = lab_slau_Cramer(A, b);



function [x, ok] = lab_slau_gauss(A, b)
ok = false;
[m, n] = size(A);
x = zeros(m, 1);

if n == m
    ok = true;
end

if ok 
    for i = 1 : 1 : n
        for k = i + 1 : 1 : m
            coefficient = (A(k, i)/A(i, i));
            A(k, :) = A(k, :) - A(i, :) * coefficient;
            b(k) = b(k) - b(i) * coefficient;
        end
    end
    for i = n : -1 : 1
        sum = 0;
        for k = i + 1 : 1 : n
            sum = sum + A(i, k) * x(k);
        end
        x(i) = (b(i) - sum)/A(i, i);
    end
end
end

function [x, ok] = lab_slau_gauss_jordan(A, b)
ok = false;
[m, n] = size(A);
x = zeros(m, 1);
C = horzcat(A, b);

if n == m
    ok = true;
end

if ok
    for i = 1 : 1 : n
        for k = i + 1 : 1 : n
            C(k, :) = C(k, :) - C(i, :) * (C(k, i)/C(i, i));
        end
    end
    for i = n : -1 : 1
        for k = i - 1 : -1 : 1
            C(k, :) = C(k, :) - C(i, :) * (C(k, i)/C(i, i));
        end
    end
    for t = 1 : 1 : m
        x(t) = C(t, m + 1)/C(t, t);
    end
    
end
end

function [x, ok] = lab_slau_Cramer(A, b)
[m,n] = size(A);
ok = false;
x = zeros(m, 1);
eps = 1e-16;

if m == n
    if det(A) < eps || det(A) > eps
        ok = true;
    end
end

if ok
    for i = 1 : 1 : n
        M = A;
        M(:, i) = b;
        x(i, 1) = (1/det(A)) * det(M);
    end
    
end
end



