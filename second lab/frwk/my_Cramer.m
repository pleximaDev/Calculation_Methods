
% Cramer's rule (Method Cramera)

function [x, ok]=my_Cramer(A, b)
[n, m]=size(A);
ok = true;
epsil = (1e-16);
x=zeros(n, 1);
if n == m % proverka matrici na to, kvadratnaya li ona
    if (-epsil) < det(A) < epsil           %det(A) ~= 0
        ok = false;
    end
    if ok 
        for i = 1 : 1 : n
            M=A;
            M(:, i) = b;
            x(i, 1) = (1/det(A)) * det(M);
        end
    end
end
end