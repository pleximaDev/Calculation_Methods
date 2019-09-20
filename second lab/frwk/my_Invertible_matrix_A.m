% method obratnoy matrici A^(-1)

function [x, ok]=my_Invertible_matrix_A(A, b)
[n, m]=size(A);
ok = false;
%x=zeros(n, 1);
C=zeros(n, m);
if det(A) ~= 0
    ok = true;
    for i = 1 : 1 : n
%         j = 1;
        %M=A;
        for j = 1 : 1 : n               %while j ~= n+1
            M=A;
            M(i, :) = [];
            M(:, j) = [];
            C(i, j) = (-1)^(i+j) * det(M);
            j = j + 1;
        end
    end
    T = (1/det(A)) * C;
    x = T * b;
end
end
