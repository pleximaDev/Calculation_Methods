clc;
clear variables;
close all force;

%2.1
A= randn(4,4);
A=A*A'
b=randn(4,1);
v=[5 4 2 ; 3 1 3; 4 2 1] %test
v=v*v'
disp(A);
disp(b);
% AA
c=A*A;
% AA'
c=A*A';
% A/A
c=A/A;
% A^4
c=A^4;
% bb' 
c=b * b';
% b'b
c=b' * b;
% b.b
c=b' * b;
%bxb

% element-wise multiplication
c=A .* A;
% element-wise division A/A
c= A ./ A;
% bitwise matrix cube A
c= A.^3;
% element-wise multiplication b b
c=b .* b;
% A^-1(Using I, and using pow)
c=A^(-1);
c=eye(4)/A;
%A b
c=A*b;
% delete A's 2nd line and 2nd column, and b's 3rd element
A(2, :)=[];
A(:, 2)=[]
b(3)=[]
B=A;
% Matrix with ones if corresponding element is positive and -5 otherwise (using logical indexing)
A(A>0)=1;
A(A<0)=-5;
C=A
A=B;
% Making all C<0 elements equal zero
C(C<0)=0

%2.2
% b · b (using dot())
c=dot(b, b);
% b x b (using cross())
c=cross(b, b);
% A's rank (using rank())
c=rank(A);
% min value of A (using min())
c=min(A);                                      
c=min(c); %only one value
% max A's value (using max())
c=max(A);
c=max(c);
% eigenvectors and eigevalues of A (using eig())
[V,D] = eig(A); % D's diagonal consists of eigenvalues; V's columns consists of eigenvectors 

% Singular values and vectors of A (using svds())
[U,S,V] = svds(A); % U - left-singular vectors, V - right-singular vectors, S - Unitary matrix with right-singular vectors in columns

% A's matrix exponent (using expm(), and using cycle);
c=expm(A);
c=0;
for k=0:1:10
    c=c+(1/factorial(k))*A^k;
end
c;
%операторные нормы ||A||1 и ||Ak||inf (как с использованием norm(),так и с помощью цикла);
c=norm(A , 1);
max=0;
for j = 1 : 1 : 3
    sm=0;
    for i =1 : 1 : 3
    c=A(i, j);
    sm=sm+c;
    end
    if max<sm
        max=sm;
    end
end
max;
c=norm(A, inf);
max=0;
for i = 1 : 1 : 3
    sm=0;
    for j =1 : 1 : 3
    c=A(i, j);
    sm=sm+c;
    end
    if max<sm
        max=sm;
    end
end
max;
% норму Фробениуса kAkF (как с использованием norm(), так и с помощью
%цикла, так и без использования цикла - через векторно-матричные операции
%и функцию sum());
c = norm(A,'fro');

c=A.^2;
c=sum(sum(c));
c=c^(1/2);

sm=0;
for i = 1 : 1 : 3
    for j =1 : 1 : 3
        sm=sm+A(i,j)^2;
    end
end
sm=sm^(1/2);
        


% след матрицы A (как с использованием trace(), так и с помощью цикла);
c=trace(A);
c=sum(diag(A));
c=0;
for i = 1 : 1 : 3
   c=c+A(i, i);
end
c;
% диагональную матрицу D, состоящую из собственным значений матрицы A (используя diag());
D=eig(A);
D=diag(D);
% верхнюю и нижнюю треугольные части матрицы A (используя tril() и triu());
c=tril(A);
c=triu(A);
% определить количество ненулевых элементов в матрице C (используя nnz()).
c=nnz(C);


%3.1
% спектральное разложение (используя eig())
A;
[V,D] = eig(A);
A = V * D * V^(-1);

% сигнулярное разложение (используя svd())
[U,S,V] = svds(A); %U - левые синг вект-ры, V правые синг вект, S синг числа
A;
A=U * S * conj(V)';

% LU-разложение (используя lu())

%L=tril(A);
%U=triu(A);
[L,U]= lu(A);
A;
N = L * U;

% LL-разложение (используя chol())
A;
L = chol(A,'lower');
R = chol(A,'upper');
A= L * R;
 
% LDL-разложение (используя ldl())
A;
[L,D] = ldl(A);
A = L * D * conj(L)';

% разложение Шура (используя schur())
A;
[U,T] = schur(A);
A = U * T * conj(U)';
    

% QR-разложение (используя qr())
A;
[Q,R] = qr(A);
A =  Q * R;


% 3.2 

% LU-decomposition

[L, U]=my_lu(A);
[L, U]=lu(A);

% LL-decomposition
L=my_chol(A);
A;
c=L*conj(L)';
L=chol(A, 'lower');


% QR-decomposition
A;
[Q,R]=my_qr(A);
[Q,R] = qr(A);


function[L,U]= my_lu(x)
U=zeros(3, 3);
L=zeros(3, 3);
for j = 1 : 1 : 3
    U( 1, j) = x(1, j);
    L(j, 1) = x(j, 1)/U(1, 1);
end
SMu=0;
SMl=0;
for i = 2 : 1 : 3
    for j = i : 1 : 3
        for k = 1 : 1 : i-1
            SMu=SMu+L(i, k) * U(k, j);
            SMl=SMl+L(j, k) * U(k, i);
        end
            U(i, j) = x(i, j) - SMu;
            L(j, i) = (1/U(i, i)) * (x(j, i) - SMl);
            SMu=0;
            SMl=0;
    end
    
end
end

function L=my_chol(x)
L=zeros(3, 3);
SMl2=0;
SMll=0;
j=1;
k=1;
t=1;
L(1, 1) = (x(1, 1))^(1/2);
for i = 2 : 1 : 3
            while j < i
                while k < j
                SMll=SMll+L(i, k) * L(j, k);
                k = k + 1;
                end
                L(i, j) = (1/L(j, j)) * (x(i, j) - SMll);
                j = j + 1;
                SMll=0;
            end
            while t < i

            SMl2=SMl2+(L(i, k))^2;
            t = t + 1;
            end
            L(i, i) = (x(i, i) - SMl2)^(1/2);
            SMl2=0;
            if j == i 
                j = 1;
            end
            t = 1;
end
end


function [Q, R]=my_qr(x)
b=zeros(3, 1);
Q=zeros(3, 3);
B=zeros(3, 3);
SMproj=zeros(3, 3);
    b = x(:, 1);
    g=b;
    g = g / ((g' * g)^(1/2));
    Q(:, 1) = -g;
    t=1;

for k = 2 : 1 : 3
    a = x(:, k);
        while t < k
            proj = (-x(:, k)' * Q(:, t)/(Q(:, t)' * Q(:, t)))*Q(:, t);
            SMproj(:, t) = proj;
            t = t + 1;
        end
        t=1;  
        l=1;
        [n, m]=size(SMproj);
        b=-a
        while l < n
            b = b - SMproj(:, l);
            l = l + 1;
        end
    g=b;
    g = g / (g' * g)^(1/2);
    B(:, k)=b;
    Q(:, k) = -g;
    SMproj=zeros(3, 3);
end
R = Q' * x;

end



% use this
%{
function Q=my_qr(x)
b=zeros(3, 1);
Q=zeros(3, 3);
B=zeros(3, 3);
SMproj=zeros(3, 3);
    b = x(:, 1);
    g=b;
    g = g / ((g' * g)^(1/2));
    Q(:, 1) = -g;
    %j=1;
    proj=0;
    t=1;

for k = 2 : 1 : 3
    a = x(:, k);
        while t < k
            testt=x(:, k)' * Q(:, t)/(Q(:, t)' * Q(:, t))
            Q(:, t)
            
            proj = (x(:, k)' * Q(:, t)/(Q(:, t)' * Q(:, t)))*Q(:, t)
            
            SMproj(:, t) = -proj
            heh=SMproj(:, t)
            t = t + 1;
        end
        t=1;  
    
        l=1;
        [n, m]=size(SMproj);
        while l < n
            SMproj
            ts=SMproj(:, l)
            b=a
            
            b = b - SMproj(:, l)
            l = l + 1;
        end
        
    g=b;
    g = g / (g' * g);
    B(:, k)=b;
    Q(:, k) = g;
    SMproj=zeros(3, 3);
end
%R = Q' * A;

end
%}

%{
function Q=my_qr(x)
b=zeros(3, 1);
Q=zeros(3, 3);
B=zeros(3, 3);
SMproj=zeros(3, 3);
    b = x(:, 1);
    g=b;
    g = g / ((g' * g)^(1/2));
    Q(:, 1) = g;
    %j=1;
    proj=0;
    t=1;

for k = 2 : 1 : 3
    a = x(:, k);
    %while j < k
        while t < k
            testt=x(:, k)' * Q(:, t)/(Q(:, t)' * Q(:, t))
            Q(:, t)
            
            proj = (x(:, k)' * Q(:, t)/(Q(:, t)' * Q(:, t)))*Q(:, t)
            
            SMproj(:, t) = proj
            heh=SMproj(:, t)
            t = t + 1;
        end
        t=1;  
        %j = j + 1;
        l=1;
        [n, m]=size(SMproj);
        while l < n
            SMproj
            ts=SMproj(:, l)
            
            b = a - SMproj(:, l)
            l = l + 1;
        end
        
    g=b;
    g = g / (g' * g);
    B(:, k)=b;
    Q(:, k) = g;
    SMproj=zeros(3, 3);
    %j=0;
    %end
end
%R = Q' * A;

end
%}



%rabochaya napolovinu
%{
function L=my_chol(x)
L=zeros(3, 3);
SMl2=0;
SMll=0;
j=1;
L(1, 1) = (x(1, 1))^(1/2);
for i = 1 : 1 : 3
    for k = 1 : 1 : i-1
        %for t= 1 : 1 : j-1
        
        SMll=SMll+L(i, k) * L(j, k);
        
            while j < i
                L(i, j) = (1/L(j, j)) * (x(i, j) - SMll);

                j = j + 1;
            end
            l=1;
            while l ~= i - 1
                SMl2=SMl2+L(i, l)^2;
                l = l + 1;
            end
            l=1;
            L(i, i) = (x(i, i) - SMl2)^(1/2);
            
            SMl2=0;
            if j == i 
                
                j = 1;
            end
               
            
        
        %end
    end
    
end
end
%}


%{
function L=my_chol(x)
L=zeros(3, 3);
SMl2=0;
SMll=0;
j=1;
L(1, 1) = (x(1, 1))^(1/2);
for i = 1 : 1 : 3
    for k = 1 : 1 : i-1
        %for t= 1 : 1 : j-1
        
        SMll=SMll+L(i, k) * L(j, k);
        SMl2=SMl2+L(i, k)^2;
            while j < i
                L(i, j) = (1/L(j, j)) * (x(i, j) - SMll);

                j = j + 1;
            end
            L(i, i) = (x(i, i) - SMl2)^(1/2);
            if j == i 
                
                j = 1;
            end
               
            
        SMl2=0;
        %end
    end
    
end
end
%}
%{

function Q=my_qr(x)
b=zeros(3, 1);
Q=zeros(3, 3);
    b = x(:, 1);
    g=b;
    g = g / (g' * g);
    Q(:, 1) = g;

for j = 2 : 1 : 3
    b = x(:, j);

    g = b - (((g'* b)/(b'* b))' * b)' * b;
    g = g / (g' * g);
    Q(:, j) = g;
end
%R = Q' * A;

end

%}




