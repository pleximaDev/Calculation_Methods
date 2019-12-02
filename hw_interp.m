clc;
clear variables;
close all force;

a = 0;
b = 1;

Na = 8;
Nb = 10 * Na;

x1 = a : (b - a)/Na : b; %interpolated grid
x0 = a : (b - a)/Nb : b; %interpolation grid

N1 = length(x1);
N0 = length(x0);

v = 1; % frequency
omega = 2 * pi * v; 
f0 = interp_function(omega, x0); % analytic function model
f1 = interp_function(omega, x1); % experimental data model

[f2, x0] = nearest_neighbor(f1, x1, x0);
[f3, x0] = linear_interpolation(f1, x1, x0);
[f4, x0] = lagrange_polynomial(f1, x1, x0);
[f5, x0, forward_diff, backward_diff] = Newton_polynomial(f1, x1, x0, 'forward');
[f6, x0, forward_diff, backward_diff] = Newton_polynomial(f1, x1, x0, 'backward');
[f7, x0] = cubic_spline(f1, x1, x0);


figure('Name','Plots','NumberTitle','off');
clf;
subplot(2, 1, 1);

hold on;
grid on;
grid minor;

stem(x0, f0, 'LineWidth', 1.5); 
stem(x1, f1, 'LineWidth', 1.5);
stairs(x0, f2, 'LineWidth', 1.5);
plot(x0, f3, 'LineWidth', 1.5);
plot(x0, f4, 'LineWidth', 1.5);

title('Interpolation');
std = {'Experimental data model','Analytical model function'};
names = horzcat(std,'Nearest-neighbor interpolation');
names = horzcat(names,'Linear interpolation','Lagrange polynomial');
legend(names);
legend('location','northeastoutside');
ylabel('f(x)');
xlabel('x');
hold off;

subplot(2, 1, 2);
hold on
plot(x0, f5, 'LineWidth', 1.5);
plot(x0, f6, 'LineWidth', 1.5);
plot(x0, f7, 'LineWidth', 1.5);
stem(x0, f0, 'LineWidth', 1.5);
stem(x1, f1, 'LineWidth', 1.5);
ylim([-1, 1]);
grid on;
grid minor;
names = {'Newton polynomial forward','Newton polynomial backward'};
names = horzcat(names, 'Cubic spline interpolation');
names = horzcat(names, std);
legend(names);
legend('location','northeastoutside');
ylabel('f(x)');
xlabel('x');
hold off;


function f = interp_function(omega, x)
f = sin(omega * x);
end

function [f, x0] = nearest_neighbor(f1, x1, x0)
k = 1;
f = zeros(length(x0), 1);
for i = 1 : 1 : length(x0)
    if k < length(x1)
        delta1 = abs(x0(i) - x1(k));
        delta2 = abs(x0(i) - x1(k + 1));
        if delta2 <= delta1
            k = k + 1;
        end
    end
    f(i, 1) = f1(k);
end
end

function [f, x0] = linear_interpolation(f1, x1, x0)
k = 1;
f = zeros(length(x0), 1);
for i = 1 : 1 : length(x0)
    if x0(i) >= x1(k + 1) && k ~= length(x1) - 1
        k = k + 1;
    end
    f(i, 1) = f1(k) + ((f1(k + 1) - f1(k))/(x1(k + 1) - x1(k))) * (x0(i) - x1(k));
end

end

function [f, x0] = lagrange_polynomial(f1, x1, x0)
f = zeros(length(x0), 1);
shift = 1;
for l = 1 : 1 : length(x0)
    for k = 0 : 1 : length(x1) - shift
        basis_polynomial = 1;
        for i = 0 : 1 : length(x1) - shift
            if(i ~= k)
                basis_polynomial = basis_polynomial * (x0(l) - x1(i + shift))/(x1(k + shift) - x1(i + shift));
            end
        end
        f(l, 1) = f(l, 1) + basis_polynomial * f1(k + shift);
        
    end
end
end


function [f, x0, forward_diff, backward_diff] = Newton_polynomial(f1, x1, x0, Newton)
% finite_diff matrix here by one loop

f=zeros(length(x0), 1);
n = length(x1) - 1;

%============================================================%
% 1  2 - 1  2 - 1  2 - 1  2 - 1  2 - 1  2 - 1  2 - 1  2 - 1 ||
% 2  3 - 2  3 - 2  3 - 2  3 - 2  3 - 2  3 - 2  3 - 2        ||
% 3  4 - 3  4 - 3  4 - 3  4 - 3  4 - 3  4 - 3               ||
% 4  5 - 4  5 - 4  5 - 4  5 - 4  5 - 4                      ||
% 5  6 - 5  6 - 5  6 - 5  6 - 5                             ||
% 6  7 - 6  7 - 6  7 - 6                                    ||
% 7  8 - 7  8 - 7                                           ||
% 8  9 - 8                                                  ||
% 9                                     forward differences ||
%============================================================%

%============================================================%
% 1                                    backward differences ||
% 2  2 - 1                                                  ||
% 3  3 - 2  3 - 2                                           ||
% 4  4 - 3  4 - 3  4 - 3                                    ||
% 5  5 - 4  5 - 4  5 - 4  5 - 4                             ||
% 6  6 - 5  6 - 5  6 - 5  6 - 5  6 - 5                      ||
% 7  7 - 6  7 - 6  7 - 6  7 - 6  7 - 6  7 - 6               ||
% 8  8 - 7  8 - 7  8 - 7  8 - 7  8 - 7  8 - 7  8 - 7        ||
% 9  9 - 8  9 - 8  9 - 8  9 - 8  9 - 8  9 - 8  9 - 8  9 - 8 ||
%============================================================%

%-------------------------------------------------------------------------%
forward_diff = f1';
for i = 2 : 1 : length(x1)
    diff = forward_diff( 2:(length(x1) - i + 2) , i - 1) - forward_diff( 1:(length(x1) - i + 1) , i - 1);
    forward_diff(:, i) = vertcat(diff, zeros(length(x1) - length(diff), 1));
end


backward_diff = f1';
for i = 2 : 1 : length(x1)
    diff = backward_diff( i:length(x1) , i - 1) - backward_diff( (i - 1):(length(x1) - 1) , i - 1);
    backward_diff(:, i) = vertcat(zeros(length(x1) - length(diff), 1), diff);
end
%-------------------------------------------------------------------------%

switch Newton
    case 'forward'
        for j = 1 : 1 : length(x0)
            polynomial=f1(1);
            h = (x1(2) - x1(1));
            q=((x0(j) - x1(1))/h);
            fact = 1;
            for i = 1 : 1 : n
                fact = fact * i;
% % %                 finite_diff = 0;
% % %                 for v = 0 : 1 : i
% % %                     N = i;
% % %                     C = factorial(N)/(factorial(v)*factorial(N - v));
% % %                     finite_diff = finite_diff + (-1)^(v) * C *f1(1 + i - v);
% % %                 end
                polynomial = polynomial + q/fact * forward_diff(1, i + 1);
                
                q = q * (((x0(j) - x1(1))/h) - (i + 1) + 1);
            end
            f(j) = polynomial;
        end
    case 'backward'
        for j = 1 : 1 : length(x0)
            polynomial=f1(n+1);  
            h = (x1(n+1)-x1(n));
            q=(x0(j)-x1(length(x1)))/h;
            fact = 1;
            for i = 1 : 1 : n
                fact = fact * i;
% % %                 finite_diff = 0;
% % %                 for v = 0 : 1 : i
% % %                     N = i;
% % %                     C = factorial(N)/(factorial(v)*factorial(N - v));
% % %                     finite_diff = finite_diff + (-1)^(v)*C*f1((n + 1) - v);
% % %                 end
                polynomial = polynomial + q/fact * backward_diff(length(x1), i + 1);
                h = (x1(n+1)-x1(n));
                q = q * ((x0(j)-x1(length(x1)))/h + i + 1 - 1);
            end
            f(j)=polynomial;
        end
end
end



function [f, x0] = cubic_spline(f1, x1, x0)
% h as variable, not vector
shift = 1;
n = length(x1) - 1;
f = zeros(length(x0), 1);
c = 0;
c(n + 1) = 0;
K = zeros(length(x1), 1);
L = zeros(length(x1), 1);

for k = 2 : 1 : n
    h_k =  x1(k + shift) - x1(k - 1 + shift);
    h_k_1 =  x1(k - 1 + shift) - x1(k - 2 + shift);

%     h(k) =  x1(k + shift) - x1(k - 1 + shift);
%     h(k - 1) =  x1(k - 1 + shift) - x1(k - 2 + shift);
    F = 3 * ((f1(k + 1) - f1(k))/h_k - (f1(k) - f1(k - 1))/h_k_1);
    V = 2 * (h_k + h_k_1);
    K(k) = (F - h_k_1 * K(k - 1))/(V - h_k_1 * L(k - 1));
    L(k) = h_k/(V - h_k_1 * L(k - 1));
end

for k = n : -1 : 2
    c(k) = K(k) - L(k) * c(k + 1);
end

for k = 1 : 1 : n
    h_k = x1(k + shift) - x1(k - 1 + shift);
    d(k) = (c(k + 1) - c(k))/(3 * h_k);
    b(k) = (f1(k + shift) - f1(k - 1 + shift))/h_k - c(k) * h_k - d(k) * (h_k)^2;
    a(k) = f1(k - 1 + shift);
end

k = 1;
for i = 1 : 1 : length(x0)
    if(x0(i) > x1(k + 1))
            k = k + 1;
    end
    h_k = x0(i) - x1(k);
    f(i)=a(k) + b(k) * h_k + c(k) * h_k^2 + d(k) * h_k^3;  
end 
end





