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


clf;
subplot(1, 1, 1);
stairs(x0, f2, 'LineWidth', 2);
hold on;
grid on;
grid minor;
% plot(x0, f3);
% plot(x0, f4);
stem(x1, f1, 'LineWidth', 2);
plot(x0, f0, 'LineWidth', 2);

title('Interpolation');
legend('Nearest-neighbor interpolation','Experimental data model','Analytical model function');
legend('location','northeastoutside');
ylabel('f(x)');
xlabel('x');
hold off;

% subplot(2, 1, 2);
% plot(x0, f5);
% hold on;
% grid on;
% grid minor;
% plot(x0, f6);
% plot(x0, f7);
% stem(x1, f1);
% plot(x0, f0);
% legend('Newton polynomial forward','Newton polynomial backward','Cubic spline interpolation', 'Experimental data model','Analytical model function')
% legend('location','eastoutside');
% ylabel('f(x)');
% xlabel('x');
% hold off;




function f = interp_function(omega, x)
f = sin(omega * x);
end

function [f, x0] = nearest_neighbor(f1, x1, x0)
k = 1;
f = zeros(length(x0));
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





