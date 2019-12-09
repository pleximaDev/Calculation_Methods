clc;
clear variables;
close all force;



fnc = @(t) hw_int_func(t);
% fnc = @(t) line_(t);

a = 0.7;
b = 1.5;
eps = 1e-6;
n = [10, 100, 1000];


[var] = graph(fnc, a, b, n, 10);
% [var] = graph(fnc, -12, 8, n, 100);




[I_l] = hw_int_analog(fnc, a, b, n(2), "left Riemann sum")
[I_r] = hw_int_analog(fnc, a, b, n(2), "right Riemann sum")
[I] = hw_int_analog(fnc, a, b, n(2), "midpoint Riemann sum")
[I] = hw_int_analog(fnc, a, b, n(2), "trapezoidal Riemann sum")

x = a:(b - a)/n(2):b;
f = fnc(x);
trapz(x, f)


function [f] = line_(t)
% f = sin(4* pi *t);
sigma = sqrt(0.5);
mu = -2;
f = 1/(sigma * sqrt(2 * pi)) * exp(-1/2 * ((t - mu)/sigma).^2)
end

function [f] = hw_int_func(t)
v = 4; % [Hz]
omega = 2 * pi * v; % [rad/s]
f = 1.16 * t + 0.13 * sin(omega * t) - 0.89 * t.^2;
end


function [I] = hw_int_analog(f, a, b, n, method)
x = a:(b - a)/n:b;
% n = n + 1;
% f = fnc(x);
I = 0;
shift = 1;

switch method
    case "left Riemann sum"
        for i = 0 : 1 : n - 1 % for i = 1 : 1 : n - 1
             I = I + f(x(i + shift)) * (x(i + 1 + shift) - x(i + shift));
        end
    case "right Riemann sum"
        for i = 1 : 1 : n % for i = 2 : 1 : n
             I = I + f(x(i + shift)) * (x(i + shift) - x(i - 1 + shift));
        end
    case "midpoint Riemann sum"
        for i = 1 : 1 : n % for i = 2 : 1 : n
             I = I + f( (x(i - 1 + shift) + x(i + shift))/2 ) * (x(i + shift) - x(i - 1 + shift));
        end
    case "trapezoidal Riemann sum"
        for i = 0 : 1 : n - 1 %for i = 1 : 1 : n - 1
            I = I + (f(x(i + shift)) + f(x(i + 1 + shift)))/2 * (x(i + 1 + shift) - x(i + shift));
        end
% % % % %         sum = 0;
% % % % %         for i = 2 : 1 : n - 1
% % % % %             sum = sum + f(x(i)) * (x(i + 1) - x(i - 1))/2;
% % % % %         end
% % % % %         I = I + f(a)/2 * (x(1 + shift) - a) + f(b)/2 * (b - x(n - 1)) + sum;
    otherwise
        fprintf("You have chosen a nonexistent method!");
end
end


function [var] = graph(fnc, a, b, n, N)
var = 'variable';
figure(1)
clf
for i = 1 : 1 : 3
    subplot(3, 1, i)
    t = a:(b - a)/n(i):b;
    f = fnc(t);
    plot(t, f)
    grid on
    grid minor
    ylabel('f(x)');
    xlabel('x');
    text = ['n = ',int2str(n(i))];
    legend(text);
end


t = a:(b - a)/10000:b;
f = fnc(t);

l = a:(b - a)/N:b;

%______________left______________
left_x = a:(b - a)/100:b;
j = 1;
for i = 1 : 1 : length(l)
    while((i <= length(l) - 1) && left_x(j) ~= l(i + 1))
        left_f(j) = fnc(l(i));
        j = j + 1;
    end
    if(i == length(l))
        left_f(end + 1) = fnc(l(end));
    end
end
%________________________________


%______________right______________
right_x = a:(b - a)/100:b;
j = 1;
for i = 1 : 1 : length(l)
    while((i <= length(l) - 1) && right_x(j) ~= l(i + 1))
        right_f(j) = fnc(l(i + 1));
        j = j + 1;
    end
    if(i == length(l))
        right_f(end + 1) = fnc(l(end));
    end
end
%_________________________________



%______________midpoint______________
mid_x = a:(b - a)/100:b;
j = 1;
for i = 1 : 1 : length(l)
    while((i <= length(l) - 1) && right_x(j) ~= l(i + 1))
        mid_f(j) = fnc((l(i) + l(i + 1))/2);
        j = j + 1;
    end
    if(i == length(l))
        mid_f(end + 1) = fnc((l(end) + l(end - 1))/2);
    end
end
%____________________________________


figure(2)
clf
subplot(2,2, 1)
plot(t, f, 'LineWidth', 1.5, 'Color', 'r');
hold on
stem(l,fnc(l), 'LineWidth', 1.5); 
stairs(left_x, left_f, 'LineWidth', 1.5, 'Color', 'b');
hold off
grid on
grid minor
title("Left Riemann sum")
ylabel('f(x)');
xlabel('x');

subplot(2,2, 2)
plot(t, f, 'LineWidth', 1.5, 'Color', 'r');
hold on
stem(l,fnc(l), 'LineWidth', 1.5); 
stairs(right_x, right_f, 'LineWidth', 1.5, 'Color', 'b');
hold off
grid on
grid minor
title("Right Riemann sum")
ylabel('f(x)');
xlabel('x');

subplot(2,2, 3)
plot(t, f, 'LineWidth', 1.5, 'Color', 'r');
hold on
stem(l,fnc(l), 'LineWidth', 1.5); 
stairs(mid_x, mid_f, 'LineWidth', 1.5, 'Color', 'b');
hold off
grid on
grid minor
title("Midpoint Riemann sum")
ylabel('f(x)');
xlabel('x');
end







