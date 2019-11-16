clc;
clear variables;
close all force;


fnc = @(t) lab_diff_f(t);

% d - derivative order
% p - finite diff order
% [A, C, b, divider] = C_coeff(d, p, method)
% divider of whole finite diff == 1/divider * (factorial(d)/1)

[A, C, b, divider, d, p] = C_coeff(3, 4, "forward")

[A, C, b, divider, d, p] = C_coeff(3, 4, "centered")
% C = C(end:-1:1)
% 1/8
clc
[A, C, b, divider, d, p] = C_coeff(1, 4, "backward")

[str] = str_finite_diff(C, d, p, divider, "backward")

A * C







%{
eqtext = '$$F_n={1 \over \sqrt{5}}';
eqtext = [eqtext '\left[\left({1+\sqrt{5}\over 2}\right)^n -'];
eqtext = [eqtext '\left({1-\sqrt{5}\over 2}\right)^n\right]$$'];


fib = zeros(1, 12);
for i = 1:12
    fib(i) = (((1+sqrt(5))/2)^i - ((1-sqrt(5))/2)^i)/sqrt(5);
end
subplot(2, 1, 1)
plot(1:12, fib, 'k^-')
text(0.5, 125, eqtext, 'Interpreter', 'Latex', 'Fontsize', 16)
%}

% subplot(1,1, 1)
% plot([0, 4], [1, 1])



function [str] = str_finite_diff(C, d, p, divider, method)

switch method
    case "forward"
        if (d ~= 1)
            start = ['$$F^' int2str(d) '(x) = '];
        else
            start = ['$$F(x) = '];
        end
        if(abs(C(1)) ~= 1)
            str = ['{ ' int2str(C(1)) ' * F(x) + '];
        else
            str = ['{F(x) +'];
        end
        
        if(C(2) > 0)
            str = [str int2str(C(2)) ' * F(x + h) + '];
        elseif (C(2) < 0)
            str(end-1:end) = [];
            str = [str int2str(C(2)) ' * F(x + h) + '];
        end
        for i = 3 : 1 : size(C, 1)
            if(C(i) < -1e-11)
                str(end-1:end) = [];
            end
            if(C(i) > 1e-11 || C(i) < -1e-11)
                if(abs(C(i)) ~= 1)
                    str = [str int2str(C(i)) ' * F(x + ' int2str(i) - 1 ' * h) + '];
                elseif(C(i) == 1)
                    str = [str ' F(x + ' int2str(i) - 1 ' * h) + '];
                else
                    str = [str ' - F(x + ' int2str(i) - 1 ' * h) + '];
                end
            end
        end
        str(end-1:end) = [];
        if(d >= 1)
            str= [str '\over' int2str(divider) ' * h^' int2str(d) '}']
        else
            str = [str '\over' int2str(divider) '}']
        end
        str = [str '$$'];
        start = [start str];
        str = start;
        
    case "backward"
        C = C(end:-1:1)
        if (d ~= 1)
            start = ['$$F^' int2str(d) '(x) = '];
        else
            start = ['$$F(x) = '];
        end
        if(abs(C(1)) ~= 1)
            str = ['{ ' int2str(C(1)) ' * F(x) + '];
        else
            str = ['{F(x) +'];
        end
        
        if(C(2) > 0)
            str = [str int2str(C(2)) ' * F(x - h) + '];
        elseif (C(2) < 0)
            str(end-1:end) = [];
            str = [str int2str(C(2)) ' * F(x - h) + '];
        end
        for i = 3 : 1 : size(C, 1)
            if(C(i) < -1e-11)
                str(end-1:end) = [];
            end
            if(C(i) > 1e-11 || C(i) < -1e-11)
                if(abs(C(i)) ~= 1)
                    str = [str int2str(C(i)) ' * F(x - ' int2str(i) - 1 ' * h) + '];
                elseif(C(i) == 1)
                    str = [str ' F(x - ' int2str(i) - 1 ' * h) + '];
                else
                    str = [str ' - F(x - ' int2str(i) - 1 ' * h) + '];
                end
            end
        end
        str(end-1:end) = [];
        if(d >= 1)
            if(d == 1)
                str = [str '\over' int2str(divider) ' * h}'];
            else
                str = [str '\over' int2str(divider) ' * h^' int2str(d) '}'];
            end
        else
            str = [str '\over' int2str(divider) '}'];
        end
        str = [str '$$'];
        start = [start str];
        str = start;
    case "centered"
end
figure('Name', 'Approximation of derivative', 'NumberTitle', 'off', ...
    'Position', [200 300 1100 300]) % [left bottom width height]
clf
plot(1, 1);
title([int2str(d) ' order derivative approximated by ' int2str(p) ...
    ' order ' char(method) ' finite difference'])
ylim([0, 3]);
xlim([0, 8]);
text(0.5, 1.5, str, 'Interpreter', 'Latex', 'Fontsize', 16)

end


function [f] = lab_diff_f(t)
v = 4 % [Hz]
omega = 2 * pi * v; % [rad/s]
f = 1.16 * t + 0.13 * sin(omega * t) - 0.89 * t^2;
end

function [df] = lab_diff_df(t)
% omega = 2 * pi * v;
omega = 4; % 4 Hz
df = 1.16 + 0.13 * omega * cos(omega * t) - 1.78 * t;
end

function [df, t] = lab_diff_do(fnc, a, b, n, k, method)
h = (b - a)/n;
t = a : h : b;
switch k 
    case 1
        switch method
            case "forward"
                df = (fnc(t + h) - fnc(t))/h;
            case "backward"
                df = (fnc(t) - fnc(t - h))/h;
            case "central"
                % Nope
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 2
        switch method
            case "forward"
                df = ( -fnc(t + 2*h) + 4 * fnc(t + h) - 3 * fnc(t))/(2 * h);
            case "backward"
                df = ( 3 * fnc(t) - 4 * fnc(t - h) + fnc(t - 2*h))/(2 * h);
            case "central"
                df = (fnc(t + h) - fnc(t - h))/(2 * h);
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 4
        switch method
            case "forward"
%                 df = 
            case "backward"
%                 df = 
            case "central"
                df = (-fnc(t + 2*h) + 8 * fnc(t + h) - 8 * fnc(t - h) + fnc(t - 2 * h))/(12 * h);
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 6
        switch method
            case "forward"
%                 df = 
            case "backward"
%                 df = 
            case "central"
%                 df = 
            otherwise
                fprintf("Error occured while entering method's name.");
        end
end

end

function [A, C, b, divider, d, p] = C_coeff(d, p, method)
format rat;
switch method
    case "forward"
        i_min = 0;
        i_max = d + p - 1;
        n = size(i_min : 1 : i_max, 2);
        i = (0 : 1 : n);
        for t = 1 : 1 : n
            A(t, :) =  [(i_min : i_max).^i(t)];
        end
        b = zeros(n, 1);
        b(d + 1) = 1;
        C = linsolve(A, b);
        [N, D] = rat(C);
        divider = max(D);
        C = C * divider;
%         A = [(i_min : i_max).^0; (i_min : i_max).^1; (i_min : i_max).^2; (i_min : i_max).^3]

    case "backward"
        i_min = -(d + p - 1);
        i_max = 0;
        n = size(i_min : 1 : i_max, 2);
        i = (0 : 1 : n);
        for t = 1 : 1 : n
            A(t, :) =  [(i_min : i_max).^i(t)];
        end
        b = zeros(n, 1);
        b(d + 1) = 1;
        C = linsolve(A, b);
        [N, D] = rat(C);
        divider = max(D);
        C = C * divider;
    case "centered"
        i_min = -fix((d + p - 1)/2); % i_max = -i_min
        i_max = fix((d + p - 1)/2);
        n = size(i_min : 1 : i_max, 2);
        i = (0 : 1 : n);
        for t = 1 : 1 : n
            A(t, :) =  [(i_min : 1 : i_max).^i(t)];
        end
        b = zeros(n, 1);
        b(d + 1) = 1;
        C = linsolve(A, b);
        C = linsolve(A, b);
        [N, D] = rat(C);
        divider = max(D);
        C = C * divider;
    otherwise
        fprintf("You have chosen a nonexistent method!");
end
format short;
end

