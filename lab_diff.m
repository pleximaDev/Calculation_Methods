clc;
clear variables;
close all force;


fnc = @(t) lab_diff_f(t);

% d - derivative order
% p - finite diff order
% [A, C, b, divider] = C_coeff(d, p, method)


[A, C, b, divider] = C_coeff(3, 1, "forward")
[A, C, b, divider] = C_coeff(3, 2, "centered")



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

function [A, C, b, divider] = C_coeff(d, p, method)
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
        b(d+1) = 1;
        C = linsolve(A, b);
        [N, D] = rat(C)
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
        b(d+1) = 1;
        C = linsolve(A, b);
        [N, D] = rat(C)
        divider = max(D);
        C = C * divider;
    case "centered"
        i_min = -(d + p - 1)/2; % i_max = -i_min
        i_max = (d + p - 1)/2;
        n = size(i_min : 1 : i_max, 2);
        i = (0 : 1 : n);
        for t = 1 : 1 : n
            A(t, :) =  [(i_min : 1 : i_max).^i(t)];
        end
        b = zeros(n, 1);
        b(d+1) = 1;
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

