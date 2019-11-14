clc;
clear variables;
close all force;


fnc = @(t) lab_diff_f(t);


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
                df = 
            case "backward"
                df = 
            case "central"
                df = (-fnc(t + 2*h) + 8 * fnc(t + h) - 8 * fnc(t - h) + fnc(t - 2 * h))/(12 * h);
            otherwise
                fprintf("Error occured while entering method's name.");
        end
    case 6
        switch method
            case "forward"
                df = 
            case "backward"
                df = 
            case "central"
                df = 
            otherwise
                fprintf("Error occured while entering method's name.");
        end
end

end
