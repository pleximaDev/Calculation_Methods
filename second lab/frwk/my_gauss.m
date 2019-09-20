
% Gaussian elimination (Method Gaussa)
function [x, ok]=my_gauss(A, b)
ok=false;
[n, m]=size(A);
A1=A;
x=zeros(n, 1);
if n == m
    ok=true;
end
% OriginMatrix = A;
% ExtendedMatrix = horzcat(A, b);
% if rank(OriginMatrix) == rank(ExtendedMatrix) % Rouché–Capelli theorem (Kronekera-Capelli)
%     ok = true;
if ok  %== true
    % pryamoy hod
    for i = 1 : 1 : n
        for k = i+1 : 1 : m
           
               coefficient=(A(k, i)/A(i, i));
               A(k,:)=A(k,:)-A(i,:)*(coefficient);
               b(k)=b(k)-b(i)*(coefficient);
% % % % %                Aki=A(k,i); % поменять коэффициент A(k, i) так, чтобы нормально считалось и для А и b
% % % % %                
% % % % %                A(k,:)=A(k,:)-A(i,:)*(Aki/A(i,i));
% % % % %                b(k)=b(k)-b(i)*(Aki/A(i,i));
               
%                A1(k,:)=A(k,:)-A(i,:)*(A(k,i)/A(i,i));
%                b(k)=b(k)-b(i)*(A(k, i)/A(i,i));
%                A(k,:)=A1(k,:);
               
               
        end
    end
    %obratny hod
    for i = n : -1 : 1
        Sum=0;
        for k = i+1 : 1 : n
            Sum=Sum + A(i,k)* x(k);
        end
    x(i)=(b(i)-Sum)/(A(i,i));
    end
end
if ok == false
    x=0;
end
end












% % function [x, ok]=my_gauss(A, b)
% % ok=false;
% % [n, m]=size(A);
% % x=zeros(n, 1);
% % epsilon = 1e-16
% % for i = 1 : 1 : n
% %     if (-epsilon) < sum(A(i, :)) < epsilon
% %         if (-epsilon) < b(i, 1) < epsilon
% %             ok = true
% %         end
% %     end
% % end
% % % OriginMatrix = A;
% % % ExtendedMatrix = horzcat(A, b);
% % % if rank(OriginMatrix) == rank(ExtendedMatrix) % Rouché–Capelli theorem (Kronekera-Capelli)
% % %     ok = true;
% % if ok == true
% %     % pryamoy hod
% %     for i = 1 : 1 : n
% %         for k = i+1 : 1 : m
% %                Aki=A(k,i);
% %                A(k,:)=A(k,:)-A(i,:)*(A(k,i)/A(i,i));
% %                b(k)=b(k)-b(i)*(Aki/A(i,i));
% %         end
% %     end
% %     %obratny hod
% %     for i = n : -1 : 1
% %         Sum=0;
% %         for k = i+1 : 1 : n
% %             Sum=Sum + A(i,k)* x(k);
% %         end
% %     x(i)=(b(i)-Sum)/(A(i,i));
% %     end
% % end
% % if ok == false
% %     x=0;
% % end
% % end


% working

% Gaussian elimination (Method Gaussa)

% % % function [x, ok]=my_gauss(A, b)
% % % ok=false;
% % % [n, m]=size(A);
% % % x=zeros(n, 1);
% % % OriginMatrix = A;
% % % ExtendedMatrix = horzcat(A, b);
% % % if rank(OriginMatrix) == rank(ExtendedMatrix) % Rouché–Capelli theorem (Kronekera-Capelli)
% % %     ok = true;
% % %     % pryamoy hod
% % %     for i = 1 : 1 : n
% % %         for k = i+1 : 1 : m
% % %                Aki=A(k,i);
% % %                A(k,:)=A(k,:)-A(i,:)*(A(k,i)/A(i,i));
% % %                b(k)=b(k)-b(i)*(Aki/A(i,i));
% % %         end
% % %     end
% % %     %obratny hod
% % %     for i = n : -1 : 1
% % %         Sum=0;
% % %         for k = i+1 : 1 : n
% % %             Sum=Sum + A(i,k)* x(k);
% % %         end
% % %     x(i)=(b(i)-Sum)/(A(i,i));
% % %     end
% % % end
% % % if ok == false
% % %     x=0;
% % % end
% % % end