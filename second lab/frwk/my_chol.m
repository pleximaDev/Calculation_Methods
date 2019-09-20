% Cholesky decomposition


    function [x,ok] = my_chol(A,b)
   [n, m] = size(A);
   ok = false;
   x = zeros(n, 1);
   epsil = 1e-16;
   %[V, De]=eigs(A);

   if ((all(eigs(A,n)>epsil)))
   %    if (issymmetric(A)) && (all(eigs(A,n)>epsil))
       ok = true;

   %if A == (A')
%        UglovoyMinor = A(1, 1);
%        epsil = 1e-16;
%        for i = 1 : 1 : m % proverka na polojitelnuu-opredelennost
%             if UglovoyMinor > epsil
%                 ok = true;
%                UglovoyMinor = det(A(1:i, 1:i));
%             else
%                 ok = false;
%                 break
%             end
%        end
%        
% % e = 1;
% % while e ~= m
% % if De(e, e) > epsil
% %     ok = true;
% % else
% %     ok = false;
% %     return
% % end
% % e = e + 1;
% % end

% if diag(De) > epsil
%     ok = true;
% else
%     ok = false;
%     return
% end
   %end
  
       y = zeros(n, 1);
       L=zeros(n,m);
        for  i = 1:1:n 
            for  j = 1:1:(i-1) 
                SUML1 = 0; 
                for   k=1:1:(j-1)
                    SUML1 = SUML1+L(i,k)*L(j,k); 
                end 
                L(i,j)=(A(i,j)-SUML1)/(L(j,j)); 
            end 
            SUML2 = 0;
            for  k = 1:1:(i-1)
                SUML2 = SUML2 + L(i,k)^2; 
            end 
            L(i,i) = sqrt(A(i,i) - SUML2); 
        end
          for  i=1:1:length(y)
            s = 0; 
                for k = 1:1:(i-1)
                    s = s + L(i,k)*y(k); 
                end 
            y(i) = (b(i) - s)/L(i,i); 
          end
          Tr=L';
        for  i=n:-1:1
            s = 0; 
            for k = (i+1):1:n
               s = s + Tr(i,k)*x(k); 
            end 
            x(i) = (y(i) - s)/Tr(i,i); 
        end
   end
end












