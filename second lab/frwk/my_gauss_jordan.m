function [x, ok]=my_gauss_jordan(A, b)
[n, m]=size(A);
ok=false;
if n==m % proverka usloviya primenimosti
    ok=true;
% if ok == true % realizaciya scripta pri vipolnenii usloviya primenimosti
    x=zeros(n, 1);
    C=horzcat(A, b);
    % pryamoy hod
    for i = 1 : 1 : n
        for k = i+1 : 1 : m
            %Cki = C(k, i);
               C(k, :) = C(k, :) - C(i, :) * C(k, i)/C(i, i); % C(k, :) = C(k, :) - C(i, :) * C(k, i)/C(i, i);
        end
    end

    % obratny hod

    for i = n : -1 : 1
    for k = i-1 : -1 : 1
        %Cki = C(k, i);
        C(k, :) = C(k, :) - C(i, :) * C(k, i)/C(i, i); %         C(k, :) = C(k, :) - C(i, :) * Cki/C(i, i);  
    end
    end
    % formirovanie otveta
    x=C(:, n+1);
    for l = 1 : 1 : n
    x(l, 1)=x(l,1)/C(l, l);
    end
end
if ok == false % funkciya vozvratit nulevoy vektor tak, kak eto trebuetsya v zadanii
    x=zeros(n, 1);
end
end