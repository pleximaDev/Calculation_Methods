
                        %%% The Beginning %%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


clc;
clear variables;
close all force;
addpath('./frwk')

% % % Testovaya = randn(3,3);
% % % testy=ones(3, 1)
% % % test=horzcat(Testovaya, testy)


                               %%% 2 %%%  

C = [1 -0.2589 -0.3093; -0.2589 1 -0.2705; -0.3093 -0.2705 1];
c = [1; 1 ; 1];
[x, ok]=my_gauss(C, c);
[x, ok]=my_gauss_jordan(C, c);
[x, ok]=my_Cramer(C, c);
[x, ok]=my_Invertible_matrix_A(C, c);
[x, ok]=my_chol(C, c);
% % test=chol(C, 'lower')
% % testY = test^(-1) * c
% % testX = (test')^(-1) * testY

                              %%% // %%%

                             
                            
                               %%% 3 %%%                             
                
              
K = 16;                             
b=randn(K, 1);

% 1) Matrica A poloj. opred. s dominir diag elem
A0 = randn(K);
I = eye(K, K);
A0 = tril(A0);
D=cell(4,1);

A = A0 * A0' + 5 * K * I;

% D(:, :,1) = horzcat(A, b); %semnadcatim stolbcom zapicivaetcya vektor b
D{1}={A,b};

% 2) matrica A simmietrichnaya, otricatelno opredelennaya, s domin diag elem
A0 = randn(K);
A0 = tril(A0);

A = A0 * A0' - 5 * K * I;

% D(:, :, 2) = horzcat(A, b);
D{2}={A,b};

% 3) matrica A Neumann - razrejennaya otric opred matreye

A = gallery('neumann', K);

% D(:, :, 3) = horzcat(A, b);
D{3}={A,b};
% 4) matrica A Neumamn razrejennaya otric opred matr v matlabe v vide POLNOY MATRICI(full)

A = full(A);

% D(:, :,4) = horzcat(A, b);
D{4}={A,b};

save('lab_slau_data.mat', 'D', '-v7');

clear variables;

load('lab_slau_data.mat');

                              %%% // %%%

                            
                            
                               %%% 4 %%%  
                               
      %%% Ocenka vremeni, zatrachivaemogo na reshenie zadachi%%%
   
% 4.1 & 4.2     
K=16
X=zeros(K, 5); 
T=zeros(4, 5);
for i = 1 : 1 : 4
% % A = D(:, :, i);     
% % b = A(:, 17);
% % A(:, 17) = [];
A = D{i}{1};
b = D{i}{2};
tic
[x, ok]=my_gauss(A, b);
t=toc;
% if ok == false
%     t = 0;
% end
T(i, 1)=t*ok; %ispravit na T(i, 1)=t*ok; chtobi izbejat uslovya
X(:, 1, i) = x;

tic
[x, ok]=my_gauss_jordan(A, b);
t=toc;
% if ok == false
%     t = 0;
% end
T(i, 2)=t*ok;
X(:, 2, i) = x;

tic
[x, ok]=my_Cramer(A, b);
t=toc;
% if ok == false
%     t = 0;
% end
T(i, 3)=t*ok;
X(:, 3, i) = x;

tic
[x, ok]=my_Invertible_matrix_A(A, b);
t=toc;
% if ok == false
%     t = 0;
% end
T(i, 4)=t*ok;
X(:, 4, i) = x;

tic
[x, ok]=my_chol(A, b);
t=toc;
% if ok == false
%     t = 0;
% end
T(i, 5)=t*ok;
X(:, 5, i) = x;
end
% % % % % % % % 

% 4.3 Opredelenie srednego vremeni

%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_gauss(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;
T=zeros(4, 5);

% % % % % 
% % % % % % % % % 
for i = 1 : 1 : 4
%     A = D(:, :, i);     
%     b = A(:, 17);    bilo
%     A(:, 17) = [];
    A = D{i}{1};
    b = D{i}{2};
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok]=my_gauss(A, b);
    timeVector(j, 1)=toc*ok;
    
    if ~ok 
        
          %vremya vipolneniya ==0
        break
    end   
%             tic
%             [x, ok]=my_gauss(A, b);
    end
    T(i, 1)= mean(timeVector); %zapicivaem srednee znachenie v T


  
%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_gauss_jordan(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;

%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
    A = D{i}{1};
    b = D{i}{2};
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok]=my_gauss_jordan(A, b);
    timeVector(j, 1)=toc*ok;
    if ~ok 
        break
    end   
     %zapicivaem srednee znachenie v T
    end
    T(i, 2)=mean(timeVector);


% % % %%%%%%%%%%%%%%%%%%%%% [x, ok]=my_gauss_jordan(A, b); %%%%%%%%%%%%%%%%%%%%%
% % % N = 1000;
% % % for i = 1 : 1 : 4
% % %     A = D(:, :, i);     
% % %     b = A(:, 17);
% % %     A(:, 17) = [];
% % %     timeVector=zeros(N, 1);
% % %     [x, ok]=my_gauss_jordan(A, b);
% % %     if ok == false
% % %         T(i, 2)=0;  %vremya vipolneniya ==0
% % %     end
% % %     if ok == true
% % %     for j = 1 : 1 : N
% % %     tic
% % %     [x, ok]=my_gauss_jordan(A, b);
% % %     timeVector(j, 1)=toc;
% % %     end
% % %     T(i, 2)=mean(timeVector); %zapicivaem srednee znachenie v T
% % %     end
% % % end

 

%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Cramer(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;

%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
    A = D{i}{1};
    b = D{i}{2};
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok]=my_Cramer(A, b);
    timeVector(j, 1)=toc*ok;
    if ~ok   %vremya vipolneniya ==0
        break
    end   
            
        
     %zapicivaem srednee znachenie v T
    
    end
    T(i, 3)=mean(timeVector);

% %%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Cramer(A, b); %%%%%%%%%%%%%%%%%%%%%
% N = 1000;
% for i = 1 : 1 : 4
%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
%     timeVector=zeros(N, 1);
%     [x, ok]=my_Cramer(A, b);
%     if ok == false
%         T(i, 3)=0;  %vremya vipolneniya ==0
%     end
%     if ok == true
%         for j = 1 : 1 : N
%             tic
%             [x, ok]=my_Cramer(A, b);
%             timeVector(j, 1)=toc;
%         end
%     T(i, 3)=mean(timeVector); %zapicivaem srednee znachenie v T
%     end
% end

%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Invertible_matrix_A(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;

%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
    timeVector=zeros(N, 1);
    for j = 1 : 1 : N
        tic
    [x, ok]=my_Invertible_matrix_A(A, b);
    timeVector(j, 1)=toc*ok;
     if ~ok 
%         T(i, 4)=0;  %vremya vipolneniya ==0
        break
    end   
            
        
    %T(i, 4)=mean(timeVector); %zapicivaem srednee znachenie v T
    end
    T(i, 4)=mean(timeVector);


% % % % % % % % % % % 
% %%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Invertible_matrix_A(A, b); %%%%%%%%%%%%%%%%%%%%%
% N = 1000;
% for i = 1 : 1 : 4
%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
%     timeVector=zeros(N, 1);
%     [x, ok]=my_Invertible_matrix_A(A, b);
%     if ok == false
%         T(i, 4)=0;  %vremya vipolneniya ==0
%     end
%     if ok == true
%         for j = 1 : 1 : N
%             tic
%             [x, ok]=my_Invertible_matrix_A(A, b);
%             timeVector(j, 1)=toc;
%         end
%     T(i, 4)=mean(timeVector); %zapicivaem srednee znachenie v T
%     end
% end


%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_chol(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;

%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
    A = D{i}{1};
    b = D{i}{2};
    timeVector=zeros(N, 1);
    
   
    for j = 1 : 1 : N
        tic
    [x, ok]=my_chol(A, b);
    timeVector(j, 1)=toc*ok;
    if ~ok 
%         T(i, 5)=0;  %vremya vipolneniya ==0
        break
    end   
            
     %zapicivaem srednee znachenie v T
    
    end
    T(i, 5)=mean(timeVector);
end

%среднее значение считать за фором, условия  timeVector(j, 1)=toc*ok;
%считать до условия





%%% Predidushaya versiya scheta vremeni bez if ~ok (but works)
%{
for i = 1 : 1 : 4
    A = D(:, :, i);     
    b = A(:, 17);
    A(:, 17) = [];
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok]=my_gauss(A, b);
    if ok == false
        T(i, 1)=0;  %vremya vipolneniya ==0
    end
    if ok == true
        
%             tic
%             [x, ok]=my_gauss(A, b);
            
            timeVector(j, 1)=toc;
            
    
    T(i, 1)= mean(timeVector); %zapicivaem srednee znachenie v T
    end
    end
end

  
%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_gauss_jordan(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;
for i = 1 : 1 : 4
    A = D(:, :, i);     
    b = A(:, 17);
    A(:, 17) = [];
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok]=my_gauss_jordan(A, b);
    if ok == false
        T(i, 2)=0;  %vremya vipolneniya ==0
    end
    if ok == true
    
    
    timeVector(j, 1)=toc;
    
    T(i, 2)=mean(timeVector); %zapicivaem srednee znachenie v T
    end
    end
end

% % % %%%%%%%%%%%%%%%%%%%%% [x, ok]=my_gauss_jordan(A, b); %%%%%%%%%%%%%%%%%%%%%
% % % N = 1000;
% % % for i = 1 : 1 : 4
% % %     A = D(:, :, i);     
% % %     b = A(:, 17);
% % %     A(:, 17) = [];
% % %     timeVector=zeros(N, 1);
% % %     [x, ok]=my_gauss_jordan(A, b);
% % %     if ok == false
% % %         T(i, 2)=0;  %vremya vipolneniya ==0
% % %     end
% % %     if ok == true
% % %     for j = 1 : 1 : N
% % %     tic
% % %     [x, ok]=my_gauss_jordan(A, b);
% % %     timeVector(j, 1)=toc;
% % %     end
% % %     T(i, 2)=mean(timeVector); %zapicivaem srednee znachenie v T
% % %     end
% % % end

 

%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Cramer(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;
for i = 1 : 1 : 4
    A = D(:, :, i);     
    b = A(:, 17);
    A(:, 17) = [];
    timeVector=zeros(N, 1);
    
    for j = 1 : 1 : N
        tic
    [x, ok]=my_Cramer(A, b);
    if ok == false
        T(i, 3)=0;  %vremya vipolneniya ==0
    end
    if ok == true
            timeVector(j, 1)=toc;
        
    T(i, 3)=mean(timeVector); %zapicivaem srednee znachenie v T
    end
    end
end
% %%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Cramer(A, b); %%%%%%%%%%%%%%%%%%%%%
% N = 1000;
% for i = 1 : 1 : 4
%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
%     timeVector=zeros(N, 1);
%     [x, ok]=my_Cramer(A, b);
%     if ok == false
%         T(i, 3)=0;  %vremya vipolneniya ==0
%     end
%     if ok == true
%         for j = 1 : 1 : N
%             tic
%             [x, ok]=my_Cramer(A, b);
%             timeVector(j, 1)=toc;
%         end
%     T(i, 3)=mean(timeVector); %zapicivaem srednee znachenie v T
%     end
% end

%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Invertible_matrix_A(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;
for i = 1 : 1 : 4
    A = D(:, :, i);     
    b = A(:, 17);
    A(:, 17) = [];
    timeVector=zeros(N, 1);
    for j = 1 : 1 : N
        tic
    [x, ok]=my_Invertible_matrix_A(A, b);
    if ok == false
        T(i, 4)=0;  %vremya vipolneniya ==0
    end
    if ok == true
            timeVector(j, 1)=toc;
        
    T(i, 4)=mean(timeVector); %zapicivaem srednee znachenie v T
    end
    end
end
% %%%%%%%%%%%%%%%%%%%%% [x, ok]=my_Invertible_matrix_A(A, b); %%%%%%%%%%%%%%%%%%%%%
% N = 1000;
% for i = 1 : 1 : 4
%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
%     timeVector=zeros(N, 1);
%     [x, ok]=my_Invertible_matrix_A(A, b);
%     if ok == false
%         T(i, 4)=0;  %vremya vipolneniya ==0
%     end
%     if ok == true
%         for j = 1 : 1 : N
%             tic
%             [x, ok]=my_Invertible_matrix_A(A, b);
%             timeVector(j, 1)=toc;
%         end
%     T(i, 4)=mean(timeVector); %zapicivaem srednee znachenie v T
%     end
% end


%%%%%%%%%%%%%%%%%%%%% [x, ok]=my_chol(A, b); %%%%%%%%%%%%%%%%%%%%%
N = 1000;
for i = 1 : 1 : 4
    A = D(:, :, i);     
    b = A(:, 17);
    A(:, 17) = [];
    timeVector=zeros(N, 1);
    
   
    for j = 1 : 1 : N
        tic
    [x, ok]=my_chol(A, b);
    
    if ok == false
        T(i, 5)=0;  %vremya vipolneniya ==0
    end
    if ok == true
            timeVector(j, 1)=toc;
   
    
    T(i, 5)=mean(timeVector); %zapicivaem srednee znachenie v T
    end
    end
end
%}
% % % % % T(2, 5)=0
% % % % % T(3, 5)=0
% % % % % T(4, 5)=0
% % % % % T(5, 5)=0

% %%%%%%%%%%%%%%%%%%%%% [x, ok]=my_chol(A, b); %%%%%%%%%%%%%%%%%%%%%
% N = 1000;
% for i = 1 : 1 : 4
%     A = D(:, :, i);     
%     b = A(:, 17);
%     A(:, 17) = [];
%     timeVector=zeros(N, 1);
%     [x, ok]=my_chol(A, b);
%     if ok == false
%         T(i, 5)=0;  %vremya vipolneniya ==0
%     end
%     if ok == true
%         for j = 1 : 1 : N
%             tic
%             [x, ok]=my_chol(A, b);
%             timeVector(j, 1)=toc;
%         end
%     
%     T(i, 5)=mean(timeVector); %zapicivaem srednee znachenie v T
%     end
% end



% % % % tic
% % % % [x, ok]=my_gauss_jordan(A, b);
% % % % t=toc
% % % % X(:, 2, i) = x
% % % % 
% % % % tic
% % % % [x, ok]=my_Cramer(A, b);
% % % % t=toc
% % % % X(:, 3, i) = x
% % % % 
% % % % tic
% % % % [x, ok]=my_Invertible_matrix_A(A, b);
% % % % t=toc
% % % % X(:, 4, i) = x
% % % % 
% % % % tic
% % % % [x, ok]=my_chol(A, b);
% % % % t=toc
% % % % X(:, 5, i) = x


                               
                            
                              %%% // %%%
                            
                              
                              
                              
                              
                              
                              
                              
                               %%% 5 %%%   
                              
              %%%%% Predstavlenie rezultatov raboty %%%%%                               
                 
              
for i = 1 : 1 : 4
% A = D(:, :, i);     
% b = A(:, 17);
% A(:, 17) = [];
A = D{i}{1};
b = D{i}{2};
fprintf('A    ');
fprintf('\r\n');
[m, n]=size(A);
fprintf('razmernost A -- %d x %d', m, n);
fprintf('\r\n');
disp(A);
fprintf('b    ');
fprintf('\r\n');
n=length(b);
fprintf('razmersnost b -- %d x 1', n);
fprintf('\r\n');
disp(b);
%methods('Gaussian', 'Gauss_Jordan', 'Cramer', 'Invertible_matrix_A', 'Cholesky');
    Gauss=X(:,1,i);
    Gauss_Jordan=X(:,2,i);
    Cramer=X(:,3,i);
    InverseA=X(:,4,i);
    Cholesky=X(:,5,i);
    Tablica=table(Gauss,Gauss_Jordan,Cramer,InverseA,Cholesky);
%Tablica=table(X(1, :, i), X(2, :, i), X(3, :, i), X(4, :, i), X(5, :, i), 'VariableNames', methods);
disp(Tablica);
end
 




title('Direct methods')
clf;
ax1 = subplot(1,2,1);
axtest=ax1;

kostil=bar(T'), ylabel('Time, ms'), title('Average time for a single computation');
% title('Average time for a single computation')
ax = gca;
ax.XTickLabel=({'Gauss','Gauss-Jordan',' Cramer',' Inverse ',' Cholesky'});
ylabel('Time, ms');
xlabel('Methods');
grid on
[m, n]=size(A);
lgd=legend(ax1,{'A>0, Symmetric','A<0, Symmetric','A<0, Sparse','A<0, Full Sparse'},'location','eastoutside');
title(lgd,'size of A is 16x16')
a=axes('position',get(gca,'position'),'visible','off');
%legend(a, kostil,{'size of A is 16x16'},'location','northeastoutside');

ax2 = subplot(1,2,2);
T2=T;
T2(:, 4)=[];
T2(3, :)=[];
kostil2=bar(T2'), ylabel('Time, ms'), xlabel('Methods');
title('Average time for a single computation  (matrix method excluded)')

ax = gca;
ax.XTickLabel=({'Gauss','Gauss-Jordan',' Cramer',' Cholesky'});
grid on
lgd2=legend(ax2,{'A>0, Symmetric','A<0, Symmetric','A<0, Full Sparse'},'location','eastoutside');
title(lgd2, 'size of A is 16x16')
a2=axes('position',get(gca,'position'),'visible','off');
% legend(a2, kostil2,{'size of A is 16x16'},'location','northeastoutside');



% % % % % % % % % left
% % % % % % % % 
% % % % % % % % figure(1);
% % % % % % % % clf;
% % % % % % % % subplot(1, 2,1)
% % % % % % % % bar(T);
% % % % % % % % ax-gca;
% % % % % % % % ax.XTicklabels=methods;
% % % % % % % % Ylabel('T, time');
% % % % % % % % names={'ad','kek',  'Neuman', 'full(Neumann)', 'lel'};
% % % % % % % % 
% % % % % % % % hold on;
% % % % % % % % subplot(1, 2, 1)=bar(T)
% % % % % % % % hold off;    
% % % % % % % % 
% % % % % % % % 
% % % % % % % % %right
% % % % % % % % subplot(1, 2,2)
% % % % % % % % bar(T);
% % % % % % % % ax-gca;
% % % % % % % % ax.XTicklabels=methods;
% % % % % % % % Ylabel('T, time');
% % % % % % % % names={'ad','kek',  'Neuman', 'full(Neumann)', 'lel'};
% % % % % % % % 
% % % % % % % % hold on;
% % % % % % % % subplot(1, 2, 1)=bar(T)
% % % % % % % % hold off;    




% % % % right 
% % % subplot(1, 2, 2)
% % % bar(T)
% % % ax=gca;
% % % ax.XTickLabels=methods;
% % % ylabel('sprava')
% % % legend(g(1, 1:4), 'Location', 'Eastoutside')
% % % rmpath('./frwk');
                              %%% // %%%
                              
   % na canva doljno bit' sleva                         

                                                     
                              
                              
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

                            %%% THE END %%%

