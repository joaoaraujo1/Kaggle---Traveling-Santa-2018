close all; clear; clc;

%%
%
%   João Araújo 2018/2019
%   
%   Tour visualization Script for the Traveling Santa 2018 contest
%
%
%%

%% Tour drawing
Data = csvread('cities.csv',1);
prime_vector = primes(size(Data,1)-1);
prime_idx = prime_vector +1;
nonprime_vector = Data(~ismember(Data(:,1),prime_vector));
nonprime_idx = nonprime_vector +1;
prime_percent = length(prime_vector)/length(nonprime_vector);
hold on
load('Path_opt_kern');
h = plot(Data(nonprime_idx,2),Data(nonprime_idx,3),'b.','MarkerSize',2);
set(gcf, 'Position', get(0, 'Screensize'));
hold on
plot(Data(Path_test +1,2),Data(Path_test +1,3),'b');
plot(Data(prime_idx,2),Data(prime_idx,3),'g.','MarkerSize',5)
plot(Data(1,2),Data(1,3),'r.','MarkerSize',15);

%% Penalty edges drawing

costInit = 0;
penalties = zeros(1,size(Path_test,1)-1);
for i = 1:length(Path_test)-1
        
    cost = dist(Data(Path_test(i)+1,2:3),Data(Path_test(i+1)+1,2:3)');
    if(rem(i,10) == 0 && ~isprime(Path_test(i-1)))
        cost = cost + .1*cost;
        penalties(i) = .1*cost;
    else
        penalties(i) = 0;
    end
    
end

for i = 1:length(Path_test)-1
    if(penalties(i) > 0)
        plot(Data(Path_test(i:i+1)+1,2),Data(Path_test(i:i+1)+1,3),'r')
    end
end
legend('Non-prime cities','Edges','Prime cities','Start/End','Penalty Edge')
