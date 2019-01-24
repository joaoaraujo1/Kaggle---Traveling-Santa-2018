clear; close all; clc
%%
%
%   João Araújo 2018/2019
%   
%   Brute-force optimization Script for the Traveling Santa 2018 contest
%
%
%%

%% Load data and get initial cost
Data = csvread('cities.csv',1);
load('Path_opt')
%Get primes in order for memberC usage
load('Primes')


costInit = 0;
pathHist = nan(1,size(Path_test,1)-1);
for i = 1:length(Path_test)-1
        
    cost = dist(Data(Path_test(i)+1,2:3),Data(Path_test(i+1)+1,2:3)');
    if(rem(i,10) == 0 && ~ismembc(Path_test(i),PrimeSet))
        cost = cost + .1*cost;
    end
    
    costInit = costInit + cost;
    pathHist(i) = costInit;
    
end


%% Brute-force

fprintf('\n\nLocal Search Optimization: Brute-force\n')
costDecrease = 0;

neighbourhood_size = [6]; %number of candidates for permutation

for p = neighbourhood_size
    
    fprintf('Brute force. Size of the path = %d\n',p+1);
    tic
    
    for dir = 1:1
        
        if dir == 1
            i_arr = 1:size(Path_test,1)-(p+2);
            %i_arr = 1 + 10-ceil(median(1:p+1)):10:size(Path_test,1)-(p+1); %Get the 10's in the middle
        else
            i_arr = size(Path_test,1)-2:-1:2;
        end

        for i = i_arr
            
            if rem(i,5000) == 0
                disp(i)
            end

            array_of_interest = nan(p+3,1);
            for j = 0:size(array_of_interest,1)-1
                array_of_interest(j+1,1) = Path_test(i+j,1)+1;
            end

            originalCost = 0;
            for j = 1:size(array_of_interest,1)-1
                cost = dist(Data(array_of_interest(j),2:3),Data(array_of_interest(j+1),2:3)');
                if(rem(j + i - 1,10) == 0 && ~ismembc(Data(array_of_interest(j),1),PrimeSet))
                    cost = cost + .1*cost;
                end

                originalCost = originalCost + cost;

            end
            edgesComb = nan(factorial(p+1),length(array_of_interest));
            edgesComb(:,2:end-1) = perms(array_of_interest(2:end-1));
            edgesComb(:,1) = array_of_interest(1);
            edgesComb(:,end) = array_of_interest(end);
            newCost = Inf(size(edgesComb,1),1);
            
            for k = 1:size(edgesComb,1)
                newCost(k) = 0;
                
                for j = 1:size(array_of_interest,1)-1
                    cost = dist(Data(edgesComb(k,j),2:3),Data(edgesComb(k,j+1),2:3)');
                    if(rem(j + i - 1,10) == 0 && ~ismembc(Data(edgesComb(k,j),1),PrimeSet))
                        cost = cost + .1*cost;
                    end

                    newCost(k) = newCost(k) + cost;

                end

            end
            
            
            [best_cost,best_k] = min(newCost);
            if best_cost < originalCost
                Path_test(i:i+length(array_of_interest)-1,1) = Data(edgesComb(best_k,:),1);
                costDecrease = costDecrease + (originalCost - best_cost);
                fprintf('\nCurrent Improvement of %.2f\n',costDecrease);
                save('Path_opt','Path_test');
            end

        end
               
    end
    toc    
    
end

fprintf('Neighbourhood Swap: Cost decreased by %.2f! Estimated new cost: %.2f\n',costDecrease,costInit-costDecrease);

%% Final Plot

estPathHist = zeros(1,length(Path_test));
for i = 1:size(Path_test,1)-1
    
    cost = dist(Data(Path_test(i)+1,2:3),Data(Path_test(i+1)+1,2:3)');
    if(rem(i,10) == 0 && ~isprime(Path_test(i)))
        cost = cost + .1*cost;
    end
    
    estPathHist(i+1) = estPathHist(i) + cost;
    
end
plot(estPathHist)
hold on
plot(pathHist)

finalCost = max(estPathHist);


