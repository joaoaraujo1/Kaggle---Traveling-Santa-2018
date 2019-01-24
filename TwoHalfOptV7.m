%%
%
%   João Araújo 2018/2019
%   
%   2.5-Opt Script for the Traveling Santa 2018 contest
%
%
%%


clear; close all; clc

%% Preferences to manage computational resources

maxPerms = false;      % Explore all possible node and chunk permutations. Higher computational cost


%% Load data

%Load cities and initial Path
Data = csvread('cities.csv',1);
%Load current best path (Path_test array)
load('Path_opt_conc')
%Get list of primes
load('Primes')
%Get our city permutations to test
load('Duplets')

% Get initial path cost, precompute edges cost and penalties
costInit = 0;
penalties = zeros(1,size(Path_test,1)-1);
baseCostPath = nan(size(Path_test,1)-1,1);
for i = 1:length(Path_test)-1
        
    cost = dist(Data(Path_test(i)+1,2:3),Data(Path_test(i+1)+1,2:3)');
    baseCostPath(i) = cost;
    if(rem(i,10) == 0 && ~ismembc(Path_test(i),PrimeSet))
        cost = cost + .1*cost;
        penalties(i) = .1*cost;
    else
        penalties(i) = 0;
    end
    costInit = costInit + cost;
    
end

%% 2.5 opt
fprintf('2.5 Opt...\n\n');
currentOpt = 0;
for index = 1:length(Comparisons)
    
    if rem(index,1000) == 0
        el = toc;
        fprintf('Permutation: %d/%d | Opt: %.2f | ETA: %.2fh\n',index,length(Comparisons),currentOpt,el*(length(Comparisons)-index)/60/60/1000);
        tic
    end
    
    curNeighbours = Comparisons(index,:);

    [~,indices] = ismember(curNeighbours,Path_test);
    [nIdx,matchIdx] = sort(indices);
    head = nIdx(1)-1;
    tail = nIdx(end)+1;

    % Get our original path cost (precomputed): pure length + penalties
    bestCost = sum(baseCostPath(head:tail-1)) + sum(penalties(head:tail-1));

    %Get our smart permutations normal
    chunk{1} = indices(1);
    chunk{2} = nIdx(1)+1:nIdx(2)-1;
    chunk{3} = indices(2);
    if maxPerms
        finalPerms = perms(chunk);
    else
        finalPerms = chunk;
    end
    
    % Array formatting
    allPerms1 = cell(size(finalPerms,1),size(finalPerms,2) - 1*isempty(chunk{1}) - 1*isempty(chunk{2}) - 1*isempty(chunk{3})); 
    for i = 1:size(finalPerms,1)
        allPerms1(i,:) = finalPerms(i,~cellfun('isempty',finalPerms(i,:)));
    end
    
    
    %Get our inter-node path reversed
    chunk{2} = flip(nIdx(1)+1:nIdx(2)-1);
    if maxPerms
        finalPerms = perms(chunk);
    else
        finalPerms = chunk;
    end
    % Array formatting
    allPerms2 = cell(size(finalPerms,1),size(finalPerms,2) - 1*isempty(chunk{1}) - 1*isempty(chunk{2}) - 1*isempty(chunk{3})); 
    for i = 1:size(finalPerms,1)
        allPerms2(i,:) = finalPerms(i,~cellfun('isempty',finalPerms(i,:)));
    end
    
    costChunk = nan(size(allPerms1,1),2);
    
    for perm = 1:2 %reverse and normal
        if perm == 1
            allPerms = allPerms1;  
        else
            allPerms = allPerms2;    
        end
    
        for k = 1:size(allPerms,1)

            curChunk = [head,cell2mat(allPerms(k,:)),tail];

            baseCostPenalty = Inf(length(head:tail-1),1);
            ptr = 1;
            for i = 1:length(allPerms(k,:))+1

                % Check head - partition relationship
                if i == 1
                    if(head == allPerms{k,i}(1)-1)
                        baseCostPenalty(ptr) = baseCostPath(head);
                    else
                        baseCostPenalty(ptr) = dist(Data(Path_test(head)+1,2:3),Data(Path_test(allPerms{k,i}(1))+1,2:3)');
                    end
                    ptr = ptr + 1;

                %On the rest of the loop check length, add elements and
                %transition   
                else
                    if perm == 1
                        baseCostPenalty(ptr:ptr + length(allPerms{k,i-1}(1:end-1)) -1) = baseCostPath(allPerms{k,i-1}(1:end-1));
                    else
                        baseCostPenalty(ptr:ptr + length(allPerms{k,i-1}(1:end-1)) -1) = baseCostPath(allPerms{k,i-1}(2:end));
                    end
                    ptr = ptr + length(allPerms{k,i-1}(1:end-1));

                    if i <= length(allPerms(k,:))
                        baseCostPenalty(ptr) = dist(Data(Path_test(allPerms{k,i-1}(end))+1,2:3),Data(Path_test(allPerms{k,i}(1))+1,2:3)');
                        ptr = ptr + 1;

                    %Transition to tail
                    else
                        baseCostPenalty(ptr) = dist(Data(Path_test(tail)+1,2:3),Data(Path_test(allPerms{k,i-1}(end))+1,2:3)');
                        ptr = ptr + 1; 
                    end


                end
            end
            % Get full chunk cost with penalty
            penaltyIdx1 = (rem(head:tail-1,10) == 0)';                             % each 10th step
            penaltyIdx2 = ~ismembc(Path_test(curChunk(1:end-1)),PrimeSet);         % non prime
            penaltyIdx = (penaltyIdx1 == penaltyIdx2) & penaltyIdx1 ~= 0;          % combine conditions
            baseCostPenalty(penaltyIdx) = 1.1*baseCostPenalty(penaltyIdx);         % set penalty to 10%
            costChunk(k,perm) = sum(baseCostPenalty);                              % len + pen

        end
    end

    [minVal,minIdx] = min(costChunk);
    if any(minVal < bestCost -0.001)
        
        % Update current optimization value
        currentOpt = currentOpt + (bestCost - min(minVal));
        
        % Update best tour
        if maxPerms
            if(minVal(1) < minVal(2))
                Path_test(head:tail) = Path_test([head,cell2mat(allPerms1(minIdx(1),:)),tail]);
            else
                Path_test(head:tail) = Path_test([head,cell2mat(allPerms2(minIdx(2),:)),tail]);
            end
        else
            if minIdx == 1
                Path_test(head:tail) = Path_test([head,cell2mat(allPerms1),tail]);
            else
                Path_test(head:tail) = Path_test([head,cell2mat(allPerms2),tail]);
            end
        end
        
        % Save our new best tour in a mat file
        save('Path_opt','Path_test');

        % Recompute pure length / penalty costs
        for i = 1:length(Path_test)-1
            cost = dist(Data(Path_test(i)+1,2:3),Data(Path_test(i+1)+1,2:3)');
            baseCostPath(i) = cost;
            if(rem(i,10) == 0 && ~ismembc(Path_test(i),PrimeSet))
                penalties(i) = .1*cost;
            else
                penalties(i) = 0;
            end
        end
    end
    
end