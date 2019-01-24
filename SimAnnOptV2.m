clear; close all; clc
%%
%
%   João Araújo 2018/2019
%   
%   Simulated Annealing 3.5 Script for the Traveling Santa 2018 contest
%
%
%%

%% Load files and
Data = csvread('cities.csv',1);
load('Path_opt')
load('Primes');
load('Neighbours')
Path_test_full = Path_test;
%% Parameters for fine-tuning

%Plots
showPlots = true;                      % display plots

%Simulated Annealing
Size = 200;                            % size of the subpath to be optimized
Inits = 5;                             % number of initializations
T = 5.0;                               % initial temperature
tk = .1;                               % cooling schedule
logCoolingSchedule = true;             % do we use log cooling schedule
logCoolingFactor = 100;                % larger = slower cooling
Mk = 500;                              % repetitions in each temperature
K = 20;                                % number of neighbour candidates for SA permutation
Runs = 3;                              % number of runs per initialization (> 1 is the equivalent of reheating)
randomStart = false;                   % start with a random solution

%Tabu list
listSize = 10;                         % size of tabu list for nodes (larger: encourages + exploration)
restrictionType = 2;                   % create a restriction list of 1- node pairs or 2-nodes

%3.5 opt
maxNeighbours = 10;                    % maximum neighbourhood size
optRounds = 1;                         % number of opt rounds
pathRevN = 2;                          % 1- original inbetween path; 2- 1+ flipped paths; 4- 2+ interchanged flipped paths



%%

currentOpt = 0;
for k = 1:Size-1:length(Path_test_full)
    
    fprintf('Optimizing cluster starting at index %d\n',k);
    
    %Get our subpath cost and assign it to our best cost so far
    Path_test = Path_test_full(k:k+Size-1);
    costInit = 0;
    for i = 1:length(Path_test)-1

        cost = dist(Data(Path_test(i)+1,2:3),Data(Path_test(i+1)+1,2:3)');
        if(rem(i+(k-1),10) == 0 && ~ismembc(Path_test(i),PrimeSet))
            cost = cost + .1*cost;
        end

        costInit = costInit + cost;
        
    end  
    bestCost = costInit;
    
    %Get the possible node permutations for 3+ opt
    Combinations = Inf(length(nchoosek(1:20,3))*(length(Path_test)-1),3);
    ptr = 1;
    for n = 1:length(Path_test)

        neighbours = NeighbourMatrix{Path_test(n)+1}(:,1);
        [~,nIdx] = ismember(neighbours,Path_test);
        neighbours = neighbours(nIdx ~= 0 & nIdx ~= 1 & nIdx ~= length(Path_test));
        if length(neighbours) > 3
            if length(neighbours) > maxNeighbours
                neighbours = neighbours(1:maxNeighbours);
            end

            combK = sort(nchoosek(neighbours,3),2);
            Combinations(ptr:ptr+size(combK,1)-1,:) = combK;
            ptr = ptr + size(combK,1);


        end
    end
    Combinations = unique(Combinations,'rows');
    Combinations(Combinations(:,1) == Inf,:) = [];
    
    for initialization = 1:Inits
        
        % Get a random initial subpath if desired
        if randomStart
            head = Path_test(1); tail = Path_test(end);
            Path_test_perm = Path_test(2:end-1);
            Path_test_perm = Path_test_perm(randperm(length(Path_test_perm)));
            Path_test = [Path_test(1);Path_test_perm;Path_test(end)];
        end
        
        % Show progression plots if desired
        if showPlots
            plot(Data((Path_test)+1,2),Data((Path_test)+1,3),'ro');
            hold on
            plot(Data((Path_test)+1,2),Data((Path_test)+1,3));
            title('Initialization');
            hold off
            drawnow
        end
        costHist = Inf(1,Runs);
        gainHist = Inf(1,Runs);
        
        if restrictionType == 1
            tabuList = Inf(listSize,2);
        else
            tabuList = Inf(listSize,1);
        end

        for r = 1:Runs

            kdtData = Data(Path_test+1,:);
            kdt = KDTreeSearcher(kdtData(:,2:3));  % kd-tree for SA permutations

            Tk = T;
            M_k = Mk;
            gainOpt = 0;

            %Iterate through all the temperature states
            Steps = 1;
            while true

                % and through all the iterations per temperature
                for m = 1:M_k

                    nodeA = Inf;
                    nodeB = Inf;

                    % Choose next swap based on the KDT and restricted by a
                    % tabu list
                    while true
                        
                        nodeA = datasample(Path_test(2:end-1),1);
                        a = find(Path_test == nodeA);

                        nodeNeighbourhood = kdtData(knnsearch(kdt,Data(nodeA+1,2:3),'K',K+1), 1:3);
                        nodeNeighbourhood(nodeNeighbourhood(:,1) == Path_test(end),:) = []; %Remove last node of the path from the neighbourhood
                        nodeNeighbourhood(nodeNeighbourhood(:,1) == Path_test(1),:) = [];   %Remove first node of the path from the neighbourhood
                        nodeB = datasample(nodeNeighbourhood(2:end,:),1);                   %Remove nodeA itself when sampling nodeB
                        
                        if restrictionType == 1 && ~(ismember([nodeA,nodeB(1,1)],tabuList,'rows') || ismember([nodeB(1,1),nodeA],tabuList,'rows'))
                            break
                        elseif restrictionType == 2 && ~(ismembc(nodeA,tabuList) || ismembc(nodeB(1,1),tabuList))
                            break
                        end

                    end

                    %Update tabu list
                    if restrictionType == 1 
                        tabuList = circshift(tabuList,1);
                        tabuList(1,:) = [nodeA,nodeB(1,1)];
                    else
                        tabuList = circshift(tabuList,2);
                        tabuList(1) = nodeA; tabuList(2) = nodeB(1,1);
                        tabuList = sort(tabuList);
                    end
                    
                    %Swap gain evaluation using KL partition algorithm

                    %Internal Cost of nodeA
                    first_edge = dist(Data(Path_test(a-1)+1,2:3),Data(Path_test(a)+1,2:3)');
                    if(rem(a-1+(k-1),10) == 0 && ~ismembc(Path_test(a-1),PrimeSet))
                        first_edge = first_edge + .1* first_edge;
                    end
                    Ia = first_edge;
                    second_edge = dist(Data(Path_test(a)+1,2:3),Data(Path_test(a+1)+1,2:3)');
                    if(rem(a+(k-1),10) == 0 && ~ismembc(Path_test(a),PrimeSet))
                        second_edge = second_edge + .1* second_edge;
                    end
                    Ia = Ia + second_edge;

                    %Internal Cost of nodeB
                    b = find(Path_test == nodeB(1,1));
                    first_edge = dist(Data(Path_test(b-1)+1,2:3),Data(Path_test(b)+1,2:3)');
                    if(rem(b-1+(k-1),10) == 0 && ~ismembc(Path_test(b-1),PrimeSet))
                        first_edge = first_edge + .1* first_edge;
                    end
                    Ib = first_edge;
                    second_edge = dist(Data(Path_test(b)+1,2:3),Data(Path_test(b+1)+1,2:3)');
                    if(rem(b+(k-1),10) == 0 && ~ismembc(Path_test(b),PrimeSet))
                        second_edge = second_edge + .1* second_edge;
                    end
                    Ib = Ib + second_edge;

                    if(a ~= b-1 && a~= b+1) %non adjacent nodes

                        %External Cost of nodeA
                        first_edge = dist(Data(Path_test(b-1)+1,2:3),Data(Path_test(a)+1,2:3)');
                        if(rem(b-1+(k-1),10) == 0 && ~ismembc(Path_test(b-1),PrimeSet))
                            first_edge = first_edge + .1* first_edge;
                        end
                        Ea = first_edge;
                        second_edge = dist(Data(Path_test(a)+1,2:3),Data(Path_test(b+1)+1,2:3)');
                        if(rem(b+(k-1),10) == 0 && ~ismembc(Path_test(a),PrimeSet))
                            second_edge = second_edge + .1* second_edge;
                        end
                        Ea = Ea + second_edge;

                        %External Cost of b
                        first_edge = dist(Data(Path_test(a-1)+1,2:3),Data(Path_test(b)+1,2:3)');
                        if(rem(a-1+(k-1),10) == 0 && ~ismembc(Path_test(a-1),PrimeSet))
                            first_edge = first_edge + .1* first_edge;
                        end
                        Eb = first_edge;
                        second_edge = dist(Data(Path_test(b)+1,2:3),Data(Path_test(a+1)+1,2:3)');
                        if(rem(a+(k-1),10) == 0 && ~ismembc(Path_test(b),PrimeSet))
                            second_edge = second_edge + .1* second_edge;
                        end
                        Eb = Eb + second_edge;

                        %Gain of nodeA swap
                        Ga = Ia - Ea;

                        %Gain of nodeB swap
                        Gb = Ib - Eb;

                        %Compute final Gain
                        Gain = Ga + Gb;
                        
                    else % Special cases: adjacent nodes

                        if a < b % x--(a-1)--a--b--(b+1)--w   >   x--(a-1)--b--a--(b+1)--w
                            minNode = a;
                            maxNode = b;
                        else     % x--(b-1)--b--a--(a+1)--w   >   x--(b-1)--a--b--(a+1)--w
                            minNode = b;
                            maxNode = a;
                        end

                        %simultaneous loop for normal and swap cost
                        C_normal = 0; C_swap = 0;
                        normalArr = [minNode-1 minNode maxNode maxNode+1];
                        swapArr   = [minNode-1 maxNode minNode maxNode+1];

                        for i = 1:3

                            costNormal = dist(Data(Path_test(normalArr(i))+1,2:3),Data(Path_test(normalArr(i+1))+1,2:3)');
                            costSwap = dist(Data(Path_test(swapArr(i))+1,2:3),Data(Path_test(swapArr(i+1))+1,2:3)');

                            if rem(minNode-2+i+(k-1),10) == 0 
                                %prime penalty for normal cost
                                if ~ismembc(Path_test(normalArr(i)),PrimeSet)
                                    costNormal = costNormal + .1*costNormal;
                                end
                                %prime penalty for swap cost
                                if ~ismembc(Path_test(swapArr(i)),PrimeSet)
                                    costSwap = costSwap + .1*costSwap;
                                end
                            end

                            C_normal = C_normal + costNormal;
                            C_swap = C_swap + costSwap;

                        end

                        %Gain is now the difference between both costs
                        Gain = C_normal - C_swap;

                    end


                    %Perform coin toss for annealing with p = exp(Gain/T)
                    Ntoss = 1000;
                    x = rand(1, Ntoss);
                    prob = (x < exp(Gain/Tk));
                    toss = datasample(prob,1);
                    
                    %Perform swap if we have gain or our coin toss was successful
                    if Gain >= 0 || toss == 1

                        tmp = Path_test(a);
                        Path_test(a) = Path_test(b);
                        Path_test(b) = tmp;
                        gainOpt = gainOpt + Gain;
                    end
                end
                
                if showPlots
                    plot(Data((Path_test)+1,2),Data((Path_test)+1,3),'ro');
                    hold on
                    plot(Data((Path_test)+1,2),Data((Path_test)+1,3));
                    title(['T = ' num2str(Tk)]);
                    hold off
                    drawnow
                end

                %if our Temperature is 0, we save our cost and proceed to the 3+opt stage
                if Tk == 0
                    
                    costA = 0;
                    for i = 1:size(Path_test,1)-1

                        cost = dist(Data(Path_test(i)+1,2:3),Data(Path_test(i+1)+1,2:3)');
                        if(rem(i+(k-1),10) == 0 && ~ismembc(Path_test(i),PrimeSet))
                            cost = cost + .1*cost;
                        end

                        costA = costA + cost;

                    end
                        
                    %% 3.5 Opt
                    if showPlots
                        title('3.5 opt');
                        drawnow
                    end
                    fprintf('I: %.2f | SA+: %.2f\n', costInit, costA);
                    
                    % Perform optimization
                    for rounds = 1:optRounds                                %across all optimization rounds
                        for forwback = 1:pathRevN                           %across all inbetween nodes path combinations
                            for c = 1:length(Combinations)                  %across all chosen combinations

                                if rem(c,1000) == 0
                                    fprintf('C = %d/%d\n',c,length(Combinations));
                                end

                                    [~,indices] = ismember(Combinations(c,:),Path_test);
                                    [nIdx,matchIdx] = sort(indices);
                                    chunk{1} = nIdx(1);
                                    if forwback == 1 || forwback == 4
                                        chunk{2} = (nIdx(1)+1:nIdx(2)-1);
                                    elseif forwback == 2 || forwback == 3
                                        chunk{2} = flip(nIdx(1)+1:nIdx(2)-1);
                                    end
                                    chunk{3} = nIdx(2);
                                    if forwback == 1 || forwback == 3
                                        chunk{4} = (nIdx(2)+1:nIdx(3)-1);
                                    elseif forwback == 2 || forwback == 4 
                                        chunk{4} = flip(nIdx(2)+1:nIdx(3)-1);
                                    end
                                    chunk{5} = nIdx(3);
                                    permNodes = repmat(num2cell(perms(nIdx)),2,1);
                                    paths{1} = chunk{2}; paths{2} = chunk{4};
                                    permPaths = repmat(perms(paths),6,1);
                                    finalPerms = cell(factorial(3)*2,5);
                                    finalPerms(:,[2 4]) = permPaths;
                                    finalPerms(:,[1 3 5]) = permNodes;
                                    allPerms = cell(size(finalPerms,1),size(finalPerms,2) - 1*isempty(chunk{2}) - 1*isempty(chunk{4}));

                                    for i = 1:length(finalPerms)
                                        allPerms(i,:) = finalPerms(i,~cellfun('isempty',finalPerms(i,:)));
                                    end

                                    Path_test_line = Path_test;

                                    for p = 1:length(allPerms)

                                        Path_test_line(nIdx(1):nIdx(3)) = Path_test_line(cell2mat(allPerms(p,:)));

                                        costOpt = 0;
                                        for i = 1:size(Path_test_line,1)-1

                                            cost = dist(Data(Path_test_line(i)+1,2:3),Data(Path_test_line(i+1)+1,2:3)');
                                            if(rem(i+(k-1),10) == 0 && ~ismembc(Path_test_line(i),PrimeSet))
                                                cost = cost + .1*cost;
                                            end

                                            costOpt = costOpt + cost;

                                        end

                                        if costOpt < costA - 0.0001
                                           Path_test = Path_test_line; 
                                           costA = costOpt;
                                           fprintf('I: %.2f | SA+: %.2f\n', costInit, costA);
                                           
                                           if showPlots
                                               plot(Data((Path_test)+1,2),Data((Path_test)+1,3),'ro');
                                               hold on
                                               plot(Data((Path_test)+1,2),Data((Path_test)+1,3));
                                               title('3.5 opt');
                                               hold off
                                               drawnow
                                           end
                                           
                                        end

                                    end


                            end
                        end
                    end


                    if costA < bestCost
                        currentOpt = currentOpt + (costInit - costA);
                        bestCost = costA;
                        Path_test_full(k:k+Size-1) = Path_test;
                        Path_test_copy = Path_test;
                        fprintf('Path Improved! Current improvement: %.2f\n',currentOpt);
                        save('Path_opt','Path_test_full');
                    end
                    break

                end

                %After each set of iterations lower temperature until T is roughly
                %equal to 0
                if ~logCoolingSchedule
                    if Tk - t_k > 0.0001
                        Tk = Tk - t_k;
                    else
                        Tk = 0;
                    end
                else
                    Tk = Tk - (1/log(Steps+1))/logCoolingFactor;
                    if Tk < 0.0001
                        Tk = 0;
                    end
                end

            end
        end
    end
end

























