close all; clear; clc;
%%
%
%   João Araújo 2018/2019
%   
%   Permutations for OPT scripts for the Traveling Santa 2018 contest
%
%
%%

%% Load data
%Load cities and coordinates
Data = csvread('cities.csv',1);


%% Define parameters
sequenceLen = 2;                                                            % This is the K in K-opt Set to 2 or 4 to use with my code
neighbourN = 4;                                                             % How many neighbours we want to use for the permutations
distanceThreshold = 10;                                                     % What is the maximum distance between a node and its neighbours. Set to 0 if irrelevant

%% Create our Permutations matrix using KNN
%Build a KDTree searcher
kdt = KDTreeSearcher(Data(:,2:3));
% Allocate memory and initialize our permutations matrix
ptrLen = length(nchoosek(1:neighbourN,sequenceLen));                
Comparisons = zeros(ptrLen*(length(Data)-1),sequenceLen);
ptr = 1;
for i = 2:length(Data)
    
    if rem(i,1000) == 0
        fprintf('%d/%d\n',i,length(Data)-1);
    end
    
    if distanceThreshold == 0
        IDX = knnsearch(kdt,Data(i,2:3),'K',neighbourN+2); %IDX: index
        Neighbours = Data(IDX(IDX ~= 1),1)';
        Neighbours = Neighbours(1:neighbourN+1);
    else
        [IDX,D] = knnsearch(kdt,Data(i,2:3),'K',100); %IDX: index; D: euclidean distance
        Neighbours = [Data(i,1);Data(IDX(D ~= 0 & D <= distanceThreshold & IDX ~= 1),1)]';
        if(length(Neighbours) > neighbourN+1)
            Neighbours = [Data(i,1);Data(Neighbours(2:neighbourN+1),1)]';
        elseif(length(Neighbours) < neighbourN+1)
            firstN = IDX(D ~= 0 & IDX ~= 1);
            Neighbours = [Data(i,1);Data(firstN(1:neighbourN),1)]';
        end
    end
    
    % Get all the possible permutations for these neighbours
    combK = nchoosek(Neighbours,sequenceLen);
    Comparisons(ptr:ptr+size(combK,1)-1,:) = combK;
    ptr = ptr + size(combK,1);

end

%Remove duplicates and we are done
fprintf('Removing duplicates...\n');
Comparisons = unique(Comparisons,'rows');
save('Dups','Comparisons');
fprintf('Done.\n')