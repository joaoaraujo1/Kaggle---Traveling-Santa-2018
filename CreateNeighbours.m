clear; close all; clc
%%
%
%   João Araújo 2018/2019
%   
%   Neighbourhood structure Script for the Traveling Santa 2018 contest
%
%
%%
Data = csvread('cities.csv',1);

%% Define the number of candidates for each node

numberOfNodes = 55;

%% Get candidate structure for each node
kdt = KDTreeSearcher(Data(:,2:3));
NeighbourMatrix{1} = NaN;
NeighbourMatrix{length(Path_test)} = NaN;

for i = 2:length(Path_test)-1
    if rem(i,1000) == 0
        disp(i)
    end
    NeighbourMatrix{i} = Data(knnsearch(kdt,Data(i,2:3),'K',numberOfNodes), 1:3);
    NeighbourMatrix{i} = NeighbourMatrix{i}(2:end,:); %Remove the same node
    NeighbourMatrix{i}(NeighbourMatrix{i}(:,1) == 0,:) = []; % Remove origin
end


save('Neighbours','NeighbourMatrix');