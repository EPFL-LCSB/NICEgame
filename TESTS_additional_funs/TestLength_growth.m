function [resultLength]=TestLength_growth(ActRxns)

% function to assign size score as described in the NICEgame workflow
% Evangelia Vayena 2022

d=1;
resultLength{d,1}='Alternative Number';
resultLength{d,2}='Alternative';
resultLength{d,3}='Number of Reactions';
d=d+1;
for k=1:size(ActRxns,1)
    lengthscore=size(ActRxns{k,1},1)-1;
    resultLength{d,1}=k;
    resultLength{d,2}=ActRxns{k,1};
    resultLength{d,3}=lengthscore;
    d=d+1;
end
end
