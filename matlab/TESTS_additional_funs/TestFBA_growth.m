function [resultFBA]=TestFBA_growth(model,ActRxns)

% function to assign FBA growth rate score as described in the NICEgame workflow
% Evangelia Vayena 2022

model.ub(ismember(model.rxns,strcat('DM_',model.mets)))=0;
model.lb(ismember(model.rxns,strcat('DM_',model.mets)))=0;
d=1;
resultFBA{d,1}='Alternative Number';
resultFBA{d,2}='Alternative';
resultFBA{d,3}='Growth Rate';
d=d+1;
% unblock alt
for k=1:size(ActRxns,1)
    alt=ActRxns{k,1};
    modelk=model;
    modelk=changeRxnBounds(modelk,alt,50,'u');
    modelk=changeRxnBounds(modelk,alt,-50,'l');
    sol=solveFBAmodelCplex(modelk);
    resultFBA{d,1}=k;
    resultFBA{d,2}=alt;
    resultFBA{d,3}=sol.f;
    
    d=d+1;
    if (mod(k,10)==0)
        save('./tmp_resultFBA.mat',resultFBA)
    end
end

save('./all_resultFBA.mat',resultFBA)
end

