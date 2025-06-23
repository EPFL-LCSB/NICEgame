function [resultTFBA]=TestTFBA_growth(model,ActRxns)
% function to assign TFA growth rate score as described in the NICEgame workflow
% Evangelia Vayena 2022

%block demand
model.var_ub(ismember(model.varNames,strcat('F_DM_',model.mets)))=0;
model.var_lb(ismember(model.varNames,strcat('R_DM_',model.mets)))=0;
d=1;
resultTFBA{d,1}='Alternative Number';
resultTFBA{d,2}='Alternative';
resultTFBA{d,3}='Growth rate TFA';
d=d+1;

% unblock alternative
for k=1:size(ActRxns,1)
    alt=ActRxns{k,1};
    modelk=model;
    modelk.var_ub(ismember(modelk.varNames, strcat( 'F_',alt)))=50;
    modelk.var_ub(ismember(modelk.varNames, strcat( 'R_',alt)))=50;
    modelk.var_ub(ismember(modelk.varNames, strcat( 'NF_',alt)))=50;
    modelk.var_lb(ismember(modelk.varNames, strcat( 'NF_',alt)))=-50;
    sol=solveTFAmodelCplex(modelk,300);
    resultTFBA{d,1}=k;
    resultTFBA{d,2}=alt;
    resultTFBA{d,3}=sol.val;
    
    d=d+1;
    if (mod(k,10)==0)
        save('./tmp_resultTFBA.mat',resultTFBA)
    end
end

save('./all_resultTFBA.mat',resultTFBA)
end

