function [rEssentiality]=Essentiality(model,rxnsToCheck,ActRxns,nonEssential,Essential)

% function to assign accuracy score as described in the NICEgame workflow
% Evangelia Vayena 2022

changeCobraSolver('cplex_direct');
model.ub(ismember(model.rxns,strcat('DM_',model.mets)))=0;
model.lb(ismember(model.rxns,strcat('DM_',model.mets)))=0;
d=1;
rEssentiality{d,1}='rxnsToCheck';
rEssentiality{d,2}='Alternative Num';
rEssentiality{d,3}='Alternative';
rEssentiality{d,4}='True Positive';
rEssentiality{d,5}='True Negative';
rEssentiality{d,6}='False Positive';
rEssentiality{d,7}='False Negative';
rEssentiality{d,8}='MCC';
d=d+1;
essThr = 0.1;
for i=1:length(rxnsToCheck)
    modeli=model;
    % unblock alt
    for k=1:size(ActRxns{i},1)
        alt=ActRxns{i}{k,1};
        modelk=modeli;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'F_',alt)))=50;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'R_',alt)))=50;
        modelk.var_ub(ismember(modelk.varNames, strcat( 'NF_',alt)))=50;
        modelk.var_lb(ismember(modelk.varNames, strcat( 'NF_',alt)))=-50;
               solt = solveTFAmodelCplex(modelk);

       % essentiality
        [grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(modelk, 'TFA', modelk.genes, 0, 0, 0, essThr);
        
        if any(isnan(grRateKO_genetfa)) %keep track of the NaN KO by comparing essential_genes_tfaNaN with essential_genes_tfa
            grRateKO_genetfaNaN = grRateKO_genetfa;
            essential_genes_tfaNaN = modelk.genes(grRateKO_genetfaNaN(:,1)<essThr*solt.val);
            yesgenetfaNaN = ismember(modelk.genes,essential_genes_tfaNaN);
        end
        
        grRateKO_genetfa(isnan(grRateKO_genetfa)) = 0; %by default NaN is considered an essential KO
        essential_genes_tfa = modelk.genes(grRateKO_genetfa(:,1)<essThr*solt.val);
        not_essential_genes_tfa = setdiff(modelk.genes,essential_genes_tfa);
        % True Positive
        TP = numel(find(ismember(not_essential_genes_tfa,nonEssential)));
        % True Negative
        TN = numel(find(ismember(essential_genes_tfa,Essential)));
        % False Positive
        FP = numel(find(ismember(not_essential_genes_tfa,Essential)));
       % False Negative
        FN = numel(find(ismember(essential_genes_tfa,nonEssential)));
        
        rEssentiality{d,1}=rxnsToCheck{i};
        rEssentiality{d,2}=k;
        rEssentiality{d,3}=alt;
        rEssentiality{d,4}= TP;
        rEssentiality{d,5}= TN;
        rEssentiality{d,6}= FP;
        rEssentiality{d,7}= FN;
        MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
        
        rEssentiality{d,8}=MCC;
        
        d=d+1;
    end
    d=d+1;
    save('./tmp_accuracy.mat','rEssentiality')
end
    save('./all_accuracy.mat','rEssentiality')

end
