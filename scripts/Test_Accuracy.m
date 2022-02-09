% script to calcualte the Accuracy score as described in the NICEgame workflow
% Evangelia Vayena 2022

clear;clc;
%%
load ( )%load model with thermodynamics with defined media ttmodel

addpath(genpath('/yorPath/mattfa/'))
addpath(genpath('/yorPath/IBM/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_osx'));
changeCobraSolver('cplex_direct','LP');
solt = solveTFAmodelCplex(ttmodel);

load('./Essential.mat') % load experimentally observed essential genes in a cell structure {'b0025';'b0029';...}
load('./nonEssential.mat')% load experimentally observed non essential genes in a cell structure {'b0038';'b0040';...}

% gene essentiality for original network
essThr = 0.1;
[grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(ttmodel, 'TFA', ttmodel.genes, 0, 0, 0, essThr);

if any(isnan(grRateKO_genetfa)) %keep track of the NaN KO by comparing essential_genes_tfaNaN with essential_genes_tfa
    grRateKO_genetfaNaN = grRateKO_genetfa;
    essential_genes_tfaNaN = ttmodel.genes(grRateKO_genetfaNaN(:,1)<essThr*solt.val);
    yesgenetfaNaN = ismember(ttmodel.genes,essential_genes_tfaNaN);
end

grRateKO_genetfa(isnan(grRateKO_genetfa)) = 0; %by default NaN is considered an essential KO
essential_genes_tfa = ttmodel.genes(grRateKO_genetfa(:,1)<essThr*solt.val);
yesgenetfa = ismember(ttmodel.genes,essential_genes_tfa);
not_essential_genes_tfa = setdiff(ttmodel.genes,essential_genes_tfa);
%% True Positive
TPwt = numel(find(ismember(not_essential_genes_tfa,nonEssential)));
TPwt1 = not_essential_genes_tfa(ismember(not_essential_genes_tfa,nonEssential));
%% True Negative
TNwt = numel(find(ismember(essential_genes_tfa,Essential)));
TNwt1 = essential_genes_tfa(ismember(essential_genes_tfa,Essential));
%% False Positive
FPwt = numel(find(ismember(not_essential_genes_tfa,Essential)));
FPwt1 = not_essential_genes_tfa(ismember(not_essential_genes_tfa,Essential));
%% False Negative
FNwt = numel(find(ismember(essential_genes_tfa,nonEssential)));
FNwt1 = essential_genes_tfa(ismember(essential_genes_tfa,nonEssential));
%% MCC
MCCwt = ((TPwt*TNwt)-(FPwt*FNwt))/sqrt((TPwt+FPwt)*(TPwt+FNwt)*(TNwt+FPwt)*(TNwt+FNwt));
ACC = (TPwt+TNwt)/(TPwt+TNwt+FPwt+FNwt);
%% load alternatives
load('') %load merged model with thermodynamics with defined media 
model = GFmodelThermo;
model.rxnGeneMat = tmodel.rxnGeneMat;
load('') % load the database

ActRxns1={};
for i=1:length(rxnsToCheck)
    load(['gapFillingBiomass_' GFmodelThermo.id '_' rxnsToCheck{i} '.mat'],'ActRxns')
    ActRxns1=[ActRxns1; ActRxns] ;
end
ActRxns=ActRxns1;
%% block all the reactions from the database
model = addNetFluxVariables(model);
indNF = getAllVar(model,{'NF'});
fromDB=(find(ismember(model.rxnFrom,{'DBmodel_1'})));
modelk.var_lb(ismember(model.varNames, strcat( 'NF_',model.rxns(fromDB))))=0;
model.var_ub(ismember(model.varNames, strcat( 'NF_',model.rxns(fromDB))))=0;
model.var_ub(ismember(model.varNames, strcat( 'F_',model.rxns(fromDB))))=0;
model.var_ub(ismember(model.varNames, strcat( 'R_',model.rxns(fromDB))))=0;
model.ub(fromDB)=0;
model.lb(fromDB)=0;
%% Accuracy test
[rEssentiality]=Essentiality(model,rxnsToCheck,ActRxns,nonEssential,Essential);