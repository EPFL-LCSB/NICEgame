
addpath(genpath('/yourPath/NICEgame/'))
addpath(genpath('/yourPath/cobratoolbox/'))

% load one of the DBs
load('/yourPath/databases/database.mat')
DBmodel = model;
% load the model
load('./yourPath/model.mat')
sourceModel = ;
organism_id = '';
% thermo? if yes flagTFA = 1, if no flagTFA = 0
flagTFA = 1;
thermo_data = yourthermoDB; %{} if no thermo

%% Check what are the biomass buidling blocks that cannot be porduced
obj = find(model.c);
bbbs_m = model.mets(find(model.S(:,find(model.c))<0));
r = length(model.rxns);
model = addDemandReaction(model,bbbs_m);
for i = 1:length(bbbs_m)
    modeli = model;
modeli.c = zeros(length(modeli.rxns),1);
modeli.c(r+i) = 1;
sol = optimizeCbModel(modeli,'max');
if sol.f < 10^-6
    sol.f = 0;
end
sol_bbb(i,1) = sol.f;
end
problematic_BBBs = bbbs_m(sol_bbb == 0); %
rxnID = strcat('DM_',problematic_BBBs);
%%
rmpath('/yourPath/cobratoolbox/')
addpath(genpath('/yourPath/mattfa/'))
addpath(genpath('/yourPath/IBM/ILOG/CPLEX_Studio1271/'))
changeCobraSolver('cplex_direct','LP')

% merge model with DB
[GFmodel, conflict] = PrepareForGapFilling(sourceModel, {DBmodel},'', 0,flagTFA,{},[],thermo_data);

%% define the media
% you can comment this part if you have already defined the media
EX = find(contains(GFmodel.varNames,'R_EX_'));
GFmodel.var_ub(EX) = 0;
media = {'EX_cl_e';'EX_ca2_e';'EX_cobalt2_e';'EX_mobd_e';...
    'EX_cu2_e';'EX_fe2_e';'EX_fe3_e';'EX_h2o_e';'EX_k_e';'EX_mg2_e';...
    'EX_mn2_e';'EX_na1_e';'EX_nh4_e';'EX_pi_e';'EX_so4_e';...
    'EX_zn2_e'}; % etc...

f = find(ismember(GFmodel.varNames,strcat('R_',media)));
GFmodel.var_ub(f) = 50;

%% if you want to gapfill for each bbb
% I suggest that you gap-fill per BBB for small scale gap-filling
% repeat for each "problematic BBB"
GFmodeli = GFmodel;

GFmodeli.var_ub(find(ismember(GFmodeli.varNames,strcat('F_DM_',bbbs_m)))) = 50;
GFmodeli.var_ub(find(ismember(GFmodeli.varNames,strcat('NF_DM_',bbbs_m)))) = 50;

GFmodeli.var_lb(find(ismember(GFmodeli.varNames,strcat('NF_',rxnID{i})))) = 10^-6;
GFmodeli.var_lb(find(ismember(GFmodeli.varNames,strcat('F_',rxnID{i})))) = 10^-6;
%% generate the alternative solutions
NumAlt = 10000;
tagMin = 1; % solutions of min size
tagMinP1 = 1;  % solutions of min +1 size
tagSave = 1;
time = 600; % time limit for the solver in seconds
filename = './gf.mat';
indUSE = GFmodeli.indUSE;
mingrRate = 0.01; % force growth rate if you want to gap-fill for each bbb then mingrRate = 0
rxnsToCheck = {};
[ActRxns, GFSumRes, DPsAll] = gapFilling(GFmodeli,indUSE,NumAlt,mingrRate,rxnsToCheck,tagMin,time,filename)
save(strcat('/youpath/,GF_',organism_id,'.mat'),'ActRxns','DPsAll')
%% generate merged model with thermodynamic constarints for evaluation of alternatives
[GFmodel, conflict] = PrepareForGapFillingThermo(sourceModel, {DBmodel},'', 0,flagTFA,{},[],thermo_data);