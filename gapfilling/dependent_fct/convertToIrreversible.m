function [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model)
%convertToIrreversible Convert model to irreversible format
%
% [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model)
%
%INPUT
% model         COBRA model structure
%
%OUTPUTS
% modelIrrev    Model in irreversible format
% matchRev      Matching of forward and backward reactions of a reversible
%               reaction
% rev2irrev     Matching from reversible to irreversible reactions
% irrev2rev     Matching from irreversible to reversible reactions
%
% Uses the reversible list to construct a new model with reversible
% reactions separated into forward and backward reactions.  Separated
% reactions are appended with '_f' and '_b' and the reversible list tracks
% these changes with a '1' corresponding to separated forward reactions.
% Reactions entirely in the negative direction will be reversed and
% appended with '_r'.
%
% written by Gregory Hannum 7/9/05
%
% Modified by Markus Herrgard 7/25/05
% Modified by Jan Schellenberger 9/9/09 for speed.

%declare variables
modelIrrev.S = spalloc(size(model.S,1),0,2*nnz(model.S));
modelIrrev.rxns = [];
modelIrrev.rev = zeros(2*length(model.rxns),1);
modelIrrev.lb = zeros(2*length(model.rxns),1);
modelIrrev.ub = zeros(2*length(model.rxns),1);
modelIrrev.c = zeros(2*length(model.rxns),1);
matchRev = zeros(2*length(model.rxns),1);

nRxns = size(model.S,2);
irrev2rev = zeros(2*length(model.rxns),1);

%loop through each column/rxn in the S matrix building the irreversible
%model
cnt = 0;
for i = 1:nRxns
    cnt = cnt + 1;
    
    %expand the new model (same for both irrev & rev rxns  
    modelIrrev.rev(cnt) = model.rev(i);
    irrev2rev(cnt) = i;

    % Reaction entirely in the negative direction
    if (model.ub(i) <= 0 && model.lb(i) < 0)
        modelIrrev.rxns{cnt} = ['F_' model.rxns{i}];
        modelIrrev.S(:,cnt) = model.S(:,i);
        modelIrrev.c(cnt) = -model.c(i);
        modelIrrev.ub(cnt) = 0;
        modelIrrev.lb(cnt) = 0;
        
        cnt = cnt + 1;
                
        % Retain original bounds but reversed
        modelIrrev.ub(cnt) = -model.lb(i);
        modelIrrev.lb(cnt) = -model.ub(i);
        % Reverse sign
        modelIrrev.S(:,cnt) = -model.S(:,i);
        modelIrrev.c(cnt) = -model.c(i);
        modelIrrev.rxns{cnt} = ['R_' model.rxns{i}];
        model.rev(i) = false;
        modelIrrev.rev(cnt) = false;
        model.rev(i) = false;
                
        if (model.lb(i) >= 0)
            modelIrrev.ub(cnt) = 0;
            modelIrrev.lb(cnt) = 0;            
        end
  
    else
        % Keep positive upper bound
        modelIrrev.ub(cnt) = model.ub(i);
        %if the lb is less than zero, set the forward rxn lb to zero 
        if model.lb(i) < 0
            modelIrrev.lb(cnt) = 0;
        else
            modelIrrev.lb(cnt) = model.lb(i);
        end
        modelIrrev.S(:,cnt) = model.S(:,i);
        modelIrrev.c(cnt) = model.c(i);
        modelIrrev.rxns{cnt} = ['F_' model.rxns{i}];

    end
    
    %if the reaction is reversible, add a new rxn to the irrev model and
    %update the names of the reactions with '_f' and '_b'
    if model.rev(i) == true
        cnt = cnt + 1;
        matchRev(cnt) = cnt - 1;
        matchRev(cnt-1) = cnt;
        modelIrrev.rxns{cnt-1} = ['F_' model.rxns{i}];
        
        if (model.ub(i) < 0)
            modelIrrev.ub(cnt-1) = 0;
        else
            modelIrrev.ub(cnt-1) = model.ub(i);
        end
        
        modelIrrev.S(:,cnt) = -model.S(:,i);
        modelIrrev.rxns{cnt} = ['R_' model.rxns{i}];
        modelIrrev.rev(cnt) = true;
        modelIrrev.lb(cnt) = 0;
        
        if (model.lb(i) < 0)
            modelIrrev.ub(cnt) = -model.lb(i);
        elseif (model.lb(i) >= 0)
            modelIrrev.ub(cnt) = 0;
        else
            modelIrrev.ub(cnt) = model.ub(i);
        end
        modelIrrev.c(cnt) = 0;
        rev2irrev{i} = [cnt-1 cnt];
        irrev2rev(cnt) = i;
    else
        matchRev(cnt) = 0;
        rev2irrev{i} = cnt;
    end
end

rev2irrev = columnVector(rev2irrev);
irrev2rev = irrev2rev(1:cnt);
irrev2rev = columnVector(irrev2rev);

% Build final structure
modelIrrev.S = modelIrrev.S(:,1:cnt);
modelIrrev.ub = columnVector(modelIrrev.ub(1:cnt));
modelIrrev.lb = columnVector(modelIrrev.lb(1:cnt));
modelIrrev.c = columnVector(modelIrrev.c(1:cnt));
modelIrrev.rev = modelIrrev.rev(1:cnt);
modelIrrev.rev = columnVector(modelIrrev.rev == 1);
modelIrrev.rxns = columnVector(modelIrrev.rxns); 
modelIrrev.mets = model.mets;
matchRev = columnVector(matchRev(1:cnt));
modelIrrev.match = matchRev;
if (isfield(model,'b'))
    modelIrrev.b = model.b;
end
if isfield(model,'description')
    modelIrrev.description = [model.description ' irreversible'];
end
% if isfield(model,'subSystems')
%     modelIrrev.subSystems = model.subSystems(irrev2rev);
% end
% if isfield(model,'genes')
%     modelIrrev.genes = model.genes;
%     genemtxtranspose = model.rxnGeneMat';
%     modelIrrev.rxnGeneMat = genemtxtranspose(:,irrev2rev)';
%     modelIrrev.rules = model.rules(irrev2rev);
% end
modelIrrev.reversibleModel = false;

