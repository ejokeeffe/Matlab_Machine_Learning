%> @brief Returns the joint probability distribution over a series of input
%> prob matrices.This is all discrete probability
%>
%> It's largely derived from coursera work and also the brml toolbox which
%> it interacts with
%> @author Eoin O'Keeffe
%>
%> @version
%> Version 1.0: Tried to do as a generic class. Didn't work
%> <br /><i>Version 1.1:</i> Redone so it is specific to the HMM. The basic
%> functioning is to read in lots of 2-D conditional distributions and then
%> return a joint contiditional distribution conditioned all over the same
%> variable. It's not very flexible.
%>
%> @date 15/10/2012
%> <br /> <i>Version 1.1:</i> 22/10/2012
%>
%> @todo make more generic - maybe try to bring in more of the BRML toolbox
classdef ProbAnalysis  < handle
    
    properties
        varnames = [];
        vars = [];
        %stores the cardinals of each variable, cardinals is the number of values
        % that the distribution holds
        card = []; 
        %stores the read in distributions, these could joint, marginal,
        %doesn't matter. Whatever variables they are marginalised over, it
        %is assumed that the variables in the joint distribution are
        %independent of that variables? - is that correct?
        dist = []; %structure of form .vars, .vals
        
        %>stores the conditional distributions read in
        %>Each input is a two dimensional conditional distribution with the
        %> second variables being the conditioned variable p(A|B), the sum
        %>over each column is 1. A is the row and B is the column.
        condDist = [];
    end %properties
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief adds another distribution to the class
        %>
        %> @param inputVals - multi-d matrix of either a distribution or a
        %> joint distribution. Doesn't need to be normalised
        %> @param array of integers/cell array of namess indicating the
        %> variables in the distribution
        function addDist(obj,inputVals,inputVars)
          if isa(inputVars,'cell')
              %it's names
              %Either resolve back or add them to list
              tmpVals = zeros(size(inputVars));
              for i=1:length(inputVars)
                  tmpVals(i) = obj.resolveVariables(inputVars(i));
                  %add the cardinality
                  obj.card(tmpVals(i)) = size(inputVals,i);
              end %for i
              inputVars=  tmpVals;
              clear 'tmpVals';
          end %if
          %so we have in the inputVars resolved back to integers
          %so now add to the joint distribution
          obj.addDistVals(inputVals,inputVars);
          
        end %addDist
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief adds another conditional distribution to the list.
        %> Currently we are assuming two variables where the second one is
        %> the condtioned variable
        %> \f[ p(A|B) \f]
        %>
        %> @param inputVals 
        %> \f[ p(A|B) \f]
        %> where sum over the cols is equal to one (A is on the rows, B is
        %> on the colums)
        %> @param array of integers/cell array of namess indicating the
        %> variables in the distribution
        function addCondDist(obj,inputVals,inputVars)
          if isa(inputVars,'cell')
              %it's names
              %Either resolve back or add them to list
              tmpVals = zeros(size(inputVars));
              for i=1:length(inputVars)
                  tmpVals(i) = obj.resolveVariables(inputVars(i));
                  %add the cardinality
                  obj.card(tmpVals(i)) = size(inputVals,i);
              end %for i
              inputVars=  tmpVals;
              clear 'tmpVals';
          end %if
          %so we have in the inputVars resolved back to integers
          %so now add to the joint distribution
          indx = length(obj.condDist)+1;
        obj.condDist(indx).vals = inputVals;
        obj.condDist(indx).vars = inputVars;
          
        end %addCondDist
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief Checks the varname against the list we have already. If
        %> there, returns the integer, if not it adds it.
        %>
        %> @param inputNames array of cell string of the variable name
        %> @retval inputVars the integer corresponding to the variable
        function inputVars = resolveVariables(obj,inputNames)
            inputVars = zeros(length(inputNames),1);
            for i=1:length(inputNames)
                if ~isempty(obj.vars)
                    indx = cellfun(@(x) strcmp(inputNames{i},x),obj.varnames);
                    if ~isempty(obj.varnames(indx,:))
                        inputVars(i) = obj.vars(indx);
                    else
                        inputVars(i) = obj.addVariable(inputNames(i));
                    end %if
                else
                    inputVars(i) = obj.addVariable(inputNames(i));
                end %if
            end %for
        end %resolveVariables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief adds a variable give a string for the name
        %>
        %> @param inputName string of the variable name
        %> @retval inputVar the integer corresponding to the variable
        function inputVar = addVariable(obj,inputName)
            obj.vars = [obj.vars;length(obj.vars)+1];
            inputVar = length(obj.vars);
            %add the name
            obj.varnames = [obj.varnames;{char(inputName)}];
            %also add a cardinality for it
            obj.card = [obj.card;0];
        end %addVariable
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief Adds a distribution - it can be joint
        %> 
        %> @param inputVals the distribution
        %> @param inputVars array of integer vars
        function addDistVals(obj,inputVals,inputVars)
                indx = length(obj.dist)+1;
                obj.dist(indx).vals = inputVals;
                obj.dist(indx).vars = inputVars;
        end %addDistVals
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief Gets the joint distribution over all the variables (it
        %> assumes that all their values have been read in) as a
        %> multi-dimensional matrix
        function [jd] = getJointDistribution(obj)
            %create the matrix to hold the jd
            jd = ones(obj.card');
            %multiply out all the distributions
            for i=1:length(obj.dist)
                %repmat this distribution over the all the dimensions of
                %the full joint matrix
                %first we need to identify the dimensions over which to
                %repeat. We are assuming the dist is in the correct order
                %of the variables- ie not like [3 1 5] for example
                %for each var, we set the corresponding value to 1- in
                %other words we don;t replicate across this variable as
                %it's distributed
                repmatDims = size(jd);
                for j=1:length(obj.dist(i).vars)
                    repmatDims(obj.dist(i).vars(j)) =  1;
                end %for
                %now repmat the distribution over these dimensions and
                %then flatten out this and jd and multiply by each other
                %and then put jd back into its correct shape
                tmpDist = repmat(obj.dist(i).vals,repmatDims);
                %now flatter
                tmpDist = tmpDist(:);
                jd = reshape(tmpDist.*jd(:),size(ones(obj.card')));
            end %for
        end %function getJointDistribution
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief This runds getJointDistribution above, but then reduces
        %it to an n-dimensional matrix in removing the dimensions for the
        %variables given
        %>
        %> @param dims the Dimensions (as variables) which we want to keep
        %> @retval n-d matrix (where n = no of dims of jd - no of vars in
        %>dims)
        function [jd] = getJointDistAndStripDimensions(obj,dims)
            %first get our joint distribution
            jd = obj.getJointDistribution();
            %Move the ones we want to the front
            remDims = [];
            reshapeDims = [];
            
            for i=1:length(obj.vars)
                if isempty(find(dims==obj.vars(i)))
                    remDims = [remDims,obj.vars(i)];
                else
                    reshapeDims = [reshapeDims,obj.card(i)];
                end %if
            end %for i
            dims = dims(:)';
            %Move the correct dimensions to the front to remove everything
            %else
            jd = permute(jd,[dims remDims]);
            %Now reshape so we only keep the dimensions requested
            jd = reshape(jd,reshapeDims,[]);
            %Now normalize
            % we need the brml toolbox
            jd = brml.condp(jd);
        end %getJointDistAndStripDimensions
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief This function gets the marginalised distribution over a 
        %> variable. any dists involving this variable should return the same
        %> marginal distribution!!!!!!!!!!
        %>
        %> @param inputVar either an integer of a name
        %> @retval d the distribution
        function [d] = getMarginal(obj,inputVar)
            if isa(inputVar,'cell')
                inputVar = obj.resolveVariable(inputVar);
            end %if
            %loop through all the dists until we reach the one we want,
            for i=1:length(obj.dist)
                if ~isempty(find(obj.dist(i).vars == inputVar))
                    idx = find(obj.dist(i).vars == inputVar);
                    %need to permute to get this to the front
                    otherIdxs = obj.dist(i).vars;
                    otherIdxs(otherIdxs==idx) = [];
                    d = obj.dist.val;
                    d = permute(d,[idx otherIdxs']);
                    %this is a pain but sum in reverse order
                    for j=length(size(d)):2
                        d = sum(d,j);
                    end %for  j
                    %that should do it.....
                    break;
                end%if
            end 
        end %getMarginal
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief gets the conditional distribution over a selected
        %> variable. Not flexible at the moment - assumes all the
        %> condDist's are conditioned on the same variable
        %> say we pass in (obj,[A B C],D)
        %> we would yield p(A,B,C|D)
        %> Matrix of size rows = card(A)*card(B)*card(C), cols = card(D)
        %>
        %> @param obj instance of Prob_eok class
        %> @param jdVars the input variable overwhich to get the
        %> distribtuion
        %> @param conditionedVar the variable over which to condition
        %> @retval jd 2-d matrix of joint distribution
        function jd = getJointConditionalDistribution(obj,jdVars,conditionedVar)
            %instantiate the matrix
            rowSize = 1;
            cards = [];
            %convert the list of vars into a sequential array with values
            %corresponding to the indexes of the condDist that these vars
            %correspond to 
            condDistIndexes = [];
            for i=1:length(jdVars)
                rowSize = rowSize*obj.card(jdVars(i));
                cards = [cards;obj.card(jdVars(i))];
                for j=1:length(obj.condDist)
                    if obj.condDist(j).vars(1) == jdVars(i)
                        condDistIndexes = [condDistIndexes j];
                    end %if
                end %for j
            end %for i
            
            jd = ones(rowSize,obj.card(conditionedVar));
            %recursively factorise across each variable
            jd = obj.recursiveFactorise(jd,condDistIndexes);
            %Now normalize on the columns
            jd = brml.condp(jd);
        end %function getJointConditionalDistribution
        
    end %public methods
    
    methods (Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief Recursive calculate a joint distribution.
        %>
        %> @param obj instance of Prob_eok
        %> @param oldPhi under construction factorised matrix
        %> @param remainingCondDists an array of indexes corresponding to
        %> the conditional distributions
        %> @retval phi returned factorisation
        function phi = recursiveFactorise(obj,oldPhi,remainingCondDists)
            totalLength = size(oldPhi,1);
            if length(remainingCondDists) > 1
                %recursively get the next phi
                oldPhi = obj.recursiveFactorise(oldPhi,remainingCondDists(2:end));
            end %if
            %take next var in line
            selIndx = remainingCondDists(1);
            % length of repition
            repLength = 1;
            for j=2:length(remainingCondDists)
                %multiply across by the cardinality
                repLength = repLength*obj.card(obj.condDist(remainingCondDists(j)).vars(1));
            end %for j
            newPhi = [];
            %rep each col over the remaining conditional distributions
            for i=1:obj.card(obj.condDist(selIndx).vars(1))
                newPhi = [newPhi;repmat(obj.condDist(selIndx).vals(i,:),repLength,1)];
            end %for i
            %finally repeat the whole matrix over the total length
            repLength = totalLength/size(newPhi,1);
            newPhi = repmat(newPhi,repLength,1);
            %now multiply by oldPhi to factorise
            phi = newPhi.*oldPhi;
        end %recursiveFactorise
    end %private methods
    
    methods (Static)
        %stores generally useful static methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief generates a transition matrix for a HMM (may have other
        %> applications, not sure)
        %>
        %> @param vals a structure of a series of states, each with an
        %>integer value, this value should be in the states array. It
        %> doesn't matter what the names of the each of the indexes of the
        %> structure are as we loop through all of them
        %> @param states a vector of states in the order that the caller wishes
        %> them to appear in the resulting matrix
        %> @retval phgh Probability of h given h where h is the hidden or
        %> state variable 
        %> \f[ p(h_t|h_{t-1}) \f]
        %> where  h_t is the row and h_t-1 is the column
        function phgh = generateTransitionMatrix(vals,states)
            noStates = length(states);
            %make sure states is the right way (ie a column vector)
            states = states(:);
            phgh = zeros(noStates,noStates);
            train = fieldnames(vals);
            %First get a factor transition matrix based on the training set above
            for i=1:length(train)
                for j=2:length(vals.(train{i}))
                    %get previous and current state index!!! (not the state as
                    %the state may not be in the index location)
                    tminus1 = find(states==vals.(train{i})(j-1));
                    t = find(states==vals.(train{i})(j));
                    %add to our matrix
                        phgh(t,tminus1)=...
                            phgh(t,tminus1)+1;
                end %for j
            end %for i
            %Now convert to a conditional matrix (cols sum to 1)
            phgh = brml.condp(phgh);
        end %generateTransitionMatrix
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief Similar to generateTransitionMatrix except this function
        %> updates and existing transition matrix based on data
        %> /f[ p(mu|D) /f]
        %>
        %> @param D New state data that is to be use dto update the
        %> transition matrix
        %> @param mu The existing transition matrix
        %> @param alpha Existing transition matrix parameters
        %> @retval mu The new transition matrix
        function [mu] = updateTransitionMatrix(D,mu,alpha)
            
        end %updateTransitionMatrix
        
    end %static methods
        
end %classdef