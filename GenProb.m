%> @brief This is a general probability class. It was created to
%> incorporate the code from Coursera PGM
%>
%> @author Eoin O'Keeffe
%>
%> @version 1.0
%>
%> @date 01/11/2012
%>
%>
classdef GenProb < handle
 
    
    properties
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % @brief FactorProduct Computes the product of two factors.
        %   C = FactorProduct(A,B) computes the product between two factors, A and B,
        %   where each factor is defined over a set of variables with given dimension.
        %   The factor data structure has the following fields:
        %       .var    Vector of variables in the factor, e.g. [1 2 3]
        %       .card   Vector of cardinalities corresponding to .var, e.g. [2 2 2]
        %       .val    Value table of size prod(.card)
        %
        %   See also FactorMarginalization.m, IndexToAssignment.m, and
        %   AssignmentToIndex.m
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function C = FactorProduct(A, B)

        % Check for empty factors
        if (isempty(A.var)), C = B; return; end;
        if (isempty(B.var)), C = A; return; end;

        % Check that variables in both A and B have the same cardinality
        [dummy iA iB] = intersect(A.var, B.var);
        if ~isempty(dummy)
            % A and B have at least 1 variable in common
            assert(all(A.card(iA) == B.card(iB)), 'Dimensionality mismatch in factors');
        end

        % Set the variables of C
        C.var = union(A.var, B.var);

        % Construct the mapping between variables in A and B and variables in C.
        % In the code below, we have that
        %
        %   mapA(i) = j, if and only if, A.var(i) == C.var(j)
        % 
        % and similarly 
        %
        %   mapB(i) = j, if and only if, B.var(i) == C.var(j)
        %
        % For example, if A.var = [3 1 4], B.var = [4 5], and C.var = [1 3 4 5],
        % then, mapA = [2 1 3] and mapB = [3 4]; mapA(1) = 2 because A.var(1) = 3
        % and C.var(2) = 3, so A.var(1) == C.var(2).

        [dummy, mapA] = ismember(A.var, C.var);
        [dummy, mapB] = ismember(B.var, C.var);

        % Set the cardinality of variables in C
        C.card = zeros(1, length(C.var));
        C.card(mapA) = A.card;
        C.card(mapB) = B.card;

        % Initialize the factor values of C:
        %   prod(C.card) is the number of entries in C
        C.val = zeros(1, prod(C.card));

        % Compute some helper indices
        % These will be very useful for calculating C.val
        % so make sure you understand what these lines are doing.
        assignments = IndexToAssignment(1:prod(C.card), C.card);
        indxA = AssignmentToIndex(assignments(:, mapA), A.card);
        indxB = AssignmentToIndex(assignments(:, mapB), B.card);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % YOUR CODE HERE:
        % Correctly populate the factor values of C
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        C.val = A.val(indxA).*B.val(indxB);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end %FactorProduct
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief FactorMarginalization Sums given variables out of a factor.
        %>   B = FactorMarginalization(A,V) computes the factor with the variables
        %>   in V summed out. The factor data structure has the following fields:
        %>       .var    Vector of variables in the factor, e.g. [1 2 3]
        %>       .card   Vector of cardinalities corresponding to .var, e.g. [2 2 2]
        %>       .val    Value table of size prod(.card)
        %>
        %>   The resultant factor should have at least one variable remaining or this
        %>   function will throw an error.
        % 
        %>   See also FactorProduct.m, IndexToAssignment.m, and AssignmentToIndex.m
        %
        %> @param A
        %> @param V
        %> @retval B 
        function B = FactorMarginalization(A, V)

        % Check for empty factor or variable list
        if (isempty(A.var) || isempty(V)), B = A; return; end;

        % Construct the output factor over A.var \ V (the variables in A.var that are not in V)
        % and mapping between variables in A and B
        [B.var, mapB] = setdiff(A.var, V);

        % Check for empty resultant factor
        if isempty(B.var)
          error('Error: Resultant factor has empty scope');
        end;

        % Initialize B.card and B.val
        B.card = A.card(mapB);
        B.val = zeros(1, prod(B.card));

        % Compute some helper indices
        % These will be very useful for calculating B.val
        % so make sure you understand what these lines are doing
        %assignments = GenProb.IndexToAssignment(1:length(A.val), A.card);
        AIndexes = find(A.val>0);
        assignments = GenProb.IndexToAssignment(AIndexes, A.card);
        indxB = GenProb.AssignmentToIndex(assignments(:, mapB), B.card);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % YOUR CODE HERE
        % Correctly populate the factor values of B
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:size(indxB,1)
                B.val(indxB(i)) = B.val(indxB(i)) + A.val(AIndexes(i));
            end %for i
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end %FactorMarginalisation
        
        %> ***************************************************************
        %>
        %> Computes the two dimension frequency matrix for factor A over
        %> the variables in V
        %>
        %> @param A Factor
        %> @param vector of vars (eg [1 3] etc) over which to generate the
        %> frequency matrix. It can only contain two values
        %>
        %> @retval B the returned frequency matrix
        function B = FrequencyMatrix(A,V)
            Marg = A;
            if length(A.var)~=length(V)
                notV = setdiff(A.var,V);
                Marg = GenProb.FactorMarginalization(A,notV);
            end
                
            indexes = find(Marg.val>0);
            B = zeros(A.card(V));
            
            ass = GenProb.IndexToAssignment(indexes, A.card(V));
            for i=1:length(indexes)
                B(ass(i,1),ass(i,2)) = Marg.val(indexes(i));
            end %for i
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> ObserveEvidence Modify a vector of factors given some evidence.
        %>   F = ObserveEvidence(F, E) sets all entries in the vector of factors, F,
        %>   that are not consistent with the evidence, E, to zero. 
        %>
        %> @param F is a vector of
        %>   factors, each a data structure with the following fields:
        %>     .var    Vector of variables in the factor, e.g. [1 2 3]
        %>     .card   Vector of cardinalities corresponding to .var, e.g. [2 2 2]
        %>     .val    Value table of size prod(.card)
        %>   @param E is an N-by-2 matrix, where each row consists of a variable/value pair. 
        %>     Variables are in the first column and values are in the second column.
        %> @retval F unnormalized
        function F = ObserveEvidence(F, E)

                % Iterate through all evidence

                for i = 1:size(E, 1),
                    v = E(i, 1); % variable
                    x = E(i, 2); % value

                    % Check validity of evidence
                    if (x == 0),
                        warning(['Evidence not set for variable ', int2str(v)]);
                        continue;
                    end;

                    for j = 1:length(F),
                          % Does factor contain variable?
                        indx = find(F(j).var == v);

                        if (~isempty(indx)),

                               % Check validity of evidence
                            if (x > F(j).card(indx) || x < 0 ),
                                error(['Invalid evidence, X_', int2str(v), ' = ', int2str(x)]);
                            end;

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % YOUR CODE HERE
                            % Adjust the factor F(j) to account for observed evidence
                            % Hint: You might find it helpful to use IndexToAssignment
                            %       and SetValueOfAssignment
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            for k=1:numel(F(j).val)
                                %get the assignment of this index
                                curAssignment = IndexToAssignment(k,F(j).card);
                                %if the value  of the variable x is not v then set to 0
                                if curAssignment(indx)~= x
                                    F(j) = SetValueOfAssignment(F(j),curAssignment,0);
                                end %if
                            end %for k
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                % Check validity of evidence / resulting factor
                            if (all(F(j).val == 0)),
                                warning(['Factor ', int2str(j), ' makes variable assignment impossible']);
                            end;

                        end;
                    end;
                end;

                end %ObserveEvidence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief  ComputeMarginal Computes the marginal over a set of given variables
        %>   M = ComputeMarginal(V, F, E) computes the marginal over variables V
        %>   in the distribution induced by the set of factors F, given evidence E
        %>
        %>  @param V is a vector containing the variables in the marginal e.g. [1 2 3] for
        %>     X_1, X_2 and X_3.
        %> @param   F is a vector of factors (struct array) containing the factors 
        %>     defining the distribution
        %> @param  E is an N-by-2 matrix, each row being a variable/value pair. 
        %>     Variables are in the first column and values are in the second column.
        %>     If there is no evidence, pass in the empty matrix [] for E.
        %>   @retval M is a factor containing the marginal over variables V
        function M = ComputeMarginal(V, F, E)

        % Check for empty factor list
        if (numel(F) == 0)
              warning('Warning: empty factor list');
              M = struct('var', [], 'card', [], 'val', []);      
              return;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % YOUR CODE HERE:
        % M should be a factor
        % Remember to renormalize the entries of M!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First use the evidence to update the probabilities
        FgE = ObserveEvidence(F,E);
        % Now get the joint distribution
        Joint = ComputeJointDistribution(FgE);
        %Now marginalize over the given variables
        %first get the variables we marginalize out
        removeV = Joint.var(~ismember(Joint.var,V));
        % Now marginalize over the selected vectors
        M = FactorMarginalization(Joint,removeV);
        % Now normalize
        M.val = M.val/sum(M.val);
        %M = struct('var', [], 'card', [], 'val', []); % Returns empty factor. Change this.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end%ComputeMarginal

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %> @brief ComputeJointDistribution Computes the joint distribution defined by a set
        %> of given factors
        %>
        %>   Joint is a factor that encapsulates the joint distribution given by F
        %> @param   F is a vector of factors (struct array) containing the factors 
        %>     defining the distribution
        %> @retval   Joint = ComputeJointDistribution(F) computes the joint distribution
        %>   defined by a set of given factors

        function Joint = ComputeJointDistribution(F)

          % Check for empty factor list
          if (numel(F) == 0)
              warning('Error: empty factor list');
              Joint = struct('var', [], 'card', [], 'val', []);      
              return;
          end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % YOUR CODE HERE:
        % Compute the joint distribution defined by F
        % You may assume that you are given legal CPDs so no input checking is required.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First find out how many variables there are, sort these and set the cardinals also
        var = [];
        for i=1:numel(F)
            var = [var F.var];
        end %for i
        var = sort(unique(var));
        card = zeros(1,numel(var));
        for i=1:numel(F)
            for j=1:numel(F(i).var)
                card(F(i).var(j)) = F(i).card(j);
            end %for j
        end %for i
        % Now get val
        val = zeros(1, prod(card));

        %now loop through val and get factorisations
        for i=1:numel(val)
            val(i) = 1;
            %get the joint assignment for this value index of the joint
            %distribution
            jointAssignment = GenProb.IndexToAssignment(i,card);
            %loop through each of the factors multiplying them out in turn
            for j=1:numel(F)
                %create vector of assignments for this factor 
                assign = zeros(1,numel(F(j).var));
                %now loop through assign and match to the variable so we can get
                %the assignment for this factor
                for k=1:numel(assign)
                    assign(k) = jointAssignment(var==F(j).var(k));
                end %k
                indexOfValue = GenProb.AssignmentToIndex(assign,F(j).card);
                val(i) = val(i)*F(j).val(indexOfValue);
            end %for j
        end %for i


        Joint = struct('var', var, 'card', card, 'val', val); % Returns empty factor. Change this.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end %computejointdistribtuion
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief AssignmentToIndex Convert assignment to index.
        %
        %>   @retval I = AssignmentToIndex(A, D) converts an assignment, A, over variables
        %>   with cardinality D to an index into the .val vector for a factor. 
        %>   If A is a matrix then the function converts each row of A to an index.
        %>
        %>   See also IndexToAssignment.m 

        function I = AssignmentToIndex(A, D)

        D = D(:)'; % ensure that D is a row vector
        if (any(size(A) == 1)),
            I = cumprod([1, D(1:end - 1)]) * (A(:) - 1) + 1;
        else
            I = sum(repmat(cumprod([1, D(1:end - 1)]), size(A, 1), 1) .* (A - 1), 2) + 1;
        end;

        end%AssignmentToIndex

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief IndexToAssignment Convert index to variable assignment.
        %>
        %> @param I index
        %> @param D
        %>   @retval A = IndexToAssignment(I, D) converts an index, I, into the .val vector
        %>   into an assignment over variables with cardinality D. If I is a vector, 
        %>   then the function produces a matrix of assignments, one assignment 
        %>   per row.
        %>
        %>   See also AssignmentToIndex.m and FactorTutorial.m

        function A = IndexToAssignment(I, D)

        D = D(:)'; % ensure that D is a row vector
        A = mod(floor(repmat(I(:) - 1, 1, length(D)) ./ repmat(cumprod([1, D(1:end - 1)]), length(I), 1)), ...
                repmat(D, length(I), 1)) + 1;

        end%IndexToAssignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief GetValueOfAssignment Gets the value of a variable assignment in a factor.
%
%   v = GetValueOfAssignment(F, A) returns the value of a variable assignment,
%   A, in factor F. The order of the variables in A are assumed to be the
%   same as the order in F.var.
%
%   v = GetValueOfAssignment(F, A, VO) gets the value of a variable assignment,
%   A, in factor F. The order of the variables in A are given by the vector VO.
%
%   See also SetValueOfAssignment.m and FactorTutorial.m

function v = GetValueOfAssignment(F, A, VO)

if (nargin == 2),
    indx = AssignmentToIndex(A, F.card);
else
    map = zeros(length(F.var), 1);
    for i = 1:length(F.var),
        map(i) = find(VO == F.var(i));
    end;
    indx = AssignmentToIndex(A(map), F.card);
end;

v = F.val(indx);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief SetValueOfAssignment Sets the value of a variable assignment in a factor.
%
%   F = SetValueOfAssignment(F, A, v) sets the value of a variable assignment,
%   A, in factor F to v. The order of the variables in A are assumed to be the
%   same as the order in F.var.
%
%   F = SetValueOfAssignment(F, A, v, VO) sets the value of a variable
%   assignment, A, in factor F to v. The order of the variables in A are given
%   by the vector VO.
%
%   Note that SetValueOfAssignment *does not modify* the factor F that is 
%   passed into the function, but instead returns a modified factor with the 
%   new value(s) for the specified assignment(s). This is why we have to  
%   reassign F to the result of SetValueOfAssignment in the code snippets 
%   shown above.
%
%   See also GetValueOfAssignment.m and FactorTutorial.m

function F = SetValueOfAssignment(F, A, v, VO)

if (nargin == 3),
    indx = GenProb.AssignmentToIndex(A, F.card);
else
    map = zeros(length(F.var), 1);
    for i = 1:length(F.var),
        map(i) = find(VO == F.var(i));
    end;
    indx = GenProb.AssignmentToIndex(A(map), F.card);
end;

F.val(indx) = v;

end


%> @brief Computes the conditional pdf for this factor. Knowing that the
%> factor is effectively the unnormalised pdf anyway conditioned on the 
%> first variable
%>
%> @param F a factor input - just one
%> @retval P - basically F normalised
%> @retval P_matrix  - if the number of variables is two then this is a two
%> dimensional matrix with the conditioning variable columnwise (ie each
%> column sums to 1)
function [P,P_matrix] = conditionalProbability(F)
    numberOfCondInstances = prod(F.card(2:end));
    P = F;
    P_matrix =[];
    for i=1:numberOfCondInstances
        % get the indexes
        ass = GenProb.IndexToAssignment(i,F.card(2:end));
        indxs = GenProb.AssignmentToIndex([[1:F.card(1)]' repmat(ass,length([1:F.card(1)]),1)],F.card);
        P.val(indxs) = F.val(indxs)/sum(F.val(indxs));
    end %for i
    
    %Now if the number of variables is two pop out as 2-d matrix
    if length(F.var)==2
        P_matrix = zeros(P.card);
        P_matrix = reshape(P.val,P.card);
    end %if
end %conditionalProbability

%-------------------------------------------------------
%> @brief Normalize a factor so it becomes a conditional probability
%>conditioned on the 2nd to last of the vars
%>
%> @param F The input factor
%> 
function [A] =  normaliseFactor(F)
        A = F;
        %Loop through the cardinality of the first var, normalising across
        %the other vars
        for i=1:prod(F.card(2:end))
            ass = GenProb.IndexToAssignment(i,F.card(2:end));
            indxs = GenProb.AssignmentToIndex([[1:F.card(1)]' repmat(ass,F.card(1),1)],F.card);
            A.val(indxs) = A.val(indxs)/sum(A.val(indxs));
        end %for i
    end %normaliseFactor
    
    function [ass] = buildAssignment(cardArray,remCard)
        if length(remCard)>1
            cardArray = GenProb.buildAssignment(cardArray,remCard(2:end));
        end
        
        %repeat the passed matrix
        if ~isempty(cardArray)
            %Create the additional column
            addedAss = repmat([1:remCard],size(cardArray,1),1);
            addedAss = addedAss(:);
            cardArray = repmat(cardArray,remCard,1);
        else
            addedAss = [1:remCard]';
        end %if
        ass = [addedAss cardArray];
    end %buildAssignment
    end %static methods
    
end

