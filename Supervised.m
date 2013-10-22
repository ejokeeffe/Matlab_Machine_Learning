%> @brief Class that contains static functions that are used for machine
%> learning- all aspects
classdef Supervised < handle
 
    
    properties
    end
    
    methods (Static)
        %> @brief gets the MSE (SSE/size of sample) from a linear regression model
        function [mse] = cost_linreg( X,y,w )
            m=size(X,1);
            mse = (1/m)*(X*w-y)'*(X*w-y);
        end
        %> @brief gets the MSE (SSE/size of sample) from a linear regression model
        %> but including a penalising paramter lambda. This is a redundant
        %>function as it is not actually the mse. Its only useful to compare
        %>cost_linreg with cost_ridgereg. The effect of lambda is accounted
        %>for when calculating w.
        function [mse] = cost_ridgereg( X,y,w,lambda )
            m=size(X,1);
            sse = (X*w-y)'*(X*w-y) + lambda*w'*w;
            mse = sse/m;
        end
        %> @brief Feature Map of the matrix X to a model of order k 
        %> The dimension of phi_X is n*k+1 where nis the dimensions of X
        function [ phi_X ] = powerBasis( X,k )
            phi_X = bsxfun(@power,X,[1:k]);
            %add the offset
            phi_X = [ones(size(phi_X,1),1),phi_X];
        end
        %> @brief Straight linear regression, no need for the function,
        %> easy to do in one line but good for explanatory purposes at least
        
        function w = linreg(X,y)
            w = X\y;
        end %linreg
        %> @brief Gets the dual cost of a kernel function
        function [ mse ] = dualcost( K,y,alpha )
            mse = (K*alpha-y)'*(K*alpha-y)/length(y);
        end
        %> @brief Calculates the gaussian kernel for a given X and sigma
        function K = rbf(X, sigma)
            %K = exp(-1/(2*sigma^2)*(X*X')^2);

            n=size(X,1);
            K=X*X'/sigma^2;
            d=diag(K);
            K=K-ones(n,1)*d'/2;
            K=K-d*ones(1,n)/2;
            K=exp(K);
        end
        %> @brief Calculates the dual weight vector (alpha) 
        %>
        %> @param K kernel matrix (Gram matrix)
        %> @param y dependent variable
        %> @param lambda Regularisation parameter
        function [ alpha ] = kridgereg( K,y,lambda )
            R = lambda*eye(size(K,1))*size(K,1); 
            alpha = (K + R)\y;
        end
        %> @brief calculates the weights, w, given a complexity parameter
        %>
        %> @param X matrix of inputs independent variables. It can also be
        %> the data following a feature mapping from a basis function
        %> @param y output dependent data
        %> @param lambda the regularisation parameter
        %>
        function w = ridgereg(X,y,lambda)

            R = lambda*eye(size(X,2))*size(X,2);
            w = ((X'*X) + R)\X'*y;
        end %ridgereg
    end %static methods
    
end

