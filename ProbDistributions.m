%> @brief used to update probability distributions based on data. Just done
%> for binomial and multinomial data for now
%>
%> @author Eoin O'Keeffe
%> 
%> @version 1.0: Just discrete distributions
%>
%> @todo add continuous distributions
%>
classdef ProbDistributions < handle

    
    properties
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief update an existing binomial distribution (the original 
        %> 'distribution of the distribution' is defined by a and b. These
        %> are the parameters used for the beta distribution which is the 
        %> conjugate of the binomial distribution
        %>
        %> @param a one of the parameters defining the distribution
        %> @param b the other parameter defining the distribution of the
        %> distribution
        %> @param m the number of positive instances in the sample of l
        %> @param l N (the sample size) = m + l. 
        %> @param plot [optional] 1|0 indicates whether we should plot the results
        %> @retval est_mu the updated distribution of the distribution
        %> @retval mle_mu the estimator of mu
        %> @retval a the updated parameter a fo the beta distrubtion
        %> @retval b the updated parameter b of the beta distribution
        function [est_mu,mle_mu,a,b]  = UpdateBinomialDistribution(a,b,m,l,plot)
            %do we want to plot?
            if nargin==4
                plot=0;
            end %if
            %first plot old mu
        mu=[0.1:0.1:0.9];
        est_mu = betapdf(mu,a,b);

        %plot it
        if plot==1
            plot(mu,est_mu,'r-');
        end %if
        %now update

        hold on
        for i=1:m
            a = a+1;
            b=b+l;
            est_mu = betapdf(mu,a,b);
            if plot==1
                plot(mu,est_mu,'b-');
            end %if
            %also get estimator of mu
            mle_mu = a/(a+b);
        end %for i
        end %UpdateBinomialDistribution
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %> @brief updates the multinomial distribution
        %> We use directlet to do this. From directlet you get a
        %> probability of a particular probability distribution. You update
        %> the alpha parameters using the new data and this in turn changes
        %> the probability of the probability distribution. However, we then
        %> need to run MLE on direchlet to determine what the new most likely
        %> probability distribution should be. No mention of this in Bishop -
        %> look elsewhere
        %> @param alpha the parameters of the dirichlet model 
        %> @param m size is the cardinality of the multinomial distribution
        %> It is observations of each value of the distribution
        %> @param N the total number of observations in the sample
        %> @param mu probability distribution over mu
        %> @param plot [optional] should we plot?
        
        function [est_mu,mle_mu,new_alpha] = updateMultinomialDistribution(alpha,m,N,mu,plot)
            %whats the cardinality
            card = length(m);
            %plotting
            if nargin==3
                plot = 0;
            end %if
            
            alpha_0 = sum(alpha);
            new_alpha = alpha(:);
            %add the sample results
            new_alpha = new_alpha + [m(:)];
            
            % Get our mu distribution
            mu = [0:0.1:1]';
            mu = repmat(mu,1,card);
            %plot the old sample
            if plot==1
                %gammas
                gamma_0 = gamma(alpha_0);
                gamma_cards = prod(gamma(alpha));
                
                dir_mu = zeros(size(mu));
                
                for i=1:size(mu,1)
                    cur_mu = mu(i,:);
                    dir_mu(i) = (gamma_0/gamma_cards)*prod(cur_mu(:).^(alpha(i)-1));
                end %for
                surf(mu(:,1),mu(:,2),dir_mu);
            end %if
        end %UpdateMultinomialDistribution
    end%methods
    
end %class

