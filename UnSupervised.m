%> @brief Class that encapsulates various clustering functions.
%>
%> @author Eoin O'Keeffe 
%>
%> @date 16/10/2013
%>
%> @detail Forked from Adaptive Modelling of complex data lecutre and 
%> tutorial series @ucl provided by 
%> the Computer Science department. It als includes some test function to
%> show how it works. 
classdef UnSupervised < handle
 
    
    properties
    end
    
    methods (Static)
        %-------------------------------------------------
        %> @brief Initialise the clusters that each datapoint
        %> is associated with
        %>
        %> @param X
        %> @param K
        %>
        %> @retval cc the cluster association vector
        %> @retval mu the average value of each cluster
        %------------------------------------------------
        function [cc, mu] = kmeans_initialise(X,K)
            mm = size(X,1);	
            % number of data items 
            dd = size(X,2);	
            % number of dimensions of data 
            cc = 1+rem(randperm(mm),K); 
            % assign randomly to clusters 
            mu = zeros(K,dd);	
            % initialize means of clusters 
            for kk = 1:K
                mu(kk,:) = mean(X(cc==kk,:),1); 
            end
        end %kmeans_initialise
        %> @brief Uses kmeans clustering to put data into clusters
        %> @param X [m x n] m is hte nmumber of observations and n is the number of
        %> dimensions [optional]
        %> @param K Scalar - Number of clusters [optional]
        %> @param showAsPlot [1/0] shows the evolution of the kmeans
        %> clustering
        %> 
        %> @retval mu [K x n] nean value of variable within cluster
        %> @retval c [m x 1] categorical indicating association of observation with
        %> cluster c
        %> @retval J Scalar Objective function value after convergence
        %> @retval X as above
        function [mu,c,J,X] = kmeans(X,K,showAsPlot)

        if nargin==0
          clf
          say Create your own dataset.
          X = getpoints;
          K = input('Enter number of clusters: ');
          showAsPlot=1;
        end
        if nargin ==1
            showAsPlot=1;
        end %if 

        numrepeats = 10; % repeat K-means 10 times.
        m = size(X,1);
        n = size(X,2);
        %if n ~= 2, error('demo only works in 2D'); end

        % figure out range of data.
        xmin = min(X,[],1);
        xmax = max(X,[],1);
        xwidth = xmax-xmin;

        % keep result of K-means for each repeat.
        allJ = cell(1,numrepeats); 
        allc = cell(1,numrepeats);
        allmu = cell(1,numrepeats);

        for repeat = 1:numrepeats

          col = colormap;
          col = col(ceil((.5:K)/K*size(col,1)),:);

          % initialize cluster assignments c and prototypes mu
          [c mu] = UnSupervised.kmeans_initialise(X,K);

          % compute initial objective function
          J = UnSupervised.kmeans_objective_ls(X,c,mu); 

          % initial plot
          if showAsPlot==1
              subplot(121);
              plotkmeans(X,c,mu);
              title(sprintf('K-means (repeat %d/%d)',repeat,numrepeats));
              subplot(122);
              plot(0,J,'x-','linewidth',3);
              title('Total distance to prototypes');
              pause
          end %if

          % iterate K-means
          for iter = 1:99
            % Run K-means for 1 iteration.  Note if a cluster gets emptied-out we 
            % reassign it to be at some randomly chosen data vector.
            [c mu] = UnSupervised.kmeans_update(X,c,mu);
            J(iter+1) = UnSupervised.kmeans_objective_ls(X,c,mu);

            % plot 
            if showAsPlot==1
                subplot(121);
                cla;
                plotkmeans(X,c,mu);
                title(sprintf('K-means (repeat %d/%d)',repeat,numrepeats));
                subplot(122);
                plot(0:iter,J,'x-','linewidth',3);
                title('Total distance to prototypes');
            end %if

            % test for convergence
            if iter>1 && all(c==oldc),
              break;
            end
            oldc = c;
            if showAsPlot==1
                        pause
            end %if 
          end
          if showAsPlot==1
              subplot(121);
              title(sprintf('K-means converged (repeat %d/%d)',repeat,numrepeats));
              pause
          end %if


          allJ{repeat} = J;
          allc{repeat} = c;
          allmu{repeat} = mu;
        end

        % best run of K-means
        J = zeros(1,numrepeats);
        if showAsPlot==1
            subplot(122);
            cla
            hold on
        end %if
        for repeat = 1:numrepeats
            if showAsPlot ==1
            plot(0:length(allJ{repeat})-1,allJ{repeat},'-');
            end %if
          J(repeat) = allJ{repeat}(end);
        end
        if showAsPlot==1
            hold off
            title('K-means total distances');
        end %if
        [bestJ repeat] = min(J);
        c = allc{repeat};
        mu = allmu{repeat};
        if showAsPlot==1
            subplot(121);
            cla;
            hold on;
            for i = 1:m
              plot([X(i,1) mu(c(i),1)],[X(i,2) mu(c(i),2)],'-','linewidth',1,...
                  'color',[.7 .7 .7]);
              plot(X(i,1),X(i,2),'.','markersize',10,'color',col(c(i),:));
            end
            for k=1:K
              plot(mu(k,1),mu(k,2),'+','markersize',20,'linewidth',3,'color',col(k,:));
            end
            hold off;
            title('Best K-means run');
        end %if
    end %kmeans
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %> @brief 
    %>
    %> @param X
    %>
    %> @param c [M x 1] Categorical associating an observation with a
    %> cluster k
    %>
    %> @param mu [K x N] mean of variables within clusters
    %>
    %> @retval cc
    %> @retval mu update mu 
    function [cc,mu] = kmeans_update(X,c,mu)

        K = size(mu,1);
        n = size(X,2);
        m = size(X,1);
        cc = zeros(size(c));
        for ii = 1:m
            % iterate over data items 
            dist = zeros(K,1); 
            for kk = 1:K
                % squared distance
                dist(kk) = sum((X(ii,:)-mu(kk,:)).^2); 
            end
            % find minimum distance cluster and assign item to it 
            [~,kk] = min(dist); 
            cc(ii) = kk;
        end
        for kk = 1:K % iterate over clusters to update means
            mu(kk,:) = mean(X(cc==kk,:),1); 
        end
    end %kmeans_update
    
        %>@brief Calculates the objective function of the current associations
        %> It uses the square of mean deviations within each cluster, the
        %> distance metric is euclidean distance.
        %>
        %> @param X [m x n] matrix of observations
        %> @param c [m x 1] vector number 1 to k assocating an observation with a
        %> cluster
        %> @param mu [k x n] 
        %> 
        %> @retval J objective function value
        function J = kmeans_objective_ls(X,c,mu)

            J=0;
            m=size(X,1);
            n=size(X,2);

            %loop through the observations
            for ii =1:m
                %get squared deviation of current observation
                sqDev = sum((X(m,:)'-mu(c(ii),:)').^2);
                J = J+sqDev;
            end %for ii
        end %kmeans_objective_ls
        
        %-------------------------------------------------
        %> @brief Plots the cluster for X (2 d only) as associated with
        %> clusters indicated in the cluster association vector c
        %>
        %> @param X [m x 2] Data matrix
        %> @param c [m x 1] Cluster association vector
        %>
        %> @retval h handle to the data in the plots 
        function h = plotClusters(X,c)

            % plot data, coloring according to cluster labels.


            if max(c) > 1
              K = max(c);
              col = colormap;
              col = col(ceil((.5:K)/K*size(col,1)),:);
              h = zeros(1,K);
              for k=1:K
                i = find(c==k);
                h(k) = plot(X(i,1),X(i,2),'.','markersize',20,'color',col(k,:));
                hold on
              end
              hold off
            else
              r = c;
              K = size(r,2);
              m = size(X,1);
              col = colormap;
              col = col(ceil((.5:K)/K*size(col,1)),:);
              h = zeros(1,m);
              for i = 1:m
                h(i) = plot(X(i,1),X(i,2),'.','markersize',10,'color',min(r(i,:)*col,1));
                hold on
              end
            end
        end%plotClusters


        %-------------------------------------------------------------
        %> @brief Perform Spectral clustering on a dataset. It uses the RBF
        %> kernel as the similarity measure
        %>
        %> @param X [m x n] dataset
        %> @param K the number of clusters
        %> @param sigma
        %>
        %> @retval cc  [m x 1] Clustering association vector
        %> @retval mu [K x n] average values of each cluster in each
        %> dimension
        %> @retval Y Normalised eigenvector matrix
        %> @retval W 
        function [cc,mu,Y,W] = speclust(X,K,sigma)
          
            %
            % Compute the W and D matrices 
            %

            %For simplicity we will use the RBF kernel as a similarity measure:
            [mm dd] = size(X); 
            % compute distances 
            nx = sum(X.^2,2); 
            % This is the squared distance between each entry
            D = nx(:,ones(1,mm)) + nx(:,ones(1,mm))' - 2*X*X'; 
            % similarity (kernel) matrix 
            W = exp(-D/(2*sigma.^2)); 
            % ignore/ self-similarity
            W = W - diag(diag(W)); 

            % Compute D and normalize W by it:
            D = diag(sum(W,2).^-.5); 
            N = D*W*D;

            %
            % Find the eigenvectors corresponding to the K largest eigenvalues of
            % D^-1/2 W D^-1/2

            %Solve the eigensystem and find K largest eigenvectors:
            [Y Lambda] = eig(N); 
            [Lambda ii] = sort(diag(Lambda),1,'descend'); 
            Y = Y(:,ii(1:K));
            %	Form the matrix of eigenvectors and normalize again: 
            ss = sum(Y.^2,2).^-.5;
            Y = Y.*ss(:,ones(1,K));
            %Finally perform K-means: 
            [mu,cc] = UnSupervised.kmeans(Y,K,0);
        end %speclust

        %-------------------------------------------
        %> @brief Plots the kmeans function with the prototype, mu, plotted
        %> and an edge linking each datapoint to the prototype. Similar to
        %> some other functions here - it was forked from Adaptive
        %> modelling of complex data class.
        %>
        %> @param X [m x n] data points
        %> @param cc indicator variable (or cluster association variable)
        %> @param mu [K x n] the prototype values for each cluster
        %>
        %> 
        function plotkmeans(X,cc,mu)

            % plot result of K-means
            mm = size(X,1);
            K  = size(mu,1);

            col = colormap;
            col = col(ceil((.5:K)/K*size(col,1)),:);

            hh = ishold;
            for ii = 1:mm
              plot([X(ii,1) mu(cc(ii),1)],[X(ii,2) mu(cc(ii),2)],'-','linewidth',1,...
                  'color',[.7 .7 .7]);
              hold on
              plot(X(ii,1),X(ii,2),'.','markersize',16,'color',col(cc(ii),:));
            end
            for kk=1:K
              plot(mu(kk,1),mu(kk,2),'+','markersize',20,'linewidth',3,'color',col(kk,:));
            end
            if ~hh, hold off; end
        end %plotkmeans
        %----------------------------------------------------------------
        %> @brief plots the linkage of clusters. Forked from Adaptive
        %> modelling of complex data course.
        %>
        %> @param X [m x n] datapoints
        %> @param C 
        %> @param K [scalar] Number of clusters
        function showlinkage(X,C,K)

            [mm dd] = size(X);
            if dd ~= 2, error('Demo only works in 2D.'); end
            if K<2||K>mm, error('K has to be between 2 and the number of data items'); end

            say Showing 2 to =K clusters obtained by a linkage algorithm.
            say Press any key to continue
            for kk = 2:K
              % find indices of clusters
              top = 2*mm-1-kk+2:2*mm-1;
              ii = C.children(:,top);
              jj = find(ii<top(1));
              ii = ii(jj);
              % plot the flat clustering at kk clusters
              cc = zeros(1,mm);
              for ll = 1:kk
                cc(C.points{ii(ll)}) = ll;
              end
              plotclusters(X,cc);
              title([num2str(kk) ' clusters']);
              pause
            end

        end %showlinkage
        
        %----------------------------------------------------------
        
    end%static methods
end %classdef

