function [xbest,fitness_hist,counter_hist,fbest] =CMA(cost_function,x0,param)
% CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
% nonlinear function minimization.
%
% This code is an excerpt from cmaes.m and implements the key parts
% of the algorithm. It is intendend to be used for READING and
% UNDERSTANDING the basic flow and all details of the CMA *algorithm*.
% Computational efficiency is sometimes disregarded.

% -------------------- Initialization --------------------------------

% User defined input parameters (need to be edited)
N = length(x0); % number of objective variables/problem dimension
[row,col] = size(x0);
if row < col
    x0 = x0.';
end
xbest = x0; % best point overall
fbest =cost_function(xbest);
xmean = x0; % objective variables initial point
fitness_hist = [];
counter_hist = [];

linear_function_options = optimoptions('fmincon','Display','off');

count = 0;

if ~isfield(param,'sigma')
    sigma = 0.5; % coordinate wise standard deviation (step-size)
else
    sigma = param.sigma;
end
if ~isfield(param,'stopfitness')
    stopfitness = -inf; % stop if fitness < stopfitness (minimization)
else
    stopfitness = param.stopfitness;
end
if ~isfield(param,'stopeval')
    stopeval = 1e3*N^2; % stop after stopeval number of function evaluations
else
    stopeval = param.stopeval;
end
if ~isfield(param,'fitnesstol')
    fitnesstol =  1e-10; % stop after fitness improvement is lower than tolerance
else
    fitnesstol = param.fitnesstol;
end

if ~isfield(param,'messages')
    messages =  0; % display intermediate messages if on
else
    messages = param.messages;
end
if ~isfield(param,'stallbreak')
    stallbreak =  1; %break if stalling
else
    stallbreak = param.stallbreak;
end
if ~isfield(param,'parallelize')
    parallelize = 1; %useparfor instead of for
else
    parallelize = param.parallelize;
end


% Strategy parameter setting: Selection
if ~isfield(param,'lambda')
    lambda = 4+floor(3*log(N)); % population size, offspring number
else
    lambda = param.lambda;
end
if ~isfield(param,'mu')
    mu = lambda/2; % lambda=12; mu=3; weights = ones(mu,1); would be (3_I,12)-ES
else
    mu = param.mu;
end
if ~isfield(param,'UB')
    UB = +inf*ones(N,1); % Upper bound in the parameter space
else
    UB = param.UB;
end
if ~isfield(param,'LB')
    LB = -inf*ones(N,1); % Lower bound in the parameter space
else
    LB = param.LB;
end
if ~isfield(param,'bound_reduction')
    bound_reduction = 0.99; % rescaling the covariance if a bound is hit
else
    bound_reduction = param.bound_reduction;
end
if ~isfield(param,'variance_reduction_counter')
    variance_reduction_counter = 100; % number of rejection before rescaling occours
else
    variance_reduction_counter = param.variance_reduction_counter;
end

if ~isfield(param,'flat_fitness_terminate')
    flat_fitness_terminate = 1; % terminate in case of flat fitness
else
    flat_fitness_terminate = param.flat_fitness_terminate;
end


weights = log(mu+1/2)-log(1:mu).'; % muXone recombination weights
mu = floor(mu); % number of parents/points for recombination


weights = weights/sum(weights); % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu

% Strategy parameter setting: Adaptation
if ~isfield(param,'cc')
    cc = (4+mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
else
    cc = param.cc;
end
if ~isfield(param,'cs')
    cs = (mueff+2)/(N+mueff+5); % t-const for cumulation for sigma control
else
    cs = param.cs;
end
if ~isfield(param,'c1')
    c1 = 2 / ((N+1.3)^2+mueff); % learning rate for rank-one update of C
else
    c1 = param.c1;
end
if ~isfield(param,'cmu')
    cmu = 2 * (mueff-2+1/mueff) / ((N+2)^2+2*mueff/2); % and for rank-mu update
else
    cmu = param.cmu;
end
if ~isfield(param,'damps')
    damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
else
    damps = param.damps;
end


% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1); % evolution paths for C and sigma
if ~isfield(param,'B')
    B = eye(N); % B defines the coordinate system
else
    B = param.B;
end
if ~isfield(param,'D')
    D = eye(N); % diagonal matrix D defines the scaling
else
    D = param.D;
end
C = B*D*(B*D).'; % covariance matrix
% chol_C = chol(C,'lower');
chol_C = B*D;
eigeneval = 0; % B and D updated at counteval == 0
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2)); % expectation of
% ||N(0,I)|| == norm(randn(N,1))

% -------------------- Generation Loop --------------------------------

counteval = 0; % the next 40 lines contain the 20 lines of interesting code
fitness_diff = inf; % fitness difference
best_fitness_prev = inf;
while counteval < stopeval
    
    % Generate and evaluate lambda offspring
    arz = zeros(N,lambda);
    arx = zeros(N,lambda);
    arfitness = zeros(lambda,1);
    if parallelize
        parfor (k=1:lambda,lambda)
            in_bound = 0;
            
            while ~in_bound
                arz(:,k) = randn(N,1);
                arx(:,k) = xmean(:) + sigma * (chol_C *arz(:,k) ); % add mutation % Eq. 40
                if (sum(arx(:,k) > LB(:)) == N) && (sum(arx(:,k) < UB(:))==N)
                    in_bound =1;
                else
                    in_bound =0;
                end
            end
            
            
            arfitness(k) = cost_function(arx(:,k)); % objective function call
            counteval = counteval+1;
        end
    else
        for k=1:lambda
            in_bound = 0;
            Resample_count = 0;
            Resample_mask = true(N,1);
            Resample_variance_reduction = ones(N,1);
            LB_stats = zeros(N,1);
            UB_stats =zeros(N,1);
            A_stats = 0;
            while ~in_bound
                arz(Resample_mask,k) = randn(sum(Resample_mask),1);
                ary(Resample_mask) = xmean(Resample_mask) + sigma * (chol_C(Resample_mask,Resample_mask) *arz(Resample_mask,k) ).*Resample_variance_reduction(Resample_mask); % add mutation % Eq. 40
                LB_check = ary(:)  < LB(:);
                LB_stats = LB_stats + LB_check;
                UB_check = ary(:)  > UB(:);
                UB_stats = UB_stats + UB_check;
                Resample_mask = logical(LB_check + UB_check);
                A_stats = A_stats + 1;
                if (sum(LB_check) == 0) && (sum(UB_check)==0)
                    arx(:,k) = ary(:) ;
                    in_bound =1;
                else
                    Resample_count = Resample_count + 1;
                    if Resample_count == variance_reduction_counter % reset and decrease variance
                        Resample_count = 0; % reset counter
                        Resampling_dimensions = find(Resample_mask);
                        for i = Resampling_dimensions % loop throught all dimensions that weren't able to be sampled
                            
                             Resample_mask(i) = 0; % because we solved it
                        end
%                         Resample_variance_reduction = Resample_variance_reduction.*bound_reduction;
                         chol_C(Resample_mask,Resample_mask) =chol_C(Resample_mask,Resample_mask).*bound_reduction;
                            in_bound =1;
                    end
                    %                      disp(['Consider increasing upper bounds:', num2str(find(~UB_check.'))])
                    %                      disp(['Consider descreasing lower bounds:', num2str(find(~LB_check.'))])
                    %                         Resample_variance_reduction = Resample_variance_reduction.*bound_reduction;
                    %                      sigma = sigma*bound_reduction;
                end
            end
            
            
            
            arfitness(k) = cost_function(arx(:,k)); % objective function call
            counteval = counteval+1;
        end
    end
    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness); % minimization
    xmean = arx(:,arindex(1:mu))*weights; % recombination % Eq. 42
    zmean = arz(:,arindex(1:mu))*weights; % == D?-1*B?*(xmean-xold)/sigma
    fitness_diff = arfitness(1) - best_fitness_prev;
    
    % Cumulation: Update evolution paths
    ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean); % Eq. 43
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4+2/(N+1);
    
    pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (chol_C*zmean); % Eq. 45
    
    % Adapt covariance matrix C
    C = (1-c1-cmu) * C ... % regard old matrix % Eq. 47
        + c1 * (pc*pc.' ... % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction
        + cmu ... % plus ra [B,D] = eig(C); % eigen decomposition, B==normalized eigenvectorsnk mu update
        * (chol_C*arz(:,arindex(1:mu))) ...
        * diag(weights) * (chol_C*arz(:,arindex(1:mu))).';
    
    % Adapt step-size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); % Eq. 44
    
    % Update B and D from C
    if counteval - eigeneval > lambda/(c1+cmu)/N/10 % to achieve O(N?2)
        eigeneval = counteval;
         [B,D] = eig(C); % eigen decomposition, B==normalized eigenvectors
%         D = diag(sqrt(diag(D))); % D contains standard deviations now
        %         C=triu(C)+triu(C,1).'; % enforce symmetry
        C=(C + C.')/2; % enforce symmetry
        
    end
    
    
    
%     
    [chol_C_new,flag] = chol(C,'lower');
    if flag == 0
        chol_C = chol_C_new;
        
    end
  [B,D] = eig(C); % eigen decomposition, B==normalized eigenvectors
%         D = diag(sqrt(diag(D))); % D contains standard deviations now
%     chol_C = (B*D);
%     flag = 0;
% %      [B,D] = eig(C);
% %      chol_C2 = B*sqrt(D);
% % flag = 0;
%     if flag > 0
%          [B,D] = eig(C); % eigen decomposition, B==normalized eigenvectors
%         D = diag(sqrt(diag(D))); % D contains standard deviations now
%     chol_C = tril(B*D);
%         
%         
%     end
 if sum(isinf(C(:))) + sum(isnan(C(:)))> 0
     disp('Numerical problems. Terminating')
%        stopeval = countereval;
     return;
 end
    while flag >0 % in this case, we need to redefine the matrix as positive definite. repeat untiel we get a better matrix
%     ehile   
if sum(isinf(C(:))) + sum(isnan(C(:)))> 0
     disp('Numerical problems. Terminating')
%        stopeval = countereval;
     return;
 end
        [B,D] = eig(C);
        eigenva = diag(D);
        n_eigen =( eigenva<0);% extract the non positive eigenvalues
        z_eigen =( eigenva==0);
        new_eigenva = zeros(size(eigenva));
        null_eigenva = eps*ones(size(eigenva));
        new_eigenva(~n_eigen) = eigenva(~n_eigen);
        new_eigenva(n_eigen) = -eigenva(n_eigen);
        new_eigenva(z_eigen) = null_eigenva(z_eigen);
        D_new = diag(new_eigenva);
        C = real(B*D_new*inv(B));
        [chol_C,flag] = chol(C,'lower');
         if sum(isinf(C(:))) + sum(isnan(C(:)))> 0
     disp('Numerical problems. Terminating')
%      stopeval = countereval;
     return;
 end

        while flag > 0 
        
            [U,S,V] = svd(C);
              [B,D] = eig(C);
            [new_eigenva] = sort(diag(S));
            TOL = max(size(S)) * eps(norm(S));
            new_eigenva(new_eigenva < TOL) = TOL;
            D_new = diag(new_eigenva);
            C = real(B*D_new*inv(B));
            [chol_C,flag] = chol(C,'lower');
             if sum(isinf(C(:))) + sum(isnan(C(:)))> 0
     disp('Numerical problems. Terminating')
     return;
 end
        end
    end
    best_fitness_prev = arfitness(1);
    fitness_hist = [fitness_hist;best_fitness_prev];
    counter_hist = [counter_hist; counteval];
    % Break, if fitness is good enough
    if arfitness(1) <= stopfitness
        disp('Fitness has reached optimal level. Terminating');
        disp(['Total function evaluations: ', num2str(counteval)])
        break;
    end
    
    % Break, if fitness stop improving
    if abs(fitness_diff) < fitnesstol
        if stallbreak
            disp('Fitness improvement below threshold level. Terminating');
            disp(['Total function evaluations: ', num2str(counteval)])
            return;
        else
            disp('Fitness improvement below threshold level. Reshaping covariance');
            disp(['Total function evaluations: ', num2str(counteval)])
            sigma = sigma * exp(0.2+cs/damps);
        end
        
        
        
    end
    
    % Escape flat fitness, or better terminate?
    if arfitness(1) == arfitness(ceil(0.7*lambda))
        disp('Flat fitness, consider reformulating the objective:');
        if flat_fitness_terminate
            disp('Terminanting.')
            return;
        else
            disp('Reshaping covariance.')
            sigma = sigma * exp(0.2+cs/damps);
        end
        %
    end
    if messages
        disp([num2str(counteval) ': ' num2str(arfitness(1))]);
    end
    %     disp(['Fitness improvement:' , num2str(fitness_diff)])
    xmin = arx(:, arindex(1)); % Return best point of last generation.
    if arfitness(1) < fbest
        xbest = xmin;
        fbest = arfitness(1);
    end
%     if arfitness(1) - fbest > 3e3
%         disp('Numerical problem: Terminating');
%         break;
%     end
end % while, end generation loop

% -------------------- Final Message ---------------------------------



% Notice that xmean is expected to be even
% better.

% ---------------------------------------------------------------
function f=felli(x)
N = size(x,1); if N < 2 error('dimension must be greater one'); end
f=1e6.^((0:N-1)/(N-1)) * x.^2; % condition number 1e6


% build a smaller cost function, by slicing the original when needed

function cost = linear_cost_function(x,i,xi,cost_function)
            x(i) = xi;
%             disp(x);
            cost = cost_function(x);
            
