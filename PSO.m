function [min_point,min_cost,min_cost_hist,f_count_hist] = PSO(cost_function, param)
%PSO Summary of this function goes here
%   Detailed explanation goes here
% References: 
% [1] A Comprehensive Survey on Particle Swarm Optimization Algorithm and Its Applications
% [2] Standard Particle Swarm Optimisation
% [3] https://en.wikipedia.org/wiki/Particle_swarm_optimization


%% PSO optimizer

min_cost_hist=[];
f_count_hist=[];

if ~isfield(param,'w')
   param.w = 1/(2*log(2)); % momemntum 
end

if ~isfield(param,'c')
   param.c = 0.5 + log(2); % exploration 
end

if ~isfield(param,'maxiter')
   param.maxiter = 1e4; % max number of iterations 
end

N_d = size(param.bounds,1); % number of dimensions
if ~isfield(param,'N_p')
   param.N_p = 4+floor(3*log(N_d)); % swarm size 
end

if ~isfield(param,'maxfeval')
   param.maxfeval = (param.maxiter +1)*(param.N_p); % max number of function evaluations
end

f_count = 0;
N_d = size(param.bounds,2);
%% Initialization
x_i = zeros(param.N_p,N_d);
p_best= zeros(param.N_p,N_d);
v_i= zeros(param.N_p,N_d);
parfor i = 1 : param.N_p
    x_i(i,:) = param.bounds(1,:) + (param.bounds(2,:)-param.bounds(1,:)).*rand(1,N_d);
    p_best(i,:) = x_i(i,:);
    v_i(i,:) = param.bounds(1,:) - x_i(i,:) + (param.bounds(2,:)-x_i(i,:)- param.bounds(1,:)-x_i(i,:)).*rand(1,N_d);
end
parfor j = 1 : param.N_p % evaluate the cost of the entire swarm
    eval_fun(j) = cost_function(x_i(j,:));
    f_count = f_count + 1;
end
eval_fun_p = eval_fun;
[val,pos] = min(eval_fun);
g_best = x_i(pos,:);
eval_fun_g_best = val;
%% iterations
m = 1;

mask_lower = false(param.N_p,N_d);
mask_upper= false(param.N_p,N_d);
while  m < param.maxiter && f_count <param.maxfeval % number of iterations
   parfor k = 1 : param.N_p % do this update for each particle
        c = (param.c)*rand(2,1);
        v_i(k,:) = param.w*v_i(k,:) + c(1)*(p_best(k,:) - x_i(k,:))+ c(2)*(g_best - x_i(k,:)); 
        x_i(k,:) = x_i(k,:) + v_i(k,:); % update position for the current particle
        mask_lower(k,:) = x_i(k,:) < param.bounds(1,:);
        if  sum(mask_lower(k,:)) >= 1 % in case we hit the bounds, we bounce back the particle - but we should do so only for the interested dimensions - 
            x_i(k,:) = x_i(k,:).*(~mask_lower(k,:)) +  param.bounds(1,:).*(mask_lower(k,:));      
            v_i(k,:) =v_i(k,:).*(~mask_lower(k,:)) -0.5*v_i(k,:).*(mask_lower(k,:));
        end
        mask_upper(k,:) = x_i(k,:) > param.bounds(2,:);
        if  sum(mask_upper(k,:)) >= 1 % in case we hit the bounds, we bounce back the particle
            x_i(k,:) = x_i(k,:).*(~mask_upper(k,:)) + param.bounds(2,:).*(mask_upper(k,:));
            v_i(k,:) = v_i(k,:).*(~mask_upper(k,:)) -0.5*v_i(k,:).*(mask_upper(k,:));
        end
        eval_fun(k) = cost_function(x_i(k,:));
         f_count = f_count + 1;
    end
    for k = 1 : param.N_p 
        if eval_fun(k) < eval_fun_p(k)
            p_best(k,:) = x_i(k,:);
            eval_fun_p(k) = eval_fun(k);
        end
        if eval_fun(k) < eval_fun_g_best
            g_best =x_i(k,:);
            eval_fun_g_best = eval_fun(k);
        end
    end
    m = m+1;   
    min_cost_hist=[min_cost_hist;eval_fun_g_best];
    f_count_hist=[f_count_hist;f_count];
    disp([num2str(f_count) ': ' num2str(eval_fun_g_best)]);
end
%% Output
min_cost = eval_fun_g_best;
min_point = g_best;


end

