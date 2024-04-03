classdef Filter < handle
    %FILTER class that defines a Kalman based filter/smoother (Extended Kalman or Unscented Kalman)
    % notation is adopted from "Bayesian filtering and smoothing" By Simon
    % Sarkka, which I strongly reccomend for understanding the code
    % For particular applications, I strongly reccomend to hardcoded certain filters to gain speed and memory
    % performance. The objects of this class accepts a widely class of
    % state-space models but are quite slow. 
    
    % Process and measurement noise are ALWAYS ZERO MEAN
    properties
        type; % the type of the filter
        f; % process model
        h; % measurement model
    end
    
    methods
        function obj = Filter(type,process_model,measurement_model)
            %FILTER create a filter object
            % ----- INPUT -----
            % type: type of filter, as string format
            % process_model: process model 
            % measurement_model: measurement model
            % ----- OUTPUT -----
            % obj: the filter object
            %
            % The process and measurement models have to return only the function evaluated at the current point. Eventually the Jacobians in case we use the EKF
            % The accept in input the state value x, the noise value q or r, the parameters value theta, and the time index value t
            %
           
            switch type
                case 'EKF' % Extended Kalman Filter
                    obj.type = type;
                case 'UKF' % Unscented Kalman Filter
                    obj.type = type;    
                otherwise
                    error('Filter not yet implemented. Please read the class description.')
            end
            % we define in our obejct the process and measurement models
            obj.f = process_model;
            obj.h = measurement_model;
            
            
        end
        
        function [m,P,yf,E] = filtering(obj,m0,P0,y,theta,param)
            %FILTERING filter out the sequence of measurements using
            %Gaussian filtering for state space models. The posterior
            %probability of the hidden states given the measurements, i.e.
            %p(x|y) is given as Gaussian approximation  
            % Nx: dimension of the hidden state
            % Ny: dimension of the measurement space
            % Nt: total number of points
            % Nq: number of independent process noise sources
            % Nr: number of independent measurement noise sources
            % ----- INPUT -----
            % obj: the current filter, containing info about measurement
            % model, process model and type of filter implemented
            % m0: initial mean estmimate for the hidden state, of dimension 1*Nx
            % P0: initial covariance estimate for the hidden state, of dimension Nx*Nx
            % y: sequence of noisy measurements, of dimension Nt*Ny
            % theta: field structure of static parameters for the state-space model. Some parameters may be common to many models, such as
            %   - Q: process noise covariance matrix, Nq*Nq
            %   - R: measurement noise covariance matrix, Nr*Nr
            %   - dt: sampling time of the measurements
            % param: field structure of parameters for the filte, which can be:
            %   - delay: tells how many predicitons we do before applying the update step. This implies that the sampling time of the measurements is actually slower of a factor "delay" compared to the sampling time for the hidden process. Useful in case we are implementing a SDE as process model and for precision issues, we are forced to use a smaller timestep 
            %   - t0: initial time index of the first measurement sample
            %   - memory_constraint: for certain filters, it avoid to initialize the whole structure P since, when Nx is large, can lead to memory overflow. If enabled, increase the comp. speed but P is returned as a mean matrix (not as a temporal structure)
            %   - range: valid range for the hidden states; as simple implementation,  for hidden state i is range_i = [lower_bound_i; upper_bound_i]. If outside of this range (can happend during the update step) we force the state to goes back within its bounds. This approach in not optimal however, more complex approach can be seen in D.Simon "Kalman filtering with state constraints:a survey of linear and nonlinear algorithms"
            % ----- OUTPUT -----      
            % m: mean time series of the Gaussian posterior state estimate, of dimension Nt*Nx 
            % P: covariance estimate time serier of the Gaussian posterior state estimate, of dimension Nx*Nx*Nt. If param.memory_constraint == 1, P is the average matrix in time, with dimensions Nx*Nx
            % yf: filtered signal, obtained by applying the measurement equation to our mean estimate m, dimension Nt*Ny. If done properly, this signal has the measurement noise removed from the original y
            % E: log prediction error aka energy function, can be viewed as the unnormalized negative log-likelihood of the measurements given the parameters, of dimension Nt*1
            
            % parameter initiazliation
            if ~isfield(param,'delay')
                param.delay = 1; % deafult delay value
            end
            if ~isfield(param,'t0')
                param.t0 = 0; % global starting time
            end
            if ~isfield(param,'memory_constraint')
                param.memory_constraint = 0; % memory constraint. if positive, try to restrict the usage of memory
            end
            dtx = theta.dt/param.delay; % effective sampling time for the process model
            theta_proc = theta;
            theta_proc.dt = dtx; % we simply use the same theta for the process model except the sampling time
            % necessary to have R and Q as function handles
            if ~isa(theta_proc.Q, 'function_handle') 
                Q_temp = theta_proc.Q;
                theta_proc.Q = @(x,q,theta_proc,t)Q_temp;
                clear Q_temp;
            end
            if ~isa(theta.R, 'function_handle')
                R_temp = theta.R;
                theta.R = @(x,q,theta,t)R_temp;
                clear R_temp;
            end
            Nx = size(m0,2);
            if ~isfield(param,'range')
                param.range = repmat([-inf;+inf],1,Nx); % deafult range value, basically the Nx-hyperspace without bounds
            end
            Ny = size(y,2);
            % we assume the process noise size isn't larger than 2 times Nx
            % (broad assumption)
            Nq = length(theta_proc.Q(zeros(1,Nx),zeros(1,2*Nx),theta,1)); %
            Nr = length(theta.R(zeros(1,Nx),zeros(1,2*Ny),theta,1));
            Ntx = size(y,1)*param.delay*(param.delay > 1) + size(y,1)*(param.delay == 1);   % effective number of prediction steps
            Nty = size(y,1); % effective number of update steps
            
            m = zeros(Ntx,Nx);
            m(1,:) = m0;
            P_prev = P0;
            if ~param.memory_constraint
                P = zeros(Nx,Nx,Ntx);
                P(:,:,1) = P0;
            else
                P = P0; % it's going to be the mean value
            end
            yf = zeros(size(y));
            E = zeros(Nty,1);
            ty = 1;
            tx = 1;
            switch obj.type
                case 'EKF' % Extended Kalman filter
                    
                    while ty <= Nty
                        % predict
                        for d = 1:param.delay
                            tx = tx+1;
                            [fm_,Fm_,Fq_] = obj.f(m(tx-1,:),zeros(1,Nq),theta_proc,tx-1+param.t0); % f, jacobian respect to x and jacobian respect to q
                            m(tx,:) = fm_;
                            P_next = Fm_*P_prev*Fm_.' + ...
                                Fq_*theta_proc.Q(m(tx-1,:),zeros(1,Nq),theta_proc,tx-1+param.t0)*Fq_.';
                            % force symmetry - necessary to avoid numerical rounding errors
                            P_next = Filter.force_symmetry(P_next);
                            P_prev = P_next;
                            if ~param.memory_constraint
                                P(:,:,tx) = P_next;
                            end
                        end
                        % update
                        [hm,Hm,Hr] = obj.h(m(tx,:),zeros(1,Nr),theta,tx+param.t0);
                        v = y(ty,:) - hm;
                        S =  Hm*P_prev*Hm.' + ...
                            Hr*theta.R(m(tx,:),zeros(1,Nr),theta,tx+param.t0)*Hr.';
                        S_inv = eye(size(S))/S;
                        K = P_prev*Hm.'*S_inv;
                        mu = ( m(tx,:).' + K*v.').';
                        Pu = P_prev - K*S*K.';
                        % condition to keep our hidden state inside the desired range
                        if sum(mu < param.range(1,:)) > 0 % lower than the lower bound
                            % so we constraint the state to be equal to the
                            % lower bound
                            c_states = mu < param.range(1,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(1,c_states);
                            Pc = Pu;
                        elseif sum(mu > param.range(2,:))>0 % higher than the upper bound
                            c_states = mu > param.range(2,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(2,c_states);
                            Pc = Pu;
                        else % we are all good
                            mc = mu;
                            Pc = Pu;
                        end
                        m(tx,:) = mc;
                        P_next = Pc;
                        % force symmetry
                        P_next = Filter.force_symmetry(P_next);
                        
                        if ~param.memory_constraint
                            P(:,:,tx) = P_next;
                        else
                            P = P + P_next;
                        end
                        P_prev = P_next;
                        try
                            yf(ty,:) = obj.h(m(tx,:),zeros(1,Nr),theta,tx);
                        catch
                            [yf(ty,:),~,~] = obj.h(m(tx,:),zeros(1,Nr),theta,tx);
                        end
                        E(ty) = 0.5*log(det(2*pi*S)) + 0.5*v*S_inv*v.';
                        ty = ty+1;
                    end
                    if param.memory_constraint
                        P = P./Nty; % get back the mean value
                    end           
                case 'UKF' % Unscented Kalman filter
                    % the following parameters are for controlling the spread of the sigma-set points. Different values may lead to more accurate sigma-point representations
                    if ~isfield(param,'kappa')
                        param.kappa = 0;
                    end
                    if ~isfield(param,'alpha')
                        param.alpha = 1e-3;
                    end
                    if ~isfield(param,'beta')
                        param.beta = 2;
                    end
                    % parameters initialization
                    Nprime = Nx + Nq;
                    Nsecond = Nx + Nr;
                    x_sigma = cell(2*Nprime+1,1);
                    y_sigma = cell(2*Nsecond+1,1);
                    lambda_ukf_prime = param.alpha^2*(Nprime + param.kappa) - Nprime;
                    lambda_ukf_second = param.alpha^2*(Nsecond + param.kappa) - Nsecond;
                    
                    W_m_prime = zeros(1,2*Nprime+1);
                    W_c_prime  = zeros(1,2*Nprime+1);
                    W_m_prime(1) = lambda_ukf_prime/(Nprime + lambda_ukf_prime);
                    W_c_prime(1) = W_m_prime(1) + (1 - param.alpha^2 + param.beta);
                    W_m_prime(2:end) = 1/2/(Nprime + lambda_ukf_prime);
                    W_c_prime(2:end) = 1/2/(Nprime + lambda_ukf_prime);
                    
                    W_m_second = zeros(1,2*Nsecond+1);
                    W_c_second  = zeros(1,2*Nsecond+1);
                    W_m_second(1) = lambda_ukf_second/(Nsecond + lambda_ukf_second);
                    W_c_second(1) = W_m_second(1) + (1 - param.alpha^2 + param.beta);
                    W_m_second(2:end) = 1/2/(Nsecond + lambda_ukf_second);
                    W_c_second(2:end) = 1/2/(Nsecond + lambda_ukf_second);
                    
                    while ty <= Nty
                        % predict
                        for d = 1:param.delay
                            tx = tx+1;
                            %1. Generate the sigma points
                            P_aug = [P(:,:,tx-1), zeros(Nx,Nq); zeros(Nq,Nx), theta_proc.Q(m(tx-1,:),zeros(1,Nq),theta_proc,tx-1+param.t0)];%
                            P_sqrt = (sqrtm(P_aug));
                            x_sigma{1} =  [m(tx-1,:), zeros(1,Nq)];
                            for i = 2:Nprime+1 % loop for sigma points
                                x_sigma{i} =  x_sigma{1} + sqrt(Nprime + lambda_ukf_prime)*P_sqrt(:,i-1).';
                            end
                            for i = Nprime+2:2*Nprime+1
                                x_sigma{i} =  x_sigma{1} - sqrt(Nprime + lambda_ukf_prime)*P_sqrt(:,i-1-Nprime).';
                            end
                            
                            %2. Propagate the points throught the dynamic model
                            for  i=1:2*Nprime+1
                                x_sigma{i} = obj.f(x_sigma{i}(1:Nx),x_sigma{i}(Nx+1:end),theta_proc,tx-1+param.t0);
                            end
                            
                            
                            %3. Compute the predicted mean and covariance
                            mx = zeros(1,Nx);
                            Px = zeros(Nx);
                            for i = 1:2*Nprime+1
                                mx = mx + W_m_prime(i)*x_sigma{i};
                            end
                            
                            for i = 1:2*Nprime+1
                                Px = Px + W_c_prime(i)*(x_sigma{i} - mx).'*(x_sigma{i} - mx);
                            end
                            m(tx,:) = mx;
                            P(:,:,tx) = Px;
                            % force symmetry
                            P(:,:,tx) = Filter.force_symmetry(P(:,:,tx));
                            
                        end
                        
                        % update
                        P_aug = [P(:,:,tx), zeros(Nx,Nr); zeros(Nr,Nx), theta.R(m(tx,:),zeros(1,Nr),theta,tx+param.t0)];
                        P_sqrt = (sqrtm(P_aug));
                        x_sigma{1} = [m(tx,:), zeros(1,Nr)];
                        for i = 2:Nsecond+1 % loop for sigma points
                            x_sigma{i} =  x_sigma{1} + sqrt(Nsecond + lambda_ukf_second)*P_sqrt(:,i-1).';
                        end
                        for i = Nsecond+2:2*Nsecond+1
                            x_sigma{i} =  x_sigma{1} - sqrt(Nsecond + lambda_ukf_second)*P_sqrt(:,i-1-Nsecond).';
                        end
                        
                        
                        for  i=1:2*Nsecond+1
                            y_sigma{i} = obj.h(x_sigma{i}(1:Nx),x_sigma{i}(Nx+1:end),theta,tx+param.t0);
                        end
                        
                        my = zeros(1,Ny);    % predicted mean
                        S = zeros(Ny);  % predicted covariance of measurements
                        C = zeros(Nx,Ny); % cross covariance between the model and the measurements
                        for i = 1:2*Nsecond+1
                            my = my + W_m_second(i)*y_sigma{i};
                        end
                        for i = 1:2*Nsecond+1
                            S = S + W_c_second(i)*(y_sigma{i} - my).'*(y_sigma{i} - my);
                            C = C + W_c_second(i)*(x_sigma{i}(1:Nx) - mx).'*(y_sigma{i} - my);
                        end
                        
                        v = y(ty,:) - my;
                        S_inv = eye(size(S))/S;
                        K = C*S_inv;
                        mu = ( m(tx,:).' + K*v.').';
                        Pu = P(:,:,tx) - K*S*K.';
                        % either one of the two will verify
                        if sum(mu < param.range(1,:)) > 0 % lower than the lower bound
                            % so we constraint the state to be equal to the
                            % lower bound
                            c_states = mu < param.range(1,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(1,c_states);
                            Pc = Pu;
                        elseif sum(mu > param.range(2,:))>0 % higher than the upper bound
                            c_states = mu > param.range(2,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(2,c_states);
                            Pc = Pu;
                        else % we are all good
                            mc = mu;
                            Pc = Pu;
                        end
                        m(tx,:) = mc;
                        P(:,:,tx) = Pc;
                        % force symmetry
                        P(:,:,tx) = Filter.force_symmetry(P(:,:,tx));
                        yf(ty,:) = obj.h(m(tx,:),zeros(1,Nr),theta,tx+param.t0);
                        E(ty) = 0.5*log(det(2*pi*S)) + 0.5*v*S_inv*v.';
                        ty = ty+1;
                    end
            end
            m = m(2:end,:); % reshape the state sequence
            if ~param.memory_constraint
                P = P(:,:,2:end); %
            end
        end
        
        function [ms,Ps,ys,Gs] = smoothing(obj,m,P,theta,param)
            %SMOOTHING smooths out the sequence of measurements using
            %Gaussian smoothing for state space models. The posterior
            %probability of the hidden states given the measurements, i.e.
            %p(x|y) is given as Gaussian approximation. The difference with filtering is that p(x|y) is the posterior given all measurements, not just up to the current time index
            % for this reason smoothing is done backward and after a prior fitlering of the sequence
            % Nx: dimension of the hidden state
            % Nt: total number of points
            % Nq: number of independent process noise sources
            % ----- INPUT -----
            % obj: the current filter/smoother, containing info about measurement model, process model and type of filter implemented
            % m: sequence of filtered means, of dimension Nt*Nx
            % m: sequence of filtered covariance matrices, of dimension Nx*Nx*Nt
            % theta: field structure of static parameters for the state-space model. Some parameters may be common to many models, such as
            %   - Q: process noise covariance matrix, Nq*Nq
            %   - dt: sampling time of the measurements
            % param: field structure of parameters for the filte, which can be:
            %   - delay: tells how many predicitons we do before applying the update step. This implies that the sampling time of the measurements is actually slower of a factor "delay" compared to the sampling time for the hidden process. Useful in case we are implementing a SDE as process model and for precision issues, we are forced to use a smaller timestep 
            %   - t0: initial time index of the first measurement sample
            %   - memory_constraint: for certain filters, it avoid to initialize the whole structure P since, when Nx is large, can lead to memory overflow. If enabled, increase the comp. speed but P is returned as a mean matrix (not as a temporal structure)
            %   - range: valid range for the hidden states; as simple implementation,  for hidden state i is range_i = [lower_bound_i; upper_bound_i]. If outside of this range (can happend during the update step) we force the state to goes back within its bounds. This approach in not optimal however, more complex approach can be seen in D.Simon "Kalman filtering with state constraints:a survey of linear and nonlinear algorithms"
            % ----- OUTPUT -----      
            % ms: smoothed mean time series of the Gaussian posterior state estimate, of dimension Nt*Nx 
            % Ps:  smoothed covariance estimate time serier of the Gaussian posterior state estimate, of dimension Nx*Nx*Nt. If param.memory_constraint == 1, P is the average matrix in time, with dimensions Nx*Nx
            % ys: smoothed signal, obtained by applying the measurement equation to our mean estimate ms, dimension Nt*Ny. If done properly, this signal has the measurement noise removed from the original y. 
            % Gs: collection of matrices of Smoother Gains, of dimension Nx*Nx*Nt. The can be later used in the Expectataion maximization algorithm
            
            if ~isa(theta.Q, 'function_handle')
                Q_temp = theta.Q;
                theta.Q = @(x,q,theta,t)Q_temp;
                clear Q_temp;
            end
            if ~isa(theta.R, 'function_handle')
                R_temp = theta.R;
                theta.R = @(x,q,theta,t)R_temp;
                clear R_temp;
            end
            Nt = size(m,1);
            Nx = size(m,2);
            if ~isfield(param,'range')
                param.range = repmat([-inf;+inf],1,Nx); % deafult range value
            end
            if ~isfield(param,'t0')
                param.t0 = 0;
            end
            Nq = length(theta.Q(zeros(1,Nx),zeros(1,Nt),theta,1));
            Nr = length(theta.R(zeros(1,Nx),zeros(1,Nt),theta,1));
            try
                y_temp = obj.h(m(end,:),zeros(1,Nr),theta,Nt);
            catch
                [y_temp,~,~] = obj.h(m(end,:),zeros(1,Nr),theta,Nt);
            end
            Ny = length(y_temp);
            
            ms = zeros(Nt,Nx);
            Ps = zeros(Nx,Nx,Nt);
            Gs = zeros(Nx,Nx,Nt);
            ys = zeros(Nt,Ny);
            t = Nt-1;
            ms(t+1,:) = m(t+1,:);
            Ps(:,:,t+1) = P(:,:,t+1);
            ys(t+1,:) = y_temp;
            
            switch obj.type
                case 'EKF' % Extended Kalman filter
                    while t > 0
                        % predict
                        [fm_,Fm_,Fq_] = obj.f(m(t,:),zeros(1,Nq),theta,t+param.t0);
                        m_ = fm_;
                        P_ = Fm_*P(:,:,t)*Fm_.' + ...
                            Fq_*theta.Q(m(t,:),zeros(1,Nq),theta,t+param.t0)*Fq_.';
                        % force symmetry
                        P_ = Filter.force_symmetry(P_);
                        % update
                        
                        G = P(:,:,t)*Fm_.'/P_;
                        
                        mu = m(t,:) + (G*(ms(t+1,:) - m_).').';
                        Pu = P(:,:,t) + G*(Ps(:,:,t+1) - P_)*(G.');
                        if sum(mu < param.range(1,:)) > 0 % lower than the lower bound
                            % so we constraint the state to be equal to the
                            % lower bound
                            c_states = mu < param.range(1,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(1,c_states);
                            Pc = Pu;
                        elseif sum(mu > param.range(2,:))>0 % higher than the upper bound
                            c_states = mu > param.range(2,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(2,c_states);
                            Pc = Pu;
                        else % we are all good
                            mc = mu;
                            Pc = Pu;
                        end
                        ms(t,:) = mc;
                        Ps(:,:,t) = Pc;
                        
                        % force symmetry
                        Ps(:,:,t) = Filter.force_symmetry(Ps(:,:,t));
                        [ys(t,:),~,~] = obj.h(ms(t,:),zeros(1,Nr),theta,t+param.t0);
                        Gs(:,:,t) = G;
                        
                        t = t-1;
                    end
                case 'UKF'
                    if ~isfield(param,'kappa')
                        param.kappa = 0;
                    end
                    if ~isfield(param,'alpha')
                        param.alpha = 1e-3;
                    end
                    if ~isfield(param,'beta')
                        param.beta = 2;
                    end
                    Nprime = Nx + Nq;
                    x_sigma = cell(2*Nprime+1,1);
                    lambda_ukf_prime = param.alpha^2*(Nprime + param.kappa) - Nprime;
                    
                    W_m_prime = zeros(1,2*Nprime+1);
                    W_c_prime  = zeros(1,2*Nprime+1);
                    W_m_prime(1) = lambda_ukf_prime/(Nprime + lambda_ukf_prime);
                    W_c_prime(1) = W_m_prime(1) + (1 - param.alpha^2 + param.beta);
                    W_m_prime(2:end) = 1/2/(Nprime + lambda_ukf_prime);
                    W_c_prime(2:end) = 1/2/(Nprime + lambda_ukf_prime);
                    
                    while t > 0
                        %1. Generate the sigma points
                        P_aug = [P(:,:,t), zeros(Nx,Nq); zeros(Nq,Nx), theta.Q(m(t,:),zeros(1,Nq),theta,t+param.t0)];%
                        P_sqrt = (sqrtm(P_aug));
                        x_sigma{1} =  [m(t,:), zeros(1,Nq)];
                        for i = 2:Nprime+1 % loop for sigma points
                            x_sigma{i} =  x_sigma{1} + sqrt(Nprime + lambda_ukf_prime)*P_sqrt(:,i-1).';
                        end
                        for i = Nprime+2:2*Nprime+1
                            x_sigma{i} =  x_sigma{1} - sqrt(Nprime + lambda_ukf_prime)*P_sqrt(:,i-1-Nprime).';
                        end
                        x_sigma_ = x_sigma;
                        %2. Propagate the points throught the dynamic model
                        for  i=1:2*Nprime+1
                            x_sigma{i} =  obj.f(x_sigma{i}(1:Nx),x_sigma{i}(Nx+1:end),theta,t+param.t0);
                        end
                        
                        %3. Compute the predicted mean, covariance and cross-covriance
                        mx = zeros(1,Nx);
                        Px = zeros(Nx);
                        D = zeros(Nx);
                        for i = 1:2*Nprime+1
                            mx = mx + W_m_prime(i)*x_sigma{i};
                        end
                        
                        for i = 1:2*Nprime+1
                            Px = Px + W_c_prime(i)*(x_sigma{i} - mx).'*(x_sigma{i} - mx);
                            D = D + W_c_prime(i)*(x_sigma_{i}(1:Nx) - m(t,:)).'*(x_sigma{i} - mx);
                        end
                        % force symmetry
                        Px = Filter.force_symmetry(Px);
                        
                        % update
                        
                        G = D/Px;
                        mu = m(t,:) + (G*(ms(t+1,:) - mx).').';
                        Pu = P(:,:,t) + G*(Ps(:,:,t+1) - Px)*(G.');
                        if sum(mu < param.range(1,:)) > 0 % lower than the lower bound
                            % so we constraint the state to be equal to the
                            % lower bound
                            c_states = mu < param.range(1,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(1,c_states);
                            Pc = Pu;
                        elseif sum(mu > param.range(2,:))>0 % higher than the upper bound
                            c_states = mu > param.range(2,:); %  check the states that have a problem
                            mc =  mu;
                            mc(c_states) = param.range(2,c_states);
                            Pc = Pu;
                        else % we are all good
                            mc = mu;
                            Pc = Pu;
                        end
                        ms(t,:) = mc;
                        Ps(:,:,t) = Pc;
                        ys(t,:) = obj.h(ms(t,:),zeros(1,Nr),theta,t+param.t0);
                        Gs(:,:,t) = G;
                        % force symmetry
                        Ps(:,:,t) = Filter.force_symmetry(Ps(:,:,t));
                        Gs(:,:,t) = Filter.force_symmetry(Gs(:,:,t));
                        t = t-1;
                    end
            end
        end
        function [LB,m0s,P0s,Qs,Rs,m,P,yf,E,ms,Ps,ys,As,Hs] = EM(obj,m0,P0,y,theta,param)
            %EM expectation maximization algortihm for state space models. Only works with EKF based filter/smoother and with additive noise models
            % It can estimate the linear static parameters of the model, as well as the intial guess of the hidden sequence (m0,P0)

            % ----- INPUT -----
            % obj: the current filter/smoother, containing info about measurement model, process model and type of filter implemented
            % m0: initial mean estimate (for filtering)
            % p0: initial covariance estimate (for filtering)
            % theta: field structure of static parameters for the state-space model. Some parameters may be common to many models, such as
            %   - Q: process noise covariance matrix, Nq*Nq
            %   - R: measurement noise covariance matrix, Nr*Nr
            %   - A: process model matrix (for linear SSM)
            %   - H: measurement model matrix (for linear SSM)
            %   - dt: sampling time of the measurements
            % param: field structure of parameters for the filte, which can be:
            %   - delay: tells how many predicitons we do before applying the update step. This implies that the sampling time of the measurements is actually slower of a factor "delay" compared to the sampling time for the hidden process. Useful in case we are implementing a SDE as process model and for precision issues, we are forced to use a smaller timestep 
            %   - t0: initial time index of the first measurement sample
            %   - memory_constraint: for certain filters, it avoid to initialize the whole structure P since, when Nx is large, can lead to memory overflow. If enabled, increase the comp. speed but P is returned as a mean matrix (not as a temporal structure)
            %   - range: valid range for the hidden states; as simple implementation,  for hidden state i is range_i = [lower_bound_i; upper_bound_i]. If outside of this range (can happend during the update step) we force the state to goes back within its bounds. This approach in not optimal however, more complex approach can be seen in D.Simon "Kalman filtering with state constraints:a survey of linear and nonlinear algorithms"
            %   - EM_range: valid range of samples used to estimate the parameters of EM. It can be used to avoid intial transients of filtering/smoothing that may lead to biases in parameter estimation. It is an heuristic feature
            % ----- OUTPUT -----      
            % LB: Lower bound that EM tries to maximize. Sometimes, due to numerical precision, it may be not so representative of the learning process. It is better to monitor the log-likelihood or the energy function
            % m0s: EM best estimate of the starting point m0
            % P0s: EM best estimate of the starting covariance matrix P0
            % Qs: EM best estimate of the process noise covariance matrix Q
            % Rs: EM best estimate of the measurement noise covariance matrix R
            % m: filtered mean time series of the Gaussian posterior state estimate, of dimension Nt*Nx 
            % P:  filtered covariance estimate time serier of the Gaussian posterior state estimate, of dimension Nx*Nx*Nt. If param.memory_constraint == 1, P is the average matrix in time, with dimensions Nx*Nx
            % yf: filtered signal, obtained by applying the measurement equation to our mean estimate ms, dimension Nt*Ny. If done properly, this signal has the measurement noise removed from the original y. 
            % E: Energy function that I suggest to monitor when learning with EM, the lower the better.
            % ms: smoothed mean time series of the Gaussian posterior state estimate, of dimension Nt*Nx 
            % Ps:  smoothed covariance estimate time serier of the Gaussian posterior state estimate, of dimension Nx*Nx*Nt. If param.memory_constraint == 1, P is the average matrix in time, with dimensions Nx*Nx
            % ys: smoothed signal, obtained by applying the measurement equation to our mean estimate ms, dimension Nt*Ny. If done properly, this signal has the measurement noise removed from the original y. 
            % As: EM best estimate of the process model matrix A
            % Hs: EM best estimate of the measurement model matrix H
            
            if ~strcmp(obj.type,'EKF')
                error('EM works only with EKF (for now)')
            end
            if ~isfield(param,'t0')
                param.t0 = 0;
            end
            if ~isfield(param,'f')
                param.f = obj.f;
            end
            if ~isfield(param,'h')
                param.h = obj.h;
            end
            % filtering and smoothing
            [m,P,yf,E] = obj.filtering(m0,P0,y,theta,param);
            [ms,Ps,ys,Gs] = obj.smoothing([m0;m],cat(3,P0,P),theta,param); 
            Nt = size(ms,1);
            Ny = size(y,2);
            Nx = size(m0,2);
            Nq = length(theta.Q);
            Nr = length(theta.R);
            [dummy,~,~] = param.h(m0,zeros(1,Nr),theta,1);
            Nyt = size(dummy,2);
            if ~isfield(theta,'A')
                theta.A = eye(Nx) ;
            end
            if ~isfield(theta,'H')
                theta.H = eye(Ny) ;
            end
            if ~isfield(param,'EM_range')
                param.EM_range = Nt:-1:2;
            end
            SIGMA = zeros(Nx);
            THETA = zeros(Nyt,Nyt);
            PHI = zeros(Nx);
            B = zeros(Ny,Nyt);
            C = zeros(Nx);
            D = zeros(Ny);
            for k = param.EM_range
                [fms_,Fms_,~] = param.f(ms(k-1,:),zeros(1,Nq),theta,k-1+param.t0);
                [hms,Hms,~] = param.h(ms(k,:),zeros(1,Nr),theta,k+param.t0);
                SIGMA = SIGMA + Ps(:,:,k) + ms(k,:).'*ms(k,:); 
                PHI = PHI + Fms_*Ps(:,:,k-1)*Fms_.'  + fms_.'*fms_; 
                THETA = THETA + Hms*Ps(:,:,k)*Hms.'  + hms.'*hms;
%                 hms
%                 y(k-1,:)
                B = B + y(k-1,:)*hms; %ok
                C = C + Ps(:,:,k)*Gs(:,:,k-1).'*Fms_.' + ms(k,:).'*fms_; % ok
                D = D + y(k-1,:).'* y(k-1,:); % ok
            end
            SIGMA = SIGMA/(length(param.EM_range));
            THETA = THETA/(length(param.EM_range));
            PHI = PHI/(length(param.EM_range));
            B = B/(length(param.EM_range));
            C = C/(length(param.EM_range));
            D = D/(length(param.EM_range));
            PHI = Filter.force_symmetry(PHI);
            THETA = Filter.force_symmetry(THETA);
            LB = real( -0.5*Filter.logdet(2*pi*P0,'chol') - length(param.EM_range)*0.5*Filter.logdet(2*pi*theta.Q,'chol') - length(param.EM_range)*0.5*Filter.logdet(2*pi*theta.R,'chol') ...
                - 0.5*trace(inv(P0)*(Ps(:,:,1) + (ms(1,:) - m0).'*(ms(1,:) - m0))) ...
                - (length(param.EM_range))*0.5*trace(inv(theta.Q)*(SIGMA - C*theta.A.' - theta.A*C.' + theta.A*PHI*theta.A.'))...
                - (length(param.EM_range))*0.5*trace(inv(theta.R)*(D - theta.H*B.' - B*theta.H.' + theta.H*THETA*theta.H.')));
            m0s = ms(1,:);
            P0s = Ps(:,:,1) + (ms(1,:) - m0).'*(ms(1,:) - m0);
            P0s = Filter.force_symmetry(P0s);
            Qs = SIGMA - C*theta.A.' - theta.A*C.' + theta.A*PHI*theta.A.';
            Rs = D - theta.H*B.' - B*theta.H.' + theta.H*THETA*theta.H.';
            As = C*(inv(PHI));
            Hs = B*(inv(THETA));
            Qs = Filter.force_symmetry(Qs);
            Rs = Filter.force_symmetry(Rs);
            % reshaping for consistent dimensions
            Ps = Ps(:,:,2:end);
            ys = ys(2:end,:);
            ms = ms(2:end,:);
        end
        function E = get_energy(obj,y,particle,fields,theta,dt,m0,P0,E0,range)
            %GET_ENERGY  simplified version of the filter to return just the energy (log-likelihood) of the current particle (point representing the new static parameter to be tested)
            % this method is best used in combination with sampling algorithms
            % ----- INPUT -----
            % obj: the filter
            % y: noisy measurements
            % particle: cointains the vector of static parameters to be evaluated in a filtering phase, with dimensions 1*Nd
            % fields: strings corresponding to the parameters stored in particle. In practice, [particle,fields] forms a table of value-keys parameters
            % theta,dt,m0,P0: same parameters of fitlering, see above
            % E0: initial energy. It is actually the prior on the current particle
            % range: range of valid values for particle; Best practice is to ensure we have a valid particle at input, making this redundant
            % ----- OUTPUT -----      
            % E: energy of the current particle
            if (sum(particle < range(1,:)) > 0) || (sum(particle > range(2,:)) > 0) % additional check on particle domain
                E = +inf;
            else
                theta = obj.particle_to_fields(particle,fields,theta);
                [~,~,~,Et] = obj.filtering(y,theta,dt,m0,P0);
                E = sum(Et)+E0;% E0 should be evaluated as the prior
                if isnan(E) || imag(E)~= 0 % in case we have nan values, or imaginary values, indicate a possible divergence of the fitler. To not risk, we set the energy to its maximum (+ infinity)
                    E = +inf;
                end
                clear theta;
            end
        end
        function [Et,dEt] = get_energy_gradient(obj,y,particle,fields,theta,dt,m0,P0,differentials)
            %GET_ENERGY_GRADIENT  similar to above but return also the gradient of the current particle. How? by using sensitivity equations:
            % in practice, the derivatives of the system w.r.t. the static parameters are also propagated in a similar fashion to what happens in the Extended Kalman filter
            % This may increase the computational time a lot, but at least it is a first order method
            % this method is best used in combination with sampling algorithms
            % ----- INPUT -----
            % obj,y,particle,fields,theta,dt,m0,P0: same as above
            % differentials: differentials of the system w.r.t. the static parameters
            % ----- OUTPUT -----      
            % E: energy of the current particle
            % dE: gradient of the current particle
             % we are going to use the 3rd dimension for the parameter
            % derivatives, in order to preserve the vector/matrix shape of
            % the elements
            % theoretically not complete. since we assume that R and Q can
            % vary with the states, therefore it would be more appropriate
            % that we apply the chain rule to find their derivatives wrt
            % the parameters (but here is never the case). So in this version we assumes that R and Q are fixed and does not depends on the static parameters.
            Nd = length(particle);
            theta = obj.particle_to_fields(particle,fields,theta);
            % still need to initialize the filter
            dtx = dt/obj.param.delay;
            if ~isa(theta.Q, 'function_handle')
                Q_temp = theta.Q;
                theta.Q = @(x,q,theta,t,dt)Q_temp;
                clear Q_temp;
            end
            if ~isa(theta.R, 'function_handle')
                R_temp = theta.R;
                theta.R = @(x,q,theta,t,dt)R_temp;
                clear R_temp;
            end
            Nq = length(theta.Q(zeros(1,3),zeros(1,3),theta,1,dtx));
            Nr = length(theta.R(zeros(1,3),zeros(1,2),theta,1,dt));
            Ny = size(y,2);
            Ntx = size(y,1)*obj.param.delay*(obj.param.delay ~= 0) + size(y,1)*(obj.param.delay == 0);
            Nty = size(y,1);
            Nx = size(m0,2);
            
            m = zeros(Ntx,Nx);
            dm = zeros(Ntx,Nx,Nd);
            P = zeros(Nx,Nx,Ntx);
            dP = zeros(Nx,Nx,Ntx,Nd);
            Et = zeros(Nty,1);
            dEt = zeros(Nty,Nd);
            m(1,:) = m0;
            P(:,:,1) = P0;
            ty = 1;
            tx = 1;
            switch obj.type
                
                case 'EKF'
                    while ty < Nty
                        % predict
                        for td = 1:obj.param.delay
                            m(tx+1,:) = obj.f(m(tx,:),zeros(1,Nq),theta,tx,dtx);
                            P(:,:,tx+1) = obj.param.Fx(m(tx,:),zeros(1,Nq),theta,tx,dtx)*P(:,:,tx)*obj.param.Fx(m(tx,:),zeros(1,Nq),theta,tx,dtx).' + ...
                                obj.param.Fq(m(tx,:),zeros(1,Nq),theta,tx,dtx)*theta.Q(m(tx,:),zeros(1,Nq),theta,tx,dtx)*obj.param.Fq(m(tx,:),zeros(1,Nq),theta,tx,dtx).';
                            
                            Fxx_temp = differentials.Fxx(m(tx,:),zeros(1,Nq),theta,tx,dt); % should be a 3D matrix. Each 2d slice is the derivative of the jacobian Fx wrt one element x
                            Fxt_temp = differentials.Fxt(m(tx,:),zeros(1,Nq),theta,tx,dt);
                            Fqx_temp = differentials.Fqx(m(tx,:),zeros(1,Nq),theta,tx,dt);
                            Fqt_temp = differentials.Fqt(m(tx,:),zeros(1,Nq),theta,tx,dt);
                            Ft_temp = differentials.Ft(m(tx,:),zeros(1,Nq),theta,tx,dt);
                            dQ_temp = differentials.dQ(m(tx,:),zeros(1,Nq),theta,tx,dt);
                            Fxxmt = zeros(Nx,Nx,Nd);
                            Fqxmt = zeros(Nx,Nq,Nd);
                            for d = 1:Nd
                                for i = 1:size(Fxx_temp,1)
                                    for j = 1:size(Fxx_temp,2)
                                        Fxxmt(i,j,d) = sum(squeeze(Fxx_temp(i,j,:)).'*dm(tx,:,d).'); % has it to be zero?
                                    end
                                    for j = 1:size(Fqx_temp,2)
                                        Fqxmt(i,j,d) = sum(squeeze(Fqx_temp(i,j,:)).'*dm(tx,:,d).'); % has it to be zero?
                                    end
                                end
                            end
                            
                            for d = 1:Nd
                                dm(tx+1,:,d) = obj.param.Fx(m(tx,:),zeros(1,Nq),theta,tx,dtx)*dm(tx,:,d).' + Ft_temp(:,d);
                                dP(:,:,tx+1,d) = (Fxxmt(:,:,d)+ Fxt_temp(:,:,d))*P(:,:,tx)*obj.param.Fx(m(tx,:),zeros(1,Nq),theta,tx,dtx).' ...
                                    + obj.param.Fx(m(tx,:),zeros(1,Nq),theta,tx,dtx)*dP(:,:,tx,d)*obj.param.Fx(m(tx,:),zeros(1,Nq),theta,tx,dtx).' ...
                                    + obj.param.Fx(m(tx,:),zeros(1,Nq),theta,tx,dtx)*P(:,:,tx)*( Fxxmt(:,:,d) + Fxt_temp(:,:,d)).' +...
                                    (Fqxmt(:,:,d) + Fqt_temp(:,:,d))*theta.Q(m(tx,:),zeros(1,Nq),theta,tx,dtx)*obj.param.Fq(m(tx,:),zeros(1,Nq),theta,tx,dtx).' +...
                                    obj.param.Fq(m(tx,:),zeros(1,Nq),theta,tx,dtx)*dQ_temp(:,:,d)*obj.param.Fq(m(tx,:),zeros(1,Nq),theta,tx,dtx).'+ ...
                                    obj.param.Fq(m(tx,:),zeros(1,Nq),theta,tx,dtx)*theta.Q(m(tx,:),zeros(1,Nq),theta,tx,dtx)*( Fqxmt(:,:,d) + Fqt_temp(:,:,d)).';
                            end
                            tx = tx+1;
                        end
                        
                        % update
                        v = y(ty,:) - obj.h(m(tx,:),zeros(1,Nr),theta,tx,dt);
                        Hxx_temp = differentials.Hxx(m(tx,:),zeros(1,Nr),theta,tx,dt); % should be a 3D matrix. Each 2d slice is the derivative of the jacobian Hx wrt one element x
                        Hxt_temp = differentials.Hxt(m(tx,:),zeros(1,Nr),theta,tx,dt);
                        Hrx_temp = differentials.Hrx(m(tx,:),zeros(1,Nr),theta,tx,dt);
                        Hrt_temp = differentials.Hrt(m(tx,:),zeros(1,Nr),theta,tx,dt);
                        Ht_temp = differentials.Ht(m(tx,:),zeros(1,Nr),theta,tx,dt);
                        Hxxmt = zeros(Ny,Nx,Nx);
                        Hrxmt = zeros(Ny,Nr,Nx);
                        S =  obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt)*P(:,:,tx)*obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt).' + ...
                            obj.param.Hr(m(tx,:),zeros(1,Nr),theta,tx,dt)*theta.R(m(tx,:),zeros(1,Nr),theta,tx,dt)*obj.param.Hr(m(tx,:),zeros(1,Nr),theta,tx,dt).';
                        S_inv = eye(size(S))/S;
                        K = P(:,:,tx)*obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt).'*S_inv;
                        m(tx,:) = ( m(tx,:).' + K*v.').';
                        P(:,:,tx) = P(:,:,tx) - K*S*K.';
                        Et(ty) = 0.5*log(det(2*pi*S)) + 0.5*v*S_inv*v.';
                        
                        for d = 1:Nd
                            for i = 1:size(Hxx_temp,1)
                                for j = 1:size(Hxx_temp,2)
                                    Hxxmt(i,j,d) = sum(squeeze(Hxx_temp(i,j,:)).'*dm(tx,:,d).');   % shouldn't be zero all times
                                end
                                for j = 1:size(Hrx_temp,2)
                                    Hrxmt(i,j,d) = sum(squeeze(Hrx_temp(i,j,:)).'*dm(tx,:,d).');     % shouldn't be zero all times
                                end
                            end
                        end
                        dv = zeros(1,Ny,Nd);
                        dS = zeros(Ny,Ny,Nd);
                        dR_temp = differentials.dR(m(tx,:),zeros(1,Nr),theta,tx,dt);
                        dK = zeros(Nx,Ny,Nd);
                        for d = 1:Nd
                            dv(:,:,d) = -obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt)*dm(tx,:,d).' + Ht_temp(:,d);% shouldn't be zero all times
                            dS(:,:,d) = (Hxxmt(:,:,d) + Hxt_temp(:,:,d))*P(:,:,tx)*obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt).' ...
                                + obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt)*dP(:,:,tx,d)*obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt).'...
                                + obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt)*P(:,:,tx)*(Hxxmt(:,:,d)+ Hxt_temp(:,:,d)).' ...
                                + (Hrxmt(:,:,d) + Hrt_temp(:,:,d))*theta.R(m(tx,:),zeros(1,Nr),theta,tx,dt)*obj.param.Hr(m(tx,:),zeros(1,Nr),theta,tx,dt).' ...
                                + obj.param.Hr(m(tx,:),zeros(1,Nr),theta,tx,dt)*dR_temp(:,:,d)*obj.param.Hr(m(tx,:),zeros(1,Nr),theta,tx,dt).' ...
                                + obj.param.Hr(m(tx,:),zeros(1,Nr),theta,tx,dt)*theta.R(m(tx,:),zeros(1,Nr),theta,tx,dt)*(Hrxmt(:,:,d) + Hrt_temp(:,:,d)).';
                            dK(:,:,d) = dP(:,:,tx,d)*obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt).'*S_inv + ...
                                P(:,:,tx)*(Hxxmt(:,:,d)+ Hxt_temp(:,:,d)).'*S_inv - ...
                                P(:,:,tx)*obj.param.Hx(m(tx,:),zeros(1,Nr),theta,tx,dt).'*S_inv*dS(:,:,d)*S_inv;
                            dm(tx,:,d) = (dm(tx,:,d).' + dK(:,:,d)*(v.') + K*dv(:,:,d).').'; % OK
                            dP(:,:,tx,d) =  dP(:,:,tx,d) - (dK(:,:,d)*S*K.' + K*dS(:,:,d)*K.' + K*S*dK(:,:,d).'); %OK
                            dEt(ty,d) = 0.5*trace(S_inv*dS(:,:,d)) + 0.5*v*S_inv*dv(:,:,d).' + 0.5*dv(:,:,d)*S_inv*v.' - 0.5*v*S_inv*dS(:,:,d)*S_inv*v.'; %OK
                        end
                        ty = ty+1;
                    end
                    
                    
                case 'UKF'
                    error('UKF has not implemented yet')
            end
        end
        
    end
    methods(Static)
        function [m,P,yf,E] = filtering_OFC(m0,P0,y,theta,param)
            %FILTERING_OFC EKF but specifically made to estimate the phase noise of Optical fequency combs. This is considerably faster than the version above, but it is hardcoded
            % ----- INPUT ----.
            % m0,P0,y,theta,param: see Filter.Filtering above
            % ----- OUTPUT -----      
            % m,P,yf,E: see Filter.Filtering above
            Nd = length(particle);
            theta = obj.particle_to_fields(particle,fields,theta);
            % still need to initialize the filter
            if ~isfield(param,'memory_constraint')
                param.memory_constraint = 0; % memory constraint. if positive, try to restrict the usage of memory
            end
            Nx = size(m0,2);
            % we assume the process noise size isn't larger than 2 times Nx
            % (borad assumption)
            Nt = size(y,1);
            m = zeros(Nt+1,Nx);
            m(1,:) = m0;
            m_prev = m0;
            P_prev = P0;
            if ~param.memory_constraint
                P = zeros(Nx,Nx,Nt+1);
                P(:,:,1) = P0;
            else
                P = P0; % dummy value
            end
            yf = zeros(size(y));
            E = zeros(Nt,1);
            for t = 1:Nt
                % predict
                P_prev = P_prev + theta.Q;
                P_prev = Filter.force_symmetry(P_prev);
                % update
                hm = (theta.H.*theta.a)*(cos((t)*theta.dt*theta.omega + m_prev)).';
                Hm = -(theta.H.*theta.a).*sin((t)*theta.dt*theta.omega +  m_prev);
                v = y(t,:) - hm;
                S =  Hm*P_prev*Hm.' + theta.R;
                S_inv = 1/S;
                K = P_prev*Hm.'*S_inv;
                m_prev = ( m_prev.' + K*v.').';
                P_prev = P_prev - K*S*K.';
                P_prev = Filter.force_symmetry(P_prev);
                m(t+1,:) = m_prev;
                if ~param.memory_constraint
                    P(:,:,t+1) = P_prev;
                else
                    P = P + P_prev;
                end
                yf(t,:) = (theta.H.*theta.a)*(cos((t)*theta.dt*theta.omega + m_prev)).';
                E(t) = 0.5*log((2*pi*S)) + 0.5*v*S_inv*v.';
            end
            if param.memory_constraint
                P = P./(Nt+1); % get back the mean value
            else
                P = P(2:end,:);
            end
            m = m(2:end,:);
            
        end
        function [ms,Ps,ys] = smoothing_OFC(m,P,theta,param)
            %SMOOTHING_OFC EKS (Extended Kalman smoother) but specifically made to estimate the phase noise of Optical fequency combs. This is considerably faster than the version above, but it is hardcoded
            % ----- INPUT ----.
            % m,P,theta,param: see Filter.Smoothing above
            % ----- OUTPUT -----      
            % ms,Ps,ys: see Filter.Smoothing above
            if ~isfield(param,'memory_constraint')
                param.memory_constraint = 0; % memory constraint. if positive, try to restrict the usage of memory
            end
            Nt = size(m,1);
            Nx = size(m,2);
            ms = zeros(Nt,Nx);
            ms(end,:)  = m(end,:);
            if ~param.memory_constraint
                Ps = zeros(Nx,Nx,Nt);
                Ps(:,:,end) = P(:,:,end);
            else
                if ismatrix(P)
                    Ps = P;
                    P_prev = P;
                else
                    Ps = P(:,:,end); % dummy value
                end
            end
            
            ys = zeros(Nt,1);
            ys(end,:) = (theta.H.*theta.a)*(cos(Nt*theta.dt*theta.omega + m(end,:))).';
            for t = (Nt-1):-1:1
                % predict
                if ~ismatrix(P)
                    P_prev = P(:,:,t);
                end
                m_ = m(t,:);
                P_ = P_prev + theta.Q;
                P_ = Filter.force_symmetry(P_);
                % update
                G = P_prev/P_;
                ms(t,:) = m_ + (G*(ms(t+1,:) - m_).').';
                if ~param.memory_constraint
                    Ps(:,:,t) = P_prev + G*(Ps(:,:,t+1) - P_)*(G.');
                else
                    Ps = P_prev + G*(Ps - P_)*(G.');
                end
                
                ys(t,:) =  (theta.H.*theta.a)*(cos(t*theta.dt*theta.omega + ms(t,:))).';
            end
            ys = ys(1:end-1,:);
        end
        function [m0s,P0s,Qs,Rs,m,P,yf,E,ms,Ps,ys] = EM_OFC(m0,P0,y,theta,param)
            %EM_OFC EM but specifically made to estimate the covariance noise matrices of Optical fequency combs. This is considerably faster than the version above, but it is hardcoded
            % ----- INPUT ----.
            % m0,P0,y,theta,param: see Filter.EM above
            % ----- OUTPUT -----      
            % m0s,P0s,Qs,Rs,m,P,yf,E,ms,Ps,ys: see Filter.EM above
            if ~isfield(param,'filtering_only')
                param.filtering_only = 0;
            end
            
            Nx = size(m0,2);
            T = size(y,1);
            m = zeros(T+1,Nx);
            m(1,:) = m0;
            m_prev = m0;
            P_prev = P0;
            if ~param.filtering_only
                P = zeros(Nx,Nx,T+1);
                P(:,:,1) = P0;
            else
                P = P0;
            end
            yf = zeros(size(y));
            E = zeros(T,1);
            for k = 1:T
                % predict
                P_prev = P_prev + theta.Q;
                P_prev = Filter.force_symmetry(P_prev);
                % update
                hm = (theta.H.*theta.a)*(cos((k)*theta.dt*theta.omega + m_prev)).';
                Hm = -(theta.H.*theta.a).*sin((k)*theta.dt*theta.omega +  m_prev);
                v = y(k,:) - hm;
        
                S =  Hm*P_prev*Hm.' + theta.R;
                S_inv = 1/S;

                K = P_prev*Hm.'*S_inv;
                m_prev = ( m_prev.' + K*v.').';
 
                P_prev = P_prev - K*S*K.';
                P_prev = Filter.force_symmetry(P_prev);
                m(k+1,:) = m_prev;
 
                if ~param.filtering_only
                    P(:,:,k+1) = P_prev;
                else
                    P = P + P_prev;
                end
                yf(k,:) = (theta.H.*theta.a)*(cos((k)*theta.dt*theta.omega + m_prev)).';
                E(k) = 0.5*log((2*pi*S)) + 0.5*v*S_inv*v.';
            end
            
            if ~param.filtering_only
                % smoothing
                THETA = 0;
                SIGMA = zeros(Nx);
                PHI = zeros(Nx);
                B = 0;
                C = zeros(Nx);
                D = sum(y.^2);
                k = T;
                
                ms = zeros(T+1,Nx);
                %             Ps = zeros(Nx,Nx,T+1);
%                 Ps = zeros(Nx,Nx);
                Ps = zeros(Nx);
                ys = zeros(T,1);
                
                ms_next = m(end,:);
                ms(end,:)  = m(end,:);
                Ps_next = P(:,:,k+1);
                %             Ps(:,:,end) = P(:,:,k+1);
                hms = (theta.H.*theta.a)*(cos((k)*theta.dt*theta.omega + ms_next)).';
                Hms = -(theta.H.*theta.a).*sin((k)*theta.dt*theta.omega +  ms_next);
%                 THETA = THETA + Hms*Ps_next*Hms.'  + hms.'*hms; % ok
%                 B = B + y(k,:)*hms; %ok
                while k > 0
                    ys(k,:) = hms;
                    % access memory
                    P_prev = P(:,:,k);
                    m_prev = m(k,:);
                    % predict
                    m_ = m_prev;
                    P_ = P_prev + theta.Q;
                    P_ = Filter.force_symmetry(P_);
                    % update
                    G = P_prev/P_;
                    ms_prev = ms_next;
                    ms_next = m_prev + (G*(ms_next - m_).').';
                    SIGMA = SIGMA + Ps_next + ms_prev.'*ms_prev; % ok
                    
                    C = C + Ps_next*G.' + ms_prev.'*ms_next; % ok
                    
                    ms(k,:) = ms_next;
                    Ps_next = P_prev + G*(Ps_next - P_)*(G.');
                    Ps_next = Filter.force_symmetry(Ps_next);
%                     Ps(:,:,k) = Ps_next;
                    Ps = Ps + Ps_next;
                    hms = (theta.H.*theta.a)*(cos((k)*theta.dt*theta.omega + ms_next)).';
                    Hms = -(theta.H.*theta.a).*sin((k)*theta.dt*theta.omega +  ms_next);
                    PHI = PHI + Ps_next  + ms_next.'*ms_next; % ok
%                     if k > 1
%                         THETA = THETA + Hms*Ps_next*Hms.'  + hms.'*hms; % ok
%                     end
                    THETA = THETA + Hms*Ps_next*Hms.'  + hms.'*hms; % ok
                    B = B + y(k,:)*hms; %ok
                    k = k-1;
                end
                THETA = THETA/(T);
                SIGMA = SIGMA/(T);
                PHI = PHI/(T);
                C = C/(T);
                D = D/(T);
                B = B/(T);
                m0s = ms(1,:);
                P0s = Ps_next + (m0s - m0).'*(m0s - m0);
                P0s = Filter.force_symmetry(P0s);
                Qs = SIGMA - C - C.' + PHI;
                Rs = D - B.' - B + THETA;
                Qs = Filter.force_symmetry(Qs);
                ms = ms(2:end,:);
                P = P(:,:,2:end);
                Ps = Ps./(T);
            else
                m0s = m0;
                P0s = P0;
                Qs = theta.Q;
                Rs = theta.R;
                ms = [];
                Ps = [];
                ys = [];
                P = P./T;
            end
            m = m(2:end,:);
        end
         function [m0s,P0s,Qs,Rs,As,Hs,m,P,yf,E,ms,Ps,ys,e,s,Ks,P_ev] = EM_OFC_ap(m0,P0,y,theta,param)
            %EM_OFC_AP EM but specifically made to estimate the covariance noise matrices of Optical fequency combs. This is considerably faster than the version above, but it is hardcoded
            % The difference from the pervious one is that we also estimate the amplitude noise of the lines
            % the first half vector are the amplitudes, the second half are the phases.
            % we must assure that we cancel the correlation induced between phases and amplitudes (is it a good thing?)
            % ----- INPUT ----.
            % m0,P0,y,theta,param: see Filter.EM above
            % ----- OUTPUT -----      
            % m0s,P0s,Qs,Rs,m,P,yf,E,ms,Ps,ys: see Filter.EM above
            % e: innovation - or prediction error - sequence
            % s: prediction error variance sequence
            % Ks: kalman gain sequences
            % P_ev: covariance coefficient evolutions
            
            if ~isfield(param,'filtering_only')
                param.filtering_only = 0;
            end
            if ~isfield(param,'monitor')
                param.monitor = 0;
            end
            % deblocking allows to decorrelate phases and amplitudes
            if ~isfield(param,'deblock')
                param.deblock = 0;
            end
            
            Nx = size(m0,2);
            p_indeces = [];
            for i  = 1:Nx
                p_indeces = [p_indeces; ((i-1)*Nx  + (i:Nx)).'];
            end

            T = size(y,1);
            if param.monitor
                e = zeros(T,1);
                s = zeros(T,1);
                Ks = zeros(Nx,T);
                P_ev = zeros(length(p_indeces),T);
            else
                e = [];
                s = [];
                Ks = [];
                P_ev = [];
            end
            m = zeros(T+1,Nx);
            m(1,:) = m0;
            m_prev = m0;
            P_prev = P0;
            if ~param.filtering_only
                P = zeros(Nx,Nx,T+1);
                P(:,:,1) = P0;
            else
                P = P0;
            end
            yf = zeros(size(y));
            E = zeros(T,1);
            for k = 1:T
                
                % predict
                m_prev = (theta.A*(m_prev.')).';
                P_prev = theta.A*P_prev*theta.A.' + theta.Q; % no need to perform deblocking because we have external control on the shape of Q
                P_prev = Filter.force_symmetry(P_prev);
                % update
                hm = (theta.H.*(theta.a + m_prev(1:Nx/2)))*(cos((k)*theta.dt*theta.omega + m_prev(Nx/2+1:end))).';
                Hm = [theta.H.*(theta.a).*cos((k)*theta.dt*theta.omega + m_prev(Nx/2+1:end)),-(theta.H.*(theta.a + m_prev(1:Nx/2))).*sin((k)*theta.dt*theta.omega +  m_prev(Nx/2+1:end))];
                v = y(k,:) - hm;
               
                S =  Hm*P_prev*Hm.' + theta.R;
                
                S_inv = 1/S;

                K = P_prev*Hm.'*S_inv;
                
                m_prev = ( m_prev.' + K*v.').';
                if param.monitor
                    e(k) = v;
                    s(k) = S;
                    Ks(:,k) = K;
                    P_ev(:,k) = P_prev(p_indeces);
                end
                
                P_prev = P_prev - K*S*K.';
                P_prev = Filter.force_symmetry(P_prev);
                % forcing deblocking in case
                if param.deblock
                   P_prev = [P_prev(1:Nx/2,1:Nx/2), zeros(Nx/2);zeros(Nx/2),P_prev(1+Nx/2:end,1+Nx/2:end)]; 
                end
                
                m(k+1,:) = m_prev;
 
                if ~param.filtering_only
                    P(:,:,k+1) = P_prev;
                else
                    P = P + P_prev;
                end
                yf(k,:) = (theta.H.*(theta.a + m_prev(1:Nx/2)))*(cos((k)*theta.dt*theta.omega + m_prev(Nx/2+1:end))).';
                E(k) = 0.5*log((2*pi*S)) + 0.5*v*S_inv*v.';
            end
            
            if ~param.filtering_only
                % smoothing
                THETA = 0;
                SIGMA = zeros(Nx);
                PHI = zeros(Nx);
                B = 0;
                C = zeros(Nx);
                D = sum(y.^2);
                k = T;
                
                ms = zeros(T+1,Nx);
                %             Ps = zeros(Nx,Nx,T+1);
                Ps = zeros(Nx,Nx);
                ys = zeros(T,1);
                
                ms_next = m(end,:);
                ms(end,:)  = m(end,:);
                Ps_next = P(:,:,k+1);
                %             Ps(:,:,end) = P(:,:,k+1);
                hms = (theta.H.*(theta.a + ms_next(1:Nx/2)))*(cos((k)*theta.dt*theta.omega + ms_next(Nx/2+1:end))).';
                Hms = [theta.H.*(theta.a).*cos((k)*theta.dt*theta.omega + ms_next(Nx/2+1:end)),-(theta.H.*(theta.a + ms_next(1:Nx/2))).*sin((k)*theta.dt*theta.omega +  ms_next(Nx/2+1:end))];
%                 THETA = THETA + Hms*Ps_next*Hms.'  + hms.'*hms; % ok
%                 hms
%                 y(k,:)
%                 B = B + y(k,:)*hms; %ok
                while k > 0 
                    
                    % access memory
                    P_prev = P(:,:,k);
                    m_prev = m(k,:);
                    % predict
                    m_ = (theta.A*(m_prev.')).';
                    P_ = theta.A*P_prev*theta.A.' + theta.Q;
                    P_ = Filter.force_symmetry(P_);
                    
                    % update
                    G = P_prev*(theta.A.')/P_;
                    ms_prev = ms_next;
                    ms_next = m_prev + (G*(ms_next - m_).').';
                    SIGMA = SIGMA + Ps_next + ms_prev.'*ms_prev; % ok
                    
                    C = C + Ps_next*G.' + ms_prev.'*ms_next; % ok
                    Ps(:,:,k) = Ps_next;
                    ms(k,:) = ms_next;
                    Ps_next = P_prev + G*(Ps_next - P_)*(G.');
                    Ps_next = Filter.force_symmetry(Ps_next);
                     if param.deblock
                         Ps_next = [Ps_next(1:Nx/2,1:Nx/2), zeros(Nx/2);zeros(Nx/2),Ps_next(1+Nx/2:end,1+Nx/2:end)]; 
                     end
                    hms = (theta.H.*(theta.a + ms_next(1:Nx/2)))*(cos((k)*theta.dt*theta.omega + ms_next(Nx/2+1:end))).';
                    Hms = [theta.H.*(theta.a).*cos((k)*theta.dt*theta.omega + ms_next(Nx/2+1:end)),-(theta.H.*(theta.a + ms_next(1:Nx/2))).*sin((k)*theta.dt*theta.omega +  ms_next(Nx/2+1:end))];
                    PHI = PHI + Ps_next  + ms_next.'*ms_next; % ok
                    if k > 1
                        THETA = THETA + Hms*Ps_next*Hms.'  + hms.'*hms; % ok  
%                         hms
%                         y(k,:)
                        
                    end
                    B = B + y(k,:)*hms; %ok
                    ys(k,:) = hms;
                    k = k-1;
                end
                THETA = THETA/(T);
                SIGMA = SIGMA/(T);
                PHI = PHI/(T);
                C = C/(T);
                D = D/(T);
                B = B/(T);
                m0s = ms(1,:);
                P0s = Ps_next + (m0s - m0).'*(m0s - m0);
                P0s = Filter.force_symmetry(P0s);
                Qs = SIGMA - C - C.' + PHI;
                Rs = D - B.' - B + THETA;
                As = C*(inv(PHI));
                Hs = B*(inv(THETA));
                Qs = Filter.force_symmetry(Qs);
                ms = ms(2:end,:);
                P = P(:,:,2:end);
            else
                m0s = m0;
                P0s = P0;
                Qs = theta.Q;
                Rs = theta.R;
                As = [];
                Hs = [];
                ms = [];
                Ps = [];
                ys = [];
                P =P./T;
            end
            m = m(2:end,:);
        end
        
        function [E,m,P,yf,theta] = filtering_rate_equations(m0,P0,y,theta,param)
            %EM_OFC EKF but specifically made to estimate the rate equation states of a laser diode. This is considerably faster than the version above, but it is hardcoded
            % applications of laser rate equations for Bayesian filtering can be found on: G brajato D zibar "Joint Learning of Laser Relative Intensity and Frequency Noise from Single Experiment and Single Detected Quadrature"
            % ----- INPUT ----.
            % m0,P0,y,theta,param: see Filter.filtering above. In this case, it is good to take a look at the Laser class to see how this differential equations are implemented
            % ----- OUTPUT -----      
            % E,m,P,yf,theta: see Filter.filtering above
            timestep_ratio = theta.dt_measure/theta.dt_predict; % has the same role of the delay featured n Filter.Filtering
            if floor(timestep_ratio) ~= timestep_ratio
                timestep_ratio = ceil(timestep_ratio);
                theta.dt_predict = theta.dt_measure/timestep_ratio;
            end
            Nx = size(m0,2);
            Nt = size(y,1);
            m = zeros(Nt+1,Nx);
            m(1,:) = m0;
            m_prev = m0;
            P_prev = P0;
            if ~param.memory_constraint
                P = zeros(Nx,Nx,Nt+1);
                P(:,:,1) = P0;
            else
                P = P0; % dummy value
            end
            yf = zeros(size(y));
            Et = zeros(Nt,1);
            for t = 1:Nt
                % predict
                for k = 1:timestep_ratio % loop executed during prediciton steps
                    F_m = [1 + (-1/theta.tau_n-theta.g*m_prev(2)/(1+ theta.eps*m_prev(2)))*theta.dt_predict ,... dfN/dN ok
                        (- theta.g*(m_prev(1) - theta.N0)/(1+ theta.eps*m_prev(2))^2 )*theta.dt_predict,... dfN/dS ok
                        0; ... dfN/dphi ok
                        ( theta.beta/theta.tau_n + theta.g*m_prev(2)/(1+ theta.eps*m_prev(2)))*theta.dt_predict ,... dfS/dN
                        1+ (theta.g*(m_prev(1) - theta.N0)/(1+ theta.eps*m_prev(2))^2 - 1/theta.tau_p)*theta.dt_predict, ... dfS/dS
                        0 ; ... dfS/dphi
                        theta.alpha/2*theta.g*theta.dt_predict,... dfphi/dN
                        0,...dfphi/dS
                        1]; % dfphi/dphi
                    Q =  [(2*m_prev(1)/theta.tau_n*(1+theta.beta*m_prev(2))), -2*theta.beta*m_prev(1)*m_prev(2)/theta.tau_n,0;
                        -2*theta.beta*m_prev(1)*m_prev(2)/theta.tau_n, 2*theta.beta*m_prev(1)*m_prev(2)/theta.tau_n, 0;
                        0,0, (theta.beta*m_prev(1)/(2*theta.tau_n*m_prev(2)))];
                    m_prev = m_prev + [theta.I/Const.q - m_prev(1)/theta.tau_n - theta.g*(m_prev(1) - theta.N0)*m_prev(2)/(1+theta.eps*m_prev(2)),...
                        theta.g*(m_prev(1) - theta.N0)*m_prev(2)/(1+theta.eps*m_prev(2)) - m_prev(2)/theta.tau_p + theta.beta*m_prev(1)/theta.tau_n,...
                        theta.alpha*theta.g/2*(m_prev(1) - theta.N_ss)].*theta.dt_predict;
                    P_prev = F_m*P_prev*F_m.' + Q.*sqrt(theta.dt_predict);
                    P_prev = Filter.force_symmetry(P_prev);
                end
                % update
                hm = theta.K*sqrt(m_prev(2)/theta.tau_p)*cos(theta.omega*t*theta.dt_measure + theta.D*(t*theta.dt_measure)^2 + m_prev(3));
                Hm = [0, theta.K/sqrt(theta.tau_p*m_prev(2))*cos(theta.omega*(t*theta.dt_measure) + theta.D*(t*theta.dt_measure)^2 + m_prev(3)),- theta.K*sqrt(m_prev(2)/theta.tau_p)*sin(theta.omega*(t*theta.dt_measure) + theta.D*(t*theta.dt_measure)^2 + m_prev(3))];
                v = y(t,:) - hm;
                S =  Hm*P_prev*Hm.' + theta.R;
                S_inv = 1/S;
                K = P_prev*Hm.'*S_inv;
                m_prev = ( m_prev.' + K*v.').';
                P_prev = P_prev - K*S*K.';
                P_prev = Filter.force_symmetry(P_prev);
                m(t+1,:) = m_prev;
                if ~param.memory_constraint
                    P(:,:,t+1) = P_prev;
                else
                    P = P + P_prev;
                end
                yf(t,:) = theta.K*sqrt(m_prev(2)/theta.tau_p)*cos(theta.omega*t*theta.dt_measure + theta.D*(t*theta.dt_measure)^2 + m_prev(3)); 
                Et(t) = 0.5*log((2*pi*S)) + 0.5*v*S_inv*v.';
            end
            if param.memory_constraint
                P = P./(Nt+1); % get back the mean value
            else
                P = P(2:end,:);
            end
            m = m(2:end,:);
            E = sum(Et);
        end
        function theta = particle_to_fields(particle,fields,theta)
            %PARTICLE_TO_FIELDS Converts the current particle (vector 1*Nd) with the associated fields (string 1*Nd) into a pre-existing structure theta, to being used in the Bayesian filter and smoother
            % the reason behind this is to emphasis the static parameters during filtering, but then have the flexibiliy to do mathematical operations on the parameter itself
            % ----- INPUT ----.
            % particle: vector of numerical values to be written in theta
            % fields: cell of strings, which contains the specific parameters associated with particle
            % theta: main parameter struct. 
            % ----- OUTPUT -----      
            % theta: new structure with the new parameters. Parameters that are not reported in fields will not be changed
            % ----- EXAMPLE -----
%             particle = randn(1,3); 
%             fields = {{'alpha'},{'1'};{'beta'},{'1'};{'gamma'},{'1'}};
%             theta = struct();
%             theta.alpha = 0;
%             theta.beta = 0;
%             theta.gamma = 0;
%             theta = Filter.particle_to_fields(particle,fields,theta);
            ind = 1;
            for i=1:size(fields,1)
                curr_field = fields{i,1};
                curr_field_elem = fields{i,2};
                eval([curr_field,' = theta.',curr_field,';'])
                for j = 1:numel(curr_field_elem) % write only the parameters that matters
                    curr_particle = particle(ind);
                    eval([curr_field,'(',num2str(curr_field_elem(j)),')=curr_particle;'])
                    ind = ind+1;
                end
                eval(['theta.',curr_field,' =', curr_field,';'  ]) % finally write in theta
                clear curr_field;
            end
        end
        function Ms = force_symmetry(M)
            Ms = (M + M')/2;
        end
        function [M,Mc] = moments(X,n)
            %MOMENTS return the empirical moments (sample moments, or raw moments) (up to n) of a set of particles (population). See https://en.wikipedia.org/wiki/Moment_(mathematics)#Sample_moments 
            % ----- INPUT ----.
            % X : population of dimension Nx * N. The shorter size is assumed to be the dimension, while the larger is assumed to be the population number
            % n : how many orders of moments we want to calculate
            % ----- OUTPUT -----      
            % M : moments
            % Mc : central moments
            %
            % ----- EXAMPLE -----
            % X = randn(2,1e5);
            % n = 3;
            % [M,Mc] = moments(X,n);
            [rows,cols] = size(X);
            if cols > rows
                X = X.';
            end
            N = length(X);
            M = cell(n,1);
            Mc = cell(n,1);
            for i = 1:n
                M{i} = 1/N*sum(X.^i);
                Mc{i} = 1/(N-1)*sum((X - mean(X)).^i);
            end
            
        end
        function v = logdet(A, op)
            %THIS IS NOT MY CODE. ALL CREDITS GOES TO DAHUA LIN, https://se.mathworks.com/matlabcentral/fileexchange/22026-safe-computation-of-logarithm-determinat-of-large-matrix?focused=5105640&tab=function
            %LOGDET logarithm of the determinant of A. It is more stable than log(det(A)), and should avoid numerical error for A large and with small eigenvalues.
            % ----- INPUT ----.
            % A: matrix
            % op: mode to use the cholesky decomposition of A
            % ----- OUTPUT -----      
            % v : log(det(A))
            % ----- EXAMPLE -----
            % B = randn(3); % random square matrix
            % v = logdet(B*B.')
            assert(isfloat(A) && ismatrix(A) && size(A,1) == size(A,2), ...
                'logdet:invalidarg', ...
                'A should be a square matrix of double or single class.');
            
            if nargin < 2
                use_chol = 0;
            else
                assert(strcmpi(op, 'chol'), ...
                    'logdet:invalidarg', ...
                    'The second argument can only be a string ''chol'' if it is specified.');
                use_chol = 1;
            end
            
            %% computation
            
            [chol_A,p] = chol(A);
            if p == 0
                use_chol = 1;
            else
                use_chol = 0;
            end
            if use_chol
                v = 2 * sum(log(diag(chol_A)));
            else
                [L, U, P] = lu(A);
                du = diag(U);
                c = det(P) * prod(sign(du));
                v = log(c) + sum(log(abs(du)));
            end
        end
        function y = brickwallFilter(x,Ts,f_limits,whitening_mode,noise_level)
        min_freq = f_limits(1);
        max_freq = f_limits(2);
        N = length(x);
        f_ext = (0:N-1)./N./Ts;
        nyquist_freq = 1/Ts/2;
        f_half = f_ext(f_ext <= nyquist_freq);
        BW_filter = zeros(N,1);
        BW_noise_conditioner = zeros(N,1);
        BW_signal_mask = (f_half > min_freq).*(f_half < max_freq);
        BW_signal_mask_ext =[BW_signal_mask;zeros(Nt-length(f),1)];
        BW_filter(logical(BW_signal_mask_ext)) = 1;
        BW_filter = [BW_filter(1:N/2);flipud(BW_filter(1:N/2))];
        BW_filter = (fft(real(ifft(BW_filter))));
        BW_noise_mask_ext = ~BW_signal_mask_ext;
        BW_noise_conditioner(logical(BW_noise_mask_ext)) = 1;
        BW_noise_conditioner = [BW_noise_conditioner(1:N/2);flipud(BW_noise_conditioner(1:N/2))];
        BW_noise_conditioner = (fft(real(ifft(BW_noise_conditioner))));
        % prefitlering or
        switch whitening_mode
            case 'off'
            y = (ifft(fft(x).*BW_filter));
            case 'classic'
            
            case 'normal'
            y = (ifft(fft(x).*(BW_filter + BW_noise_conditioner.*noise_level)));    
            case 'random'
            y = (ifft(fft(x).*(BW_filter + BW_noise_conditioner.*(noise_level + 0.1*randn(N,1)))));     
        end
        end
    end
end

