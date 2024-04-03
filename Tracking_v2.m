classdef Tracking_v2 < matlab.mixin.SetGetExactNames
    %TRACKING_V2 class to define the tracker obejct using Kalman based
    %filters. For more detailed explanation, check out the online book:
    %Bayesian filtering and smoothing by Simo Sarkka
    
    properties
        A; % state transition matrix - which is the same as the Jacobian for linear systems
        Q; % covariance matrix of the process noise
        R; % covariance matrix of the measurement noise
        x_ini; % initial mean of the hidden state
        P_ini; % initial covariance of the hidden state
        h; % measurement function (for Kalman filter and extended Kalman filter)
        H_k; % Jacobian of the measurement function with respect to the system state
        meas_fun; % measurement function (for Unscented Kalman filter)
        state_fun; % state transition function (for Unscented Kalman filter)
        alpha; % alpha parameter of the sigma set for Unscented Kalman filter
        beta; % beta parameter of the sigma set for Unscented Kalman filter
        kappa; % kappa parameter of the sigma set for Unscented Kalman filter
        version; % version to use of the Unscented Kalman filter
        df; % NOT used
        kernel_bw; % the kernel size for the Kernel Unscented Kalman filter
        step_size; % the step size for the Kernel Unscented Kalman filter
        amp; % NOT used
        %   Ts;
    end
    
    methods
        function obj = Tracking_v2(param)
            % initialize the tracking object
            obj.Q = param.Q;
            obj.R = param.R;
            obj.x_ini = param.x_ini;
            %            obj.Ts = param.Ts;
            if isfield(param,'df')
                obj.df = param.df;
            else
                obj.df = [];
            end
            if isfield(param,'A')
                obj.A = param.A;
            else
                obj.A = [];
            end
            if isfield(param,'P_ini')
                obj.P_ini = param.P_ini;
            else
                obj.P_ini =[];
            end
            if isfield(param,'h')
                obj.h = param.h;
            else
                obj.h = [];
            end
            if isfield(param,'H_k')
                obj.H_k = param.H_k;
            else
                obj.H_k = [];
            end
            if isfield(param,'state_fun') && isfield(param,'meas_fun')
                obj.state_fun = param.state_fun;
                obj.meas_fun = param.meas_fun;
            else
                obj.state_fun = [];
                obj.meas_fun = [];
            end
            if isfield(param,'alpha') && isfield(param,'beta') && isfield(param,'kappa')
                obj.alpha = param.alpha;
                obj.beta = param.beta;
                obj.kappa = param.kappa;
            else
                obj.alpha = [];
                obj.beta = [];
                obj.kappa = [];
            end
            if isfield(param,'version')
                obj.version = param.version;
            else
                obj.version = [];
            end
            if isfield(param,'kernel_bw')
                obj.kernel_bw = param.kernel_bw;
            else
                obj.kernel_bw = [];
            end
            if isfield(param,'step_size')
                obj.step_size = param.step_size;
            else
                obj.step_size = [];
            end
            
            if isfield(param,'amp')
                obj.amp = param.amp;
            else
                obj.amp = [];
            end
        end
        
        function [x_est, P_est, log_pred_error] = EKF(obj,y)
            % Extened Kalman Filter             
            
            % this code works only for linear process models, where the
            % transition matrix and its Jacobian coincide.
            % for more general process models, a transition function and a
            % jacobian matrix are two different parameters
            x_up = obj.x_ini;
            P_up = obj.P_ini;
            L = length(y);
            x_est = zeros(L,length(obj.x_ini));
            P_est = zeros(length(obj.x_ini),length(obj.x_ini),L);
            S_pred = eye(size(y,2));
            p_yk_yk_1 = zeros(L,1);
            
            for n = 1 : L
                % predict
                x_pred = obj.A*x_up; 
                P_pred = obj.A*P_up*obj.A.' + obj.Q;
                % update
                v =  y(n,:)' - obj.h(x_pred,n);
                S_up = obj.H_k(x_pred,n)*P_pred*obj.H_k(x_pred,n)' + obj.R;
                [~,p] = chol(S_up); % checking for matrix singularities - takes time, so maybe to speed things up this step could be removed?
                if p~=0
                   S_up = S_pred; % just use the matrix at the previous step if the current one is singular
                end
                
                S_up_inv = 1/S_up; % doing the inverse only once per time step
                K = (P_pred*obj.H_k(x_pred,n)')*S_up_inv;
                x_up =  x_pred + K*v;
                P_up =  P_pred - K*S_up*K';
                % negative log likelihood     
                p_yk_yk_1(n) = log(mvnpdf(y(n,:),obj.h(x_pred,n)',S_up));
                % memorize the estimated sequence
                x_est(n,:) = x_up;
                P_est(:,:,n) = P_up;
                S_pred = S_up;
            end
            log_pred_error = sum( p_yk_yk_1 );
        end
         function [xs_est, Ps_est, G_est] = EKS(obj,x_est,P_est)
            % Extened Kalman smoother             
            
            % this code works only for linear process models, where the
            % transition matrix and its Jacobian coincide.
            % for more general process models, a transition function and a
            % jacobian matrix are two different parameters
            L = length(x_est);           
            
            xs_est = zeros(L,length(obj.x_ini));
            Ps_est = zeros(length(obj.x_ini),length(obj.x_ini),L);
            G_est = zeros(length(obj.x_ini),length(obj.x_ini),L);
            
            xs_est(end,:) = x_est(end,:);
            Ps_est(:,:,end) = P_est(:,:,end);
            for n = L-1:-1:1
                % predict
                x_pred = (obj.A*x_est(n,:).').'; 
                P_pred = obj.A*P_est(:,:,n)*obj.A.' + obj.Q;
                % update
                v = xs_est(n+1,:) - x_pred;
                S = Ps_est(:,:,n+1) - P_pred;
                
                G = P_est(:,:,n)*obj.A.'/P_pred;
                x_up =  x_pred + (G*v.').';
                P_up =  P_pred + G*S*G';
                % memorize the estimated sequence
                xs_est(n,:) = x_up;
                Ps_est(:,:,n) = P_up;
                G_est(:,:,n) = G;
            end
        end
        function [x_est, P_est,log_pred_error] = KF(obj,y)
            % Kalman filter
            
            x_up = obj.x_ini;
            P_up = obj.P_ini;
            L = length(y);
            x_est = zeros(L,length(obj.x_ini));
            p_yk_yk_1 = zeros(L,1);
            if length(obj.x_ini)  ~= 1
                P_est = zeros(length(obj.x_ini),length(obj.x_ini),L);
            else
                P_est = zeros(length(obj.x_ini),L);
            end
            for n = 1 : L
                % predict
                x_pred = obj.A*x_up;
                P_pred = obj.A*P_up*obj.A.' + obj.Q;
                % update
                v =  y(n) - obj.H_k*x_pred;
                S = obj.H_k*P_pred*obj.H_k' + obj.R;
                K = (P_pred*obj.H_k')/S;
                x_up =  x_pred + K*v;
                P_up =  P_pred - K*S*K';
                % log prediction error
                p_yk_yk_1(n) = log(mvnpdf(y(n),obj.H_k*x_pred,S));
                x_est(n,:) = x_up;
                P_est(n,:) = P_up;
            end
            log_pred_error = sum( p_yk_yk_1 );
        end
        function [x_est,P_est,log_pred_error] = UKF(obj,y)
            %Unscented Kalman Filter
            
            if  ~iscolumn(y)
                y = y';
            end
            L = length(y);
            x_est = zeros(L,length(obj.x_ini));
            P_est = zeros(length(obj.x_ini),length(obj.x_ini),L);
            P = obj.P_ini;
            x = obj.x_ini;
            
            switch obj.version
                case 'MATLAB' % use the Matlab UKF
                    obj_ukf = unscentedKalmanFilter(obj.state_fun,obj.meas_fun,obj.x_ini,'ProcessNoise',obj.Q,'MeasurementNoise',obj.R);
                    obj_ukf.Alpha = obj.alpha;
                    obj_ukf.Beta = obj.beta;
                    obj_ukf.Kappa = obj.kappa;
                    for n = 1 : L
                        [PredictedState,PredictedStateCovariance] = predict(obj_ukf);
                        [CorrectedState,CorrectedStateCovariance] = correct(obj_ukf,y(n,:));
                        x_est(n,:) = CorrectedState;
                        P_est(:,:,n) = CorrectedStateCovariance;
                    end
                    log_pred_error = nan; % it is necessary to check matlab documentation to get a calculation of the negative log likelihood
                    
                case 'DTU' % use our version UKF
                    [row,col] = size(y);
                    
                    if col < row
                        y = y';
                        N = col;
                    end
                    
                    M = length(obj.x_ini);
                    
                    %sigma weights for the mean estimation
                    Wm_0 = 1 - M/(obj.alpha^2*(M+obj.kappa));
                    Wm_i = 1/(2*obj.alpha^2*(M+obj.kappa));
                    Wm = [Wm_0,repmat(Wm_i,1,2*M)];
                   
                    %sigma weights for the covariance estimation
                    Wc_0 = (2 - obj.alpha^2 + obj.beta) - M/(obj.alpha^2*(M+obj.kappa));
                    Wc_i = 1/(2*obj.alpha^2*(M+obj.kappa));
                    Wc = [Wc_0,repmat(Wc_i,1,2*M)];
                    
                    g = zeros(N,2*M+1); % sigma set for prediction
                    f = zeros(M,2*M+1); % sigma set for update
                    
                    p_yk_yk_1 = zeros(length(y),1);
                    
                    for n = 1 : length(y)
                        
                        % PREDICT
                        Psc = obj.alpha^2*(M + obj.kappa)*P;
                        sPsc = (sqrtm(Psc)); % doing the matrix squareroot (expensive computation) only once
                        dx_i1 = sPsc;
                        dx_i2 = -sPsc;
                        x_i = [x, x + dx_i1, x + dx_i2];
                        for i = 1 : 2*M+1 % propagation of the sigma set throught the process model
                            f(:,i) = obj.state_fun(x_i(:,i));
                        end
                        
                        fm = sum(Wm.*f,2);
                        P = Wc.*(f - fm)*(f - fm)' + obj.Q;
                        x = fm;
                        
                        % UPDATE
                        Psc = obj.alpha^2*(M + obj.kappa)*P;
                        sPsc = (sqrtm(Psc)); % doing the matrix squareroot (expensive computation) only once
                        dx_i1 = sPsc;
                        dx_i2 = -sPsc;
                        x_i = [x, x + dx_i1, x + dx_i2];
                        
                        for i = 1 : 2*M+1 % propagation of the sigma set throught the measurement model
                            g(:,i) = obj.meas_fun(x_i(:,i),n);
                        end
                        
                        ym = sum(Wm.*g,2);
                        Py = Wc.*(g - ym)*(g - ym)' + obj.R;
                        Pxy = Wc.*(x_i - x)*(g - ym)';
                        S_inv = 1/(Py);
                        K = Pxy*S_inv;
                        v = (y(:,n) - ym);
                        x = x + K*v;
                        P = P - K*Py*K';
                        
                        x_est(n,:) = x;
                        P_est(:,:,n) = P;
                        p_yk_yk_1(n) = 0.5*log(det(2*pi*Py)) + 0.5*v*S_inv*v.';
                        
                    end
                    log_pred_error = sum( p_yk_yk_1 );
                otherwise
                    error('method not implemented')
                    
            end
        end
        function [x_est,P_est] = K_UKF(obj,y,X_tr)
            L = length(y);
            e_k = zeros(1,L);
            e_k(1) = X_tr(2);
            y_pred = zeros(1,L);
            x_est = zeros(L,length(obj.x_ini));
            P_est = zeros(length(obj.x_ini),length(obj.x_ini),L);
            P = obj.P_ini;
            x = obj.x_ini;
            
            [row,col] = size(y);
            
            if col < row
                y = y';
                N = col;
            end
            
            M = length(obj.x_ini);
            
            Wm_0 = 1 - M/(obj.alpha^2*(M+obj.kappa));
            Wm_i = 1/(2*obj.alpha^2*(M+obj.kappa));
            Wm = [Wm_0,repmat(Wm_i,1,2*M)];
            
            Wc_0 = (2 - obj.alpha^2 + obj.beta) - M/(obj.alpha^2*(M+obj.kappa));
            Wc_i = 1/(2*obj.alpha^2*(M+obj.kappa));
            Wc = [Wc_0,repmat(Wc_i,1,2*M)];
            
            g = zeros(N,2*M+1);
            f = zeros(M,2*M+1);
            
            for n = 2 : length(y)
                
                Psc = obj.alpha^2*(M + obj.kappa)*P;
                dx_i1 = (sqrtm(Psc));
                dx_i2 = -(sqrtm(Psc));
                x_i = [x, x + dx_i1, x + dx_i2];
                
                for i = 1 : 2*M+1
                    g(:,i) = obj.meas_fun(x_i(:,i),n);
                end
                
                ym = sum(Wm.*g,2);
                Py = Wc.*(g - ym)*(g - ym)' + obj.R;
                Pxy = Wc.*(x_i - x)*(g - ym)';
                
                K = Pxy/Py;
                x = x + K*(y(:,n) - ym);
                P = P - K*Py*K';
                
                Psc = obj.alpha^2*(M + obj.kappa)*P;
                dx_i1 = (sqrtm(Psc));
                dx_i2 = -(sqrtm(Psc));
                x_i = [x, x + dx_i1, x + dx_i2];
                for i = 1 : 2*M+1
                    f(:,i) = obj.state_fun(x_i(:,i),e_k,X_tr,n); %sqrt(0.5e-8)*e_k(1:n-1)*(exp(-((x_i(:,i) - X_tr(:,1:n-1)).^2)))';%;
                end
                
                fm = sum(Wm.*f,2);
                P = Wc.*(f - fm)*(f - fm)' + obj.Q;
                
                x_est(n,:) = x;
                P_est(:,:,n) = P;
                
                kk = (exp(-obj.kernel_bw*((x - X_tr(:,1:n-1)).^2)))';
                y_pred(n) = obj.step_size*e_k(1:n-1)*kk;
                e_k(n) = X_tr(n) - y_pred(n);
            end
        end
    end
    
    methods(Static)
        function [e] = mse(x,y)
            
            if  ~iscolumn(x)
                x = x';
            end
            
            if  ~iscolumn(y)
                y = y';
            end
            
            e = mean((x-y).^2);
        end
    end
end


