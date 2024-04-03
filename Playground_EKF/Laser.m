classdef Laser < handle
    %LASER
    % The following class contains properties and methods
    % to solve a laser diode rate equation system, the
    % electric field, the output power, the intensity noise
    % and the frequency noise. The laser is considered with constant input
    % current, defined as static parameter in the model itself.
    
    
    % - Different methods can be used for simualting a laser system
    % (numerical methods to solve differential equations - the rate
    % equations). The theoretical model is taken from "Numerical Simulation of Intensity and Phase Noise
    % From Extracted Parameters for CW DFB Lasers " ,Irshaad Fatadin, Member, IEEE, David Ives, and Martin Wicks, IEEE JOURNAL OF QUANTUM ELECTRONICS, VOL. 42, NO. 9, SEPTEMBER 2006
    
    % This is the ideal simulation for heterodyne coherent detection scheme
    properties
        % PRIMARY PARAMETERS
        % LASER RATE EQUATION STATIC PARAMETERS
        laser_primary_param;
        %         I; % input current, [A]
        %         g; % gain slope [?]
        %         tau_n;  % carrier lifetime [s]
        %         tau_p;  % photon lifetime [s]
        %         eps; % (epsilon) non-linear gain compression [?]
        %         beta; % spontaneous recombination factor [-] - possible range [0,1]?
        %         alpha; % line-width enhancement factor [?]
        %         N0; % carriers at transparency [carrier]
        %         lambda; % optical wavelength [m]
        %         nu;% optical frequency
        
        % DETECTOR STATIC PARAMETERS
        detector_primary_param;
        %         dt; %sampling time [s]
        %         df; %sampling frequency [Hz]
        %         T; % absolute temperature [K]
        %         R_L; % load resistance [Ohm]
        %         P_LO; % local oscillator power [W]
        %         ph_LO; % local oscillator phase
        %         nu_LO; %  frequency local oscillator
        %         lambda_LO;
        %		  Fn; % receiver noise figure [-]
        %         eta_d; % Quantum efficiency [-] - possible range [0,1]?
        %         D; % drift
        
        % SECONDARY PARAMETERS
        % OTHER PARAMETERS CALCULATED ON THE PREVIOUS ONE
        laser_secondary_param;
        %         S_ss; %photons steady state
        %         N_ss; %carriers steady state
        %         P_s; % average signal power
        %         I_th; % threshold current
        
        
        % OTHER PARAMETERS CALCULATED ON THE PREVIOUS ONE
        detector_secondary_param;
        %         Resp; % responsivity of photo-diode
        %         sigma_s; % shot noise variance
        %         sigma_th; % thermal noise variance
        %         R; %measurement noise covariance matrix
        %         SNR ; % thermal and shot noise contribution
        theta; % keep all parameters for convenience
        theta_primary; % keep all primary parameters.
        theta_secondary; % keep all secondary parameters
    end
    
    methods
        function obj = Laser(param)
            % LASER Construct an instance of this class, by passing the
            %parameters "param" containing all the properties described
            %above . If param is missing, it uses a default set of
            %parameters (can lately be changed freely)
            % ----- INPUT ----.
            % param: the primary parameters described above
            
            % SET PRIMARY LASER PARAMETERS
            if isfield(param,'I')
                obj.laser_primary_param.I = param.I;
            else
                obj.laser_primary_param.I = 12.9e-3;
            end
            if isfield(param,'g')
                obj.laser_primary_param.g = param.g;
            else
                obj.laser_primary_param.g = 1.13e4;
            end
            if isfield(param,'tau_n')
                obj.laser_primary_param.tau_n = param.tau_n;
            else
                obj.laser_primary_param.tau_n = 0.33e-9;
            end
            if isfield(param,'tau_p')
                obj.laser_primary_param.tau_p = param.tau_p;
            else
                obj.laser_primary_param.tau_p = 7.15e-12;
            end
            if isfield(param,'eps')
                obj.laser_primary_param.eps = param.eps;
            else
                obj.laser_primary_param.eps = 4.58e-8;
            end
            if isfield(param,'beta')
                obj.laser_primary_param.beta = param.beta;
            else
                obj.laser_primary_param.beta = 3.54e-5;
            end
            if isfield(param,'alpha')
                obj.laser_primary_param.alpha = param.alpha;
            else
                obj.laser_primary_param.alpha = 4.55;
            end
            if isfield(param,'N0')
                obj.laser_primary_param.N0 = param.N0;
            else
                obj.laser_primary_param.N0 = 8.20e6;
            end
            
            
            try
                try
                    obj.laser_primary_param.lambda = param.lambda;
                    obj.laser_primary_param.nu = Const.c/obj.laser_primary_param.lambda;
                catch
                    obj.laser_primary_param.nu = param.nu;
                    obj.laser_primary_param.lambda = Const.c/obj.laser_primary_param.nu;
                end
            catch
                obj.laser_primary_param.lambda = 1550e-9;
                obj.laser_primary_param.nu = Const.c/obj.laser_primary_param.lambda;
            end
            
            
            % DETECTOR
            
            try
                try
                    obj.detector_primary_param.dt = param.dt;
                    obj.detector_primary_param.df = 1/obj.detector_primary_param.dt;
                catch
                    obj.detector_primary_param.df = param.df;
                    obj.detector_primary_param.dt = 1/param.detector_primary_param.df;
                end
            catch
                obj.detector_primary_param.dt = 1e-11; %sampling time [s]
                obj.detector_primary_param.df = 1e11; %sampling frequency [Hz]
            end
            
            
            
            if isfield(param,'T')
                obj.detector_primary_param.T = param.T;
            else
                obj.detector_primary_param.T = 293; % absolute temperature [K]
            end
            
            if isfield(param,'R_L')
                obj.detector_primary_param.R_L = param.R_L;
            else
                obj.detector_primary_param.R_L = 50; % load resistance [Ohm]
            end
            
            if isfield(param,'P_LO')
                obj.detector_primary_param.P_LO = param.P_LO;
            else
                obj.detector_primary_param.P_LO = 10^(7/10)*1e-3; % local oscillator power [W]
            end
            
            if isfield(param,'ph_LO')
                obj.detector_primary_param.ph_LO = param.ph_LO;
            else
                obj.detector_primary_param.ph_LO = 0;
            end
            
            try
                try
                    obj.detector_primary_param.lambda_LO = param.lambda_LO;
                    obj.detector_primary_param.nu_LO = Const.c/obj.detector_primary_param.lambda_LO;
                catch
                    obj.detector_primary_param.nu_LO = param.nu_LO;
                    obj.detector_primary_param.lambda_LO = Const.c/obj.detector_primary_param.nu_LO;
                end
            catch
                obj.detector_primary_param.lambda_LO = 1555e-9;
                obj.detector_primary_param.nu_LO = Const.c/obj.detector_primary_param.lambda_LO;
            end
            
            
            if isfield(param,'Fn')
                obj.detector_primary_param.Fn = param.Fn;
            else
                obj.detector_primary_param.Fn = 1; % receiver noise figure [-]
            end
            
            if isfield(param,'eta_d')
                obj.detector_primary_param.eta_d = param.eta_d;
            else
                obj.detector_primary_param.eta_d = 0.5; % Quantum efficiency [-] - possible range [0,1]?
            end
            
            
            if isfield(param,'D')
                obj.detector_primary_param.D = param.D;
            else
                obj.detector_primary_param.D = 0; % drift
            end
            
            % function for calculating the secondary parameters and the
            % rest
            
            [obj.laser_secondary_param,obj.detector_secondary_param,obj.theta,obj.theta_primary,obj.theta_secondary] = Laser.get_secondary_param(obj.laser_primary_param,obj.detector_primary_param);
            
            
        end
        function obj = set(obj,new_param) % set new parameters for the laser - maybe useless function?
            %SET set the new parameters on the laser object
            % ----- INPUT ----.
            % new_param: new parameters for the laser object
            
            
            try
                obj.laser_primary_param.I = new_param.I;
            catch
            end
            try
                obj.laser_primary_param.g = new_param.g;
            catch
            end
            try
                obj.laser_primary_param.tau_n = new_param.tau_n;
            catch
            end
            try
                obj.laser_primary_param.tau_p = new_param.tau_p;
            catch
            end
            try
                obj.laser_primary_param.eps = new_param.eps;
            catch
            end
            try
                obj.laser_primary_param.beta = new_param.beta;
            catch
            end
            try
                obj.laser_primary_param.alpha = new_param.alpha;
            catch
            end
            try
                obj.laser_primary_param.N0 = new_param.N0;
            catch
            end
            try
                try
                    obj.laser_primary_param.lambda = new_param.lambda;
                    obj.laser_primary_param.nu = Const.c/obj.laser_primary_param.lambda;
                catch
                    obj.laser_primary_param.nu = new_param.nu;
                    obj.laser_primary_param.lambda = Const.c/obj.laser_primary_param.nu;
                end
            catch
            end
            
            
            % DETECTOR
            
            try
                try
                    obj.detector_primary_param.dt = new_param.dt;
                    obj.detector_primary_param.df = 1/obj.detector_primary_param.dt;
                catch
                    obj.detector_primary_param.df = new_param.df;
                    obj.detector_primary_param.dt = 1/obj.detector_primary_param.df;
                end
            catch
            end
            
            
            try
                obj.detector_primary_param.T = new_param.T;
            catch
            end
            
            try
                obj.detector_primary_param.R_L = new_param.R_L;
            catch
            end
            
            try
                obj.detector_primary_param.P_LO = new_param.P_LO;
            catch
            end
            
            try
                obj.detector_primary_param.ph_LO = new_param.ph_LO;
            catch
            end
            
            try
                try
                    obj.detector_primary_param.lambda_LO = new_param.lambda_LO;
                    obj.detector_primary_param.nu_LO = Const.c/obj.detector_primary_param.lambda_LO;
                catch
                    obj.detector_primary_param.nu_LO = new_param.nu_LO;
                    obj.detector_primary_param.lambda_LO = Const.c/obj.detector_primary_param.nu_LO;
                end
            catch
            end
            
            try
                obj.detector_primary_param.Fn = new_param.Fn;
            catch
            end
            
            try
                obj.detector_primary_param.eta_d = new_param.eta_d;
            catch
            end
            
            try
                obj.detector_primary_param.D = new_param.D;
            catch
            end
            
            % function for calculating the secondary parameters
            
            [obj.laser_secondary_param,obj.detector_secondary_param,obj.theta,obj.theta_primary,obj.theta_secondary] = Laser.get_secondary_param(obj.laser_primary_param,obj.detector_primary_param);
            
            
            
        end %
        function[f,RIN,FN] = analytic_spectrums(obj,dt,Nt) % (from small signal analysis)
            % ANALYTIC_SPECTRUMS return the analytic spectrums from small signal analysis of the rate equations. This theoretical spectrums are then an approximation of the spectra calculated on the generated internal states
            % ----- INPUT ----.
            % dt: sampling time
            % Nt: number of points in the spectra
            % ----- OUTPUT -----
            % f: frequency vector
            % RIN: relative intensity noise
            % FN: frequency noise
            
            fs = 1/dt; % sampling frequency
            f = (0:Nt-1).'/Nt*fs; %frequency vector
            omega = 2*pi*f; % angular frequency
            DNN = 2*obj.theta.N_ss*(1+obj.theta.beta*obj.theta.S_ss)/obj.theta.tau_n;
            DSS = 2*obj.theta.beta*obj.theta.N_ss*obj.theta.S_ss/obj.theta.tau_n;
            DSN = -DSS;
            DPP = obj.theta.beta*obj.theta.N_ss/(2*obj.theta.tau_n*obj.theta.S_ss);
            Aw = @(w)[1i.*w + 1/obj.theta.tau_n + obj.laser_primary_param.g*obj.theta.S_ss/(1+obj.theta.eps*obj.theta.S_ss),...
                obj.theta.g*(obj.theta.N_ss - obj.theta.N0)/((1+obj.theta.eps*obj.theta.S_ss)^2);...
                -obj.theta.g*obj.theta.S_ss/(1+obj.theta.eps*obj.theta.S_ss)- obj.theta.beta/obj.theta.tau_n, ...
                1i.*w + 1/obj.theta.tau_p - obj.theta.g*(obj.theta.N_ss - obj.theta.N0)/((1 + obj.theta.eps*obj.theta.S_ss)^2)];
            FN = zeros(Nt,1);
            RIN = zeros(Nt,1);
            for k = 1:length(omega)
                wk = omega(k);
                Bwk = inv(Aw(wk));
                FN(k) = 1/(4*pi^2)*((obj.theta.alpha*obj.theta.g/2)^2*(Bwk(1,1)*conj(Bwk(1,1))*DNN + Bwk(1,2)*conj(Bwk(1,2))*DSS + DSN*(Bwk(1,1)*conj(Bwk(1,2)) + Bwk(1,2)*conj(Bwk(1,1)))) + DPP);
                RIN(k) = 10*log10((Bwk(2,1)*conj(Bwk(2,1))*DNN + Bwk(2,2)*conj(Bwk(2,2))*DSS + DSN*(Bwk(2,1)*conj(Bwk(2,2)) + Bwk(2,2)*conj(Bwk(2,1))))/(obj.theta.S_ss^2));
            end
        end
        
        function [X,Y,Yc,P] = simulating(obj,Nt,X0,solver,q)
            %SIMULATING a laser system according to a process and a measurement model,
            %using Euler-Mayruama or S-ROCK sotchastic Chebyschev method
            % X are the laser's states while Y are the measurements
            % proc_model and meas_model are structures containing:
            % - F(x,theta): deterministic transition function for one step approximation
            % - G(x,theta,q): stochastic transition function for one step
            % approximation, theta are the static parameters, q are the
            % noise components
            % - Q/R: covariance matrices of process/measurement noise
            % - Nx/Ny: number of states/measurements
            % static parameters of the model should be already included in
            % the current obj
            %
            % ----- INPUT ----.
            % Nt: length of sequence to be created
            % X0: initial point
            % solver: type of solver implemented
            %   - type: can be Euler-Mayurama or S-ROCK . For S-ROCK, read Abdulle and Yi : "S-Rock Methods For Stiff Ito Sdes"
            %   - m / eta : parameters for S-ROCK solver
            % q: optional - vector of Nt*3 of white Gaussian samples. Useful to compare different solver on the same sequence
            % ----- OUTPUT -----
            % X: laser states
            % Y: detected photocurrent with measurement noise
            % Yc: detected photocurrent wthout measurement noise
            % P: optical power output
            
            
            % we assume the model to be fixed
            Nx = 3;
            Nr = 2;
            Nq = 3;
            if nargin < 5 % in case we supply with our random vector
                q = randn(Nt,Nq);
            end
            X = zeros(Nt,Nx); % states
            Y = zeros(Nt,2); % measures
            Yc = zeros(Nt,2); % clean measures
            P = zeros(Nt,1); % optical power (from the laser source only, no local oscillator)
            
            X(1,:) = X0; %starting state
            r = randn(Nt,Nr).*repmat(sqrt(diag(obj.theta.R)).',Nt,1);
            switch solver.type
                case 'EM'
                    for k = 1:Nt
                        % states update
                        X(k+1,:) = X(k,:) +  obj.theta.dt.*obj.drift(X(k,:)) + sqrt(obj.theta.dt).*(obj.diffusion(X(k,:))*(q(k,:).')).';
                        if sum(X(k+1,1:2) < 0) > 0
                            disp('Convergence problems. Suggested a time-step reduction')
                        end
                    end
                case 'MILSTEIN'
                    error('Not implemented')
                case 'HEUN'
                    error('Not implemented')
                case 'RK4'
                    error('Not implemented')
                case 'CUSTOM'
                    error('Not implemented')
                case  'SROCK'
                    if ~isfield(solver,'m')
                        solver.m = 10;
                    end
                    if ~isfield(solver,'eta')
                        solver.eta = 100;
                    end
                    
                    omega0 = 1+solver.eta/solver.m^2;
                    Cp = zeros(solver.m,1);
                    for im = 1:solver.m
                        Cp(im) = chebyshevT(im,omega0);
                    end
                    syms x
                    Cp_prime_m_s =diff(chebyshevT(im,x),x);
                    x = omega0;
                    Cp_prime_m = double(subs(Cp_prime_m_s));
                    omega1 = Cp(im)/Cp_prime_m;
                    for k = 1:Nt
                        % states update
                        % deterministic Chebyschev stages
                        K = zeros(solver.m,Nx);
                        K(1,:) = X(k,:);
                        K(2,:) = X(k,:) + 2*obj.theta.dt*omega1/omega0.*obj.drift(K(1,:));
                        for j = 3:solver.m
                            K(j,:) = 2*obj.theta.dt*omega1*(Cp(j-1)/Cp(j))*obj.drift(K(j-1,:)) + 2*omega0*(Cp(j-1)/Cp(j))*K(j-1,:) - (Cp(j-2)/Cp(j))*K(j-2,:);
                        end
                        % stochastic finishing stage
                        X(k+1,:) = K(solver.m,:) + (sqrt(obj.theta.dt).*obj.diffusion(K(solver.m-1,:))*(q(k,:).')).';
                        if sum(X(k+1,1:2) < 0) > 0
                            disp('Convergence problems. Suggested a time-step reduction')
                        end
                    end
                    
                otherwise
                    error('Incorrect solver method.')
                    
            end
            X = X(2:end,:); % reshape the state sequence
            for k = 1:Nt
                % measurements update
                P(k) = X(k,2)*Const.h*obj.theta.nu/obj.theta.tau_p;
                Yc(k,:) = 2*(obj.theta.eta_d*Const.q/(Const.h*obj.theta.nu))*sqrt(P(k)*obj.theta.P_LO).*[cos(2*pi*(obj.theta.nu - obj.theta.nu_LO)*k*obj.theta.dt + obj.theta.D*(k*obj.theta.dt)^2 + (X(k,3) - obj.theta.ph_LO)),sin(2*pi*(obj.theta.nu - obj.theta.nu_LO)*k*obj.theta.dt + obj.theta.D*(k*obj.theta.dt)^2 + (X(k,3) - obj.theta.ph_LO))];
                Y(k,:) = Yc(k,:) +[r(k,1),r(k,2)];
            end
            
        end
        function f_x = drift(obj,x)
            % DRIFT the drift function (deterministic part) of the rate equations
            % ----- INPUT ----.
            % x: current state
            % obj: laser object
            % ----- OUTPUT -----
            % f_x: drift value for state x
            f_x = [obj.theta.I/Const.q - x(1)/obj.theta.tau_n - obj.theta.g*(x(1) - obj.theta.N0)*x(2)/(1+obj.theta.eps*x(2)),...
                obj.theta.g*(x(1) - obj.theta.N0)*x(2)/(1+obj.theta.eps*x(2)) - x(2)/obj.theta.tau_p + obj.theta.beta*x(1)/obj.theta.tau_n ,...
                obj.theta.alpha*obj.theta.g/2*(x(1) - obj.theta.N_ss) ];
        end
        function f_x = f_drift_N(obj,x)
            % F_DRIFT_N the drift function (deterministic part) of the rate equations, only the carrier drift.
            % ----- INPUT ----.
            % x: current state
            % obj: laser object
            % ----- OUTPUT -----
            % f_x: drift value for state x
            f_x = obj.theta.I/Const.q - x(1)/obj.theta.tau_n - obj.theta.g*(x(1) - obj.theta.N0)*x(2)/(1+obj.theta.eps*x(2));
        end
        function f_x = f_drift_S(obj,x)
            % F_DRIFT_S the drift function (deterministic part) of the rate equations, only the photons drift.
            % ----- INPUT ----.
            % x: current state
            % obj: laser object
            % ----- OUTPUT -----
            % f_x: drift value for state x
            f_x = obj.theta.g*(x(1) - obj.theta.N0)*x(2)/(1+obj.theta.eps*x(2)) - x(2)/obj.theta.tau_p + obj.theta.beta*x(1)/obj.theta.tau_n;
        end
        function Jf_x = jacobian_drift(obj,x)
            % JACOBIAN_DRIFT the  jacobian of the drift function (deterministic part) of the rate equations. Useful for stability checks
            % ----- INPUT ----.
            % x: current state
            % obj: laser object
            % ----- OUTPUT -----
            % Jf_x: jacobian drift value for state x
            Jf_x = [(-1/obj.theta.tau_n-obj.theta.g*x(2)/(1+ obj.theta.eps*x(2))), (- obj.theta.g*(x(1) - obj.theta.N0)/(1+ obj.theta.eps*x(2))^2 ),0;
                ( obj.theta.beta/obj.theta.tau_n + obj.theta.g*x(2)/(1+ obj.theta.eps*x(2))),(obj.theta.g*(x(1) - obj.theta.N0)/(1+ obj.theta.eps*x(2))^2 - 1/obj.theta.tau_p),0;
                obj.theta.alpha/2*obj.theta.g,0,0];
        end
        function Qc = diffusion(obj,x)
            % DIFFUSION The diffusion function of the rate equation (stochastic part). From the way the langevin noise is generated, this corresponds to the cholesky lower triangular decomposition of langevin covariance noise. Qc*Qc.' =Q
            % in this way the diffusion process can be generated by a
            % linear function, Gamma = Qc*q where Gamma are the Langevin
            % noises, q is multivariate normal distributed vector
            % ----- INPUT ----.
            % x: current state
            % obj: laser object
            % ----- OUTPUT -----
            % Q_c: diffusion matrix
            
            Qc = [sqrt(2*x(1)/obj.theta.tau_n*(1+obj.theta.beta*x(2))), 0,0;
                - sqrt(2*obj.theta.beta^2*x(1)*x(2)^2/(obj.theta.tau_n*(1+obj.theta.beta*x(2)))),  sqrt(2*obj.theta.beta*x(1)*x(2)/(obj.theta.tau_n*(1+obj.theta.beta*x(2)))), 0;
                0,0, sqrt(obj.theta.beta*x(1)/(2*obj.theta.tau_n*x(2)))];
        end
        function Qc = diffusion2(obj,x)
            % DIFFUSION2 The diffusion function of the rate equation (stochastic part). From the way the langevin noise is generated, this corresponds to the cholesky upper triangular decomposition of langevin covariance noise. Qc.'*Qc =Q
            % in this way the diffusion process can be generated by a
            % linear function, Gamma = q.'*Qc where Gamma are the Langevin
            % noises, q is multivariate normal distributed vector
            % ----- INPUT ----.
            % x: current state
            % obj: laser object
            % ----- OUTPUT -----
            % Q_c: diffusion matrix
            
            Qc = [sqrt(2*x(1)/obj.theta.tau_n*(1+obj.theta.beta*x(2))), - sqrt(2*obj.theta.beta^2*x(1)*x(2)^2/(obj.theta.tau_n*(1+obj.theta.beta*x(2)))),0;
                0,  sqrt(2*obj.theta.beta*x(1)*x(2)/(obj.theta.tau_n*(1+obj.theta.beta*x(2)))), 0;
                0,0, sqrt(obj.theta.beta*x(1)/(2*obj.theta.tau_n*x(2)))];
        end
        
        
        function Q = langevin_cov(obj,x) % langevin covariance noise
            % LANGEVIN_COV The langevin covariance matrix of the noise sources
            % ----- INPUT ----.
            % x: current state
            % obj: laser object
            % ----- OUTPUT -----
            % Q: covariance matrix
            Q = [(2*x(1)/obj.theta.tau_n*(1+obj.theta.beta*x(2))), -2*obj.theta.beta*x(1)*x(2)/obj.theta.tau_n,0;
                -2*obj.theta.beta*x(1)*x(2)/obj.theta.tau_n, 2*obj.theta.beta*x(1)*x(2)/obj.theta.tau_n, 0;
                0,0, (obj.theta.beta*x(1)/(2*obj.theta.tau_n*x(2)))];
        end
    end
    methods(Static)
        function Q = langevin_covariance(x,theta)
            % LANGEVIN_COVARIANCE The langevin covariance matrix of the noise sources
            % ----- INPUT ----.
            % x: current state
            % theta: rate equation parameters
            % ----- OUTPUT -----
            % Q: covariance matrix
            Q = [(x(1)/theta.tau_n*(1+theta.beta*x(2))), -theta.beta*x(1)*x(2)/theta.tau_n;
                -theta.beta*x(1)*x(2)/theta.tau_n, theta.beta*x(1)*x(2)/theta.tau_n];
        end
        function dX = get_derivatives(x,theta)
            % GET_DERIVATIVES Calculate the deterministic derivatives (drift) of the rate equations (without noise), for carriers and photons (no phase)
            % ----- INPUT ----.
            % x: current state
            % theta: rate equation parameters
            % ----- OUTPUT -----
            % dX: derivative (drift)
            dX = [theta.I/Const.q - x(1)/theta.tau_n - theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)),...
                theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)) - x(2)/theta.tau_p + theta.beta*x(1)/theta.tau_n ];
        end
        function [f,Fx,Fq] = predict(x,q,theta,t,param)
            % PREDICT Prediciton function for the Kalman filter setup, specifically it implements the rate equations as part of the Kalman (extended/unscented) prediciton step
            % ----- INPUT ----.
            % x: current state
            % q: process noise state
            % theta: rate equation parameters
            % t: time as index
            % param: parameters for the solver method
            %   - solver: can be Euler-Mayurama or S-ROCK . For S-ROCK, read Abdulle and Yi : "S-Rock Methods For Stiff Ito Sdes"
            %   - m / eta : parameters for S-ROCK solver
            % ----- OUTPUT -----
            % f: the prediction of the next point according to the rate equations
            % Fx: Jacobian of prediction wrt the states
            % Fq: Jacobian of prediction wrt the noise
            switch param.solver
                case 'EM'
                    f = x + [theta.I/Const.q - x(1)/theta.tau_n - theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)),...
                        theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)) - x(2)/theta.tau_p + theta.beta*x(1)/theta.tau_n,...
                        theta.alpha*theta.g/2*(x(1) - theta.N_ss)].*theta.dt + ... drift process
                        ...
                        ([sqrt(2*x(1)/theta.tau_n*(1+theta.beta*x(2))), 0,0;
                        - sqrt(2*theta.beta^2*x(1)*x(2)^2/(theta.tau_n*(1+theta.beta*x(2)))),  sqrt(2*theta.beta*x(1)*x(2)/(theta.tau_n*(1+theta.beta*x(2)))), 0;
                        0,0, sqrt(theta.beta*x(1)/(2*theta.tau_n*x(2)))].*sqrt(theta.dt)*(q.')).'; % diffusion process
                    f = [f(1:2).*(f(1:2) > 0),f(3)]; % if the first 2 components are below zero, we will just modify it forcing them to 0
                    if nargout > 1
                        Fx =  [1 + (-1/theta.tau_n-theta.g*x(2)/(1+ theta.eps*x(2)))*theta.dt +            sqrt((1+theta.beta*x(2))/(2*x(1)*theta.tau_n))*sqrt(theta.dt)*q(1),... dfN/dN ok
                            (- theta.g*(x(1) - theta.N0)/(1+ theta.eps*x(2))^2 )*theta.dt + theta.beta*sqrt(x(1)/((1 + theta.beta*x(2))*2*theta.tau_n))*sqrt(theta.dt)*q(1)  ,... dfN/dS ok
                            0; ... dfN/dphi ok
                            ( theta.beta/theta.tau_n + theta.g*x(2)/(1+ theta.eps*x(2)))*theta.dt +      sqrt(theta.beta*x(2)/(2*x(1)*theta.tau_n*(1+theta.beta*x(2))))*sqrt(theta.dt)*q(2) - theta.beta*x(2)/sqrt(x(1)*2*theta.tau_n*(1+theta.beta*x(2)))*sqrt(theta.dt)*q(1) ,... dfS/dN
                            1+ (theta.g*(x(1) - theta.N0)/(1+ theta.eps*x(2))^2 - 1/theta.tau_p)*theta.dt + sqrt(theta.beta*x(1)/(2*x(2)*theta.tau_n*(1+theta.beta*x(2))^3))*sqrt(theta.dt)*q(2) - (2*x(2) + theta.beta*x(2)^2)*sqrt(theta.beta^2*x(1)/(2*theta.tau_n*x(2)^2*(1+theta.beta*x(2))^3))*sqrt(theta.dt)*q(1), ... dfS/dS
                            0 ; ... dfS/dphi
                            theta.alpha/2*theta.g*theta.dt + sqrt(theta.beta/(8*x(1)*x(2)*theta.tau_n))*sqrt(theta.dt)*q(3),... dfphi/dN
                            -sqrt(theta.beta*x(1)/(8*x(2)^3*theta.tau_p))*sqrt(theta.dt)*q(3)       ,...dfphi/dS
                            1]; % dfphi/dphi
                        
                        Fq =  [sqrt(2*x(1)/theta.tau_n*(1+theta.beta*x(2))), 0,0;
                            - sqrt(2*theta.beta^2*x(1)*x(2)^2/(theta.tau_n*(1+theta.beta*x(2)))),  sqrt(2*theta.beta*x(1)*x(2)/(theta.tau_n*(1+theta.beta*x(2)))), 0;
                            0,0, sqrt(theta.beta*x(1)/(2*theta.tau_n*x(2)))].*sqrt(theta.dt);
                    end
                case 'SROCK'
                    if ~isfield(param,'m')
                        param.m = 10;
                    end
                    if ~isfield(param,'eta')
                        param.eta = 100;
                    end
                    if nargout == 1
                        dr = @(x)[theta.I/Const.q - x(1)/theta.tau_n - theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)),...
                            theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)) - x(2)/theta.tau_p + theta.beta*x(1)/theta.tau_n,...
                            theta.alpha*theta.g/2*(x(1) - theta.N_ss)];
                        
                        di = @(x)[sqrt(2*x(1)/theta.tau_n*(1+theta.beta*x(2))), 0,0;
                            - sqrt(2*theta.beta^2*x(1)*x(2)^2/(theta.tau_n*(1+theta.beta*x(2)))),  sqrt(2*theta.beta*x(1)*x(2)/(theta.tau_n*(1+theta.beta*x(2)))), 0;
                            0,0, sqrt(theta.beta*x(1)/(2*theta.tau_n*x(2)))];
                        
                        k = zeros(param.m,3);
                        k(1,:) = x;
                        k(2,:) = x + 2*theta.dt*param.omega1/param.omega0.*dr(k(1,:));
                        for j = 3:param.m
                            k(j,:) = 2*theta.dt*param.omega1*(param.Cp(j-1)/param.Cp(j))*dr(k(j-1,:)) + 2*param.omega0*(param.Cp(j-1)/param.Cp(j))*k(j-1,:) - (param.Cp(j-2)/param.Cp(j))*k(j-2,:);
                        end
                        f = k(param.m,:) +  (sqrt(theta.dt).*di(k(param.m-1,:))*(q.')).';
                        
                    else
                        dr = @(x)[theta.I/Const.q - x(1)/theta.tau_n - theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)),...
                            theta.g*(x(1) - theta.N0)*x(2)/(1+theta.eps*x(2)) - x(2)/theta.tau_p + theta.beta*x(1)/theta.tau_n,...
                            theta.alpha*theta.g/2*(x(1) - theta.N_ss)];
                        
                        di = @(x)[sqrt(2*x(1)/theta.tau_n*(1+theta.beta*x(2))), 0,0;
                            - sqrt(2*theta.beta^2*x(1)*x(2)^2/(theta.tau_n*(1+theta.beta*x(2)))),  sqrt(2*theta.beta*x(1)*x(2)/(theta.tau_n*(1+theta.beta*x(2)))), 0;
                            0,0, sqrt(theta.beta*x(1)/(2*theta.tau_n*x(2)))];
                        
                        Jdr_x = @(x)[  -1/theta.tau_n-theta.g*x(2)/(1+ theta.eps*x(2))    , - theta.g*(x(1) - theta.N0)/(1+ theta.eps*x(2))^2         ,    0;
                            theta.beta/theta.tau_n + theta.g*x(2)/(1+ theta.eps*x(2)), theta.g*(x(1) - theta.N0)/(1+ theta.eps*x(2))^2 - 1/theta.tau_p , 0  ;
                            theta.alpha/2*theta.g  ,       0       ,   0];
                        
                        Jdi_Nx =  @(x)[sqrt((1+theta.beta*x(2))/(2*x(1)*theta.tau_n)),theta.beta*sqrt(x(1)/((1 + theta.beta*x(2))*2*theta.tau_n)),0;
                            - theta.beta*x(2)/sqrt(x(1)*2*theta.tau_n*(1+theta.beta*x(2))), - (2*x(2) + theta.beta*x(2)^2)*sqrt(theta.beta^2*x(1)/(2*theta.tau_n*x(2)^2*(1+theta.beta*x(2))^3)),0;
                            0,0,0];
                        
                        Jdi_Sx =  @(x)[0,0,0;
                            sqrt(theta.beta*x(2)/(2*x(1)*theta.tau_n*(1+theta.beta*x(2)))),sqrt(theta.beta*x(1)/(2*x(2)*theta.tau_n*(1+theta.beta*x(2))^3)),0
                            0,0,0];
                        
                        Jdi_phix =  @(x)[0,0,0;
                            0,0,0;
                            sqrt(theta.beta/(8*x(1)*x(2)*theta.tau_n)),-sqrt(theta.beta*x(1)/(8*x(2)^3*theta.tau_p)),0];
                        
                        k = zeros(param.m,3);
                        K = zeros(3,3,param.m);
                        k(1,:) = x;
                        K(:,:,1) = eye(3);
                        k(2,:) = x + 2*theta.dt*param.omega1/param.omega0.*dr(k(1,:));
                        K(:,:,2) = eye(3) + 2*theta.dt*param.omega1/param.omega0.*Jdr_x(k(1,:));
                        for j = 3:param.m
                            k(j,:) = 2*theta.dt*param.omega1*(param.Cp(j-1)/param.Cp(j))*dr(k(j-1,:))                 + 2*param.omega0*(param.Cp(j-1)/param.Cp(j))*k(j-1,:)              - (param.Cp(j-2)/param.Cp(j))*k(j-2,:);
                            K(:,:,j) = 2*theta.dt*param.omega1*(param.Cp(j-1)/param.Cp(j))*Jdr_x(k(j-1,:))*K(:,:,j-1) + 2*param.omega0*(param.Cp(j-1)/param.Cp(j))*K(:,:,j-1)             - (param.Cp(j-2)/param.Cp(j))*K(:,:,j-2);
                        end
                        % stochastic finishing stage
                        f = k(param.m,:) +  (sqrt(theta.dt).*di(k(param.m-1,:))*(q.')).';
                        Fx = K(:,:,param.m) + (Jdi_Nx(k(param.m-1,:)).*q(1) + Jdi_Sx(k(param.m-1,:)).*q(2) + Jdi_phix(k(param.m-1,:)).*q(3))*K(:,:,param.m-1).*sqrt(theta.dt) ;
                        Fq = sqrt(theta.dt).*di(k(param.m-1,:)) ;
                    end
                    f = [f(1:2).*(f(1:2) > 0),f(3)]; % if the first 2 components are below zero, we will just modify it forcing them to 0
                otherwise
                    error('Numerical solver not implemented yet')
            end
        end
        function [h,Hx,Hr] = measure(x,r,theta,t,param)
            % MEASURE Measurement function for the Kalman filter setup, specifically it implements the heterodyne measurement equation as part of the Kalman (extended/unscented) update step
            % it can be used for both quadratures or for a single quadrature
            % ----- INPUT ----.
            % x: current state
            % r: measurement noise state
            % theta: heterodyne receiver parameters
            % t: time as index
            % param: parameters for the equation
            %   - RE: rate equation based
            %   - RW: random walk based
            % ----- OUTPUT -----
            % h: the measured photocurrent given the laser state x
            % Hx: Jacobian of the measured photocurrent wrt the states
            % Hr: Jacobian of the measured photocurrent of prediction wrt the measurement noise
            Ny = length(r);
            switch Ny
                case 1
                    switch param.model
                        case 'RE' % rate equation based
                            h =  2*(theta.eta_d*Const.q/(Const.h*theta.nu))*sqrt(x(2)*Const.h*theta.nu/theta.tau_p*theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO)+r(1);
                            if nargout > 1
                                Hx = [0, (theta.eta_d*Const.q/(Const.h*theta.nu))/(sqrt(x(2)))*sqrt(Const.h*theta.nu/theta.tau_p*theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO), - 2*(theta.eta_d*Const.q/(Const.h*theta.nu))*sqrt(x(2)*Const.h*theta.nu/theta.tau_p*theta.P_LO)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO)];
                                Hr = 1;
                            end
                        case 'RW' % amplitude-phase based: x(1) amplitude, x(2) phase
                            h =  2*(theta.eta_d*Const.q/(Const.h*theta.nu))*abs(x(1))*sqrt(theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))+r(1);
                            if nargout > 1
                                Hx = [2*sign(x(1))*(theta.eta_d*Const.q/(Const.h*theta.nu))*sqrt(theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)), - 2*(theta.eta_d*Const.q/(Const.h*theta.nu))*x(1)*sqrt(theta.P_LO)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))];
                                Hr = 1;
                            end
                    end
                case 2
                    switch param.model
                        case 'RE' % rate equation based
                            h =  2*(theta.eta_d*Const.q/(Const.h*theta.nu))*sqrt(x(2)*Const.h*theta.nu/theta.tau_p*theta.P_LO)*[cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO),sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO)]+[r(1),r(2)];
                            if nargout > 1
                                Hx = [  0, (theta.eta_d*Const.q/(Const.h*theta.nu))/(sqrt(x(2)))*sqrt(Const.h*theta.nu/theta.tau_p*theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO), - 2*(theta.eta_d*Const.q/(Const.h*theta.nu))*sqrt(x(2)*Const.h*theta.nu/theta.tau_p*theta.P_LO)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO);
                                    0, (theta.eta_d*Const.q/(Const.h*theta.nu))/(sqrt(x(2)))*sqrt(Const.h*theta.nu/theta.tau_p*theta.P_LO)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO),  2*(theta.eta_d*Const.q/(Const.h*theta.nu))*sqrt(x(2)*Const.h*theta.nu/theta.tau_p*theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) - theta.ph_LO)];
                                Hr = eye(2);
                            end
                        case 'RW'  % amplitude-phase based: x(1) amplitude, x(2) phase
                            h =  2*(theta.eta_d*Const.q/(Const.h*theta.nu))*abs(x(1))*sqrt(theta.P_LO)*[cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)),sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))]+[r(1),r(2)];
                            if nargout > 1
                                Hx = [  2*sign(x(1))*theta.eta_d*Const.q/(Const.h*theta.nu)*sqrt(theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)), - 2*theta.eta_d*Const.q/(Const.h*theta.nu)*x(1)*sqrt(theta.P_LO)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2));
                                    2*sign(x(1))*theta.eta_d*Const.q/(Const.h*theta.nu)*sqrt(theta.P_LO)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)),  2*theta.eta_d*Const.q/(Const.h*theta.nu)*x(1)*sqrt(theta.P_LO)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))];
                                Hr = eye(2);
                            end
                    end
                otherwise
                    error('Model not yet implemented.')
            end
        end
        function [h,Hx,Hr] = measure_exp_data(x,r,theta,t,param) % this is used for tracking experimental data
            % MEASURE_EXP_DATA similar to MEASURe (above), but it is designed for experimental data, where most of the constants in the equations are marged together
            % ----- INPUT ----.
            % x: current state
            % r: measurement noise state
            % theta: heterodyne receiver parameters
            % t: time as index
            % param: parameters for the equation
            %   - RE: rate equation based
            %   - RW: random walk based
            % ----- OUTPUT -----
            % h: the measured photocurrent given the laser state x
            % Hx: Jacobian of the measured photocurrent wrt the states
            % Hr: Jacobian of the measured photocurrent of prediction wrt the measurement noise
            Ny = length(r);
            switch Ny
                case 1
                    switch param.model
                        case 'RE' % rate equation based
                            h =  theta.K*sqrt(x(2)/theta.tau_p)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3))+r(1);
                            if nargout > 1
                                Hx = [0, theta.K/sqrt(theta.tau_p*x(2))*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3)),- theta.K*sqrt(x(2)/theta.tau_p)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3))];
                                Hr = 1;
                            end
                        case 'RW' % amplitude-phase based: x(1) amplitude, x(2) phase
                            h =  sqrt(x(1))*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))+r(1);
                            if nargout > 1
                                Hx = [1/sqrt(x(1))*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)), - sqrt(x(1))*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))];
                                Hr = 1;
                            end
                    end
                case 2
                    switch param.model
                        case 'RE' % rate equation based
                            h =  theta.K*sqrt(x(2)/theta.tau_p)*[cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3)),sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3) )]+[r(1),r(2)];
                            if nargout > 1
                                Hx = [  0, theta.K/sqrt(theta.tau_p*x(2))*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3)),- theta.K*sqrt(x(2)/theta.tau_p)*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3));
                                    0, theta.K/sqrt(theta.tau_p*x(2))*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3)), theta.K*sqrt(x(2)/theta.tau_p)*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(3))];
                                Hr = eye(2);
                            end
                        case 'RW'  % amplitude-phase based: x(1) amplitude, x(2) phase
                            h =  sqrt(x(1))*[cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)),sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))]+[r(1),r(2)];
                            if nargout > 1
                                Hx = [ 1/sqrt(x(1))*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)), - sqrt(x(1))*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2));
                                    1/sqrt(x(1))*sin(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2)),   sqrt(x(1))*cos(theta.omega*(t*theta.dt) + theta.D*(t*theta.dt)^2 + x(2))];
                                Hr = eye(2);
                            end
                    end
                otherwise
                    error('Model not yet implemented.')
            end
        end
        function [RIN,f] = RIN(P,dt,W,s)
            % RIN relative intensity noise spectrum, from the optical power P and the sampling time dt
            % ----- INPUT ----.
            % P: power time series of the signal we want to calculate the RIN
            % dt: sampling time of P
            % W: (optional) window size. If no specified, the fft is calculated on the whole length of P. If specified, the fft is calculated in a moving window W that moves s steps
            % s: (optional) step-length for the moving window (above). If not specified, s = 1 (this can be quite slow)
            % ----- OUTPUT -----
            % RIN: relative intensity noise spectrum
            % f: frequency vector that matches the RIN dimensions
            Nt = length(P);
            fs = 1/dt; % sampling frequency
            
            if nargin < 3
                f = (0:Nt-1).'/Nt*fs; %frequency vector
                RIN = ((1/(mean(P).^2*Nt*dt)).*abs(fft(P - mean(P)).*dt).^2);
            else
                if W > Nt
                    W = Nt;
                end
                K = floor((Nt-W)/s)+1;
                f = (0:W-1).'/W*fs; %frequency vector
                RIN = zeros(W,1);
                for k = 1:K
                    Pk = P((k-1)*s + (1:W));
                    RIN = (k-1)/k.*RIN + ((1/(mean(Pk).^2*W*dt)).*abs(fft(Pk - mean(Pk)).*dt).^2)./k;
%                     disp(['RIN calculation: ',num2str(k/K*100),' % completed.'])
                end
            end
            % return the single side density
            RIN = RIN(1:floor(Nt/2)+1,:);
            %             f_ss = (0:(1/Ts)/N:(1/Ts)/2).';
            f = (0:Nt-1).'/Nt/dt;
            f = f(f <= 1/2/dt);
            RIN(2:end-1,:) = 2*RIN(2:end-1,:);
            
        end
        function [FN,f] = FN(phi,dt,W,s)
            % FN frequency noise spectrum, from the optical phase phi and the sampling time dt
            % ----- INPUT ----.
            % phi: phase time series of the signal we want to calculate the FN
            % dt: sampling time of phi
            % W: (optional) window size. If no specified, the fft is calculated on the whole length of phi. If specified, the fft is calculated in a moving window W that moves s steps
            % s: (optional) step-length for the moving window (above). If not specified, s = 1 (this can be quite slow)
            % ----- OUTPUT -----
            % FN: frequency noise spectrum
            % f: frequency vector that matches the FN dimensions
            fs = 1/dt; % sampling frequency
            d_ni = diff(phi)./(2*pi*dt);
            Nt = length(d_ni);
            
            if nargin < 3
                f = (0:Nt-1).'/Nt*fs; %frequency vector
                FN = (1/(Nt*dt)).*abs(fft(d_ni)*dt).^2;
            else
                if W > Nt
                    W = Nt;
                end
                K = floor((Nt-W)/s)+1;
                f = (0:W-1).'/W*fs; %frequency vector
                FN = zeros(W,1);
                for k = 1:K
                    d_ni_k = d_ni((k-1)*s + (1:W));
                    FN = (k-1)/k.*FN +  ((1/(W*dt)).*abs(fft(d_ni_k)*dt).^2)./k;
%                     disp(['FN calculation: ',num2str(k/K*100),' % completed.'])
                end
            end
            % return the single side density
            FN = FN(1:floor(Nt/2)+1,:);
            %             f_ss = (0:(1/Ts)/N:(1/Ts)/2).';
            f = (0:Nt-1).'/Nt/dt;
            f = f(f <= 1/2/dt);
            FN(2:end-1,:) = 2*FN(2:end-1,:);
            
        end
        function [omega0,omega1,Cp] = get_srock_param(m,eta)
            % GET_SROCK_PARAM get parameters for the srock method. This function is normally called once before the filtering, since the input paramters are fixed.
            % ----- INPUT ----.
            % m,eta: s-rock method input parameters
            % ----- OUTPUT -----
            % omega0, omega1,Cp : s-rock method output parameters
            omega0 = 1+eta/m^2;
            Cp = zeros(m,1);
            for i = 1:m
                Cp(i) = chebyshevT(i,omega0);
            end
            syms x
            Cp_prime_m_s =diff(chebyshevT(i,x),x);
            x = omega0;
            Cp_prime_m = double(subs(Cp_prime_m_s));
            omega1 = Cp(i)/Cp_prime_m;
            clear x;
        end
        function [laser_secondary_param,detector_secondary_param,theta,theta_primary, theta_secondary] = get_secondary_param(l_param,d_param)
            % GET_SECONDARY_PARAM get secondary parameters for the laser-detector model. Those parameters are dependent from the primary parameters
            % ----- INPUT ----.
            % l_param: laser primary parameters
            % d_param: detector primary parameters
            % ----- OUTPUT -----
            % laser_secondary_param: laser secondary parameters
            % detector_secondary_param: detector secondary parameters
            % theta: all parameters updates
            % theta_primary: all primary parameters
            % theta_secondary: all secondary parameters
            
            % STEADY STATES
            % Laser secondary parameters
            aa = Const.q*(l_param.eps + l_param.g*l_param.tau_n);
            bb = Const.q + l_param.g*l_param.tau_p*((1-l_param.beta)*l_param.N0*Const.q-l_param.tau_n*l_param.I) - l_param.beta*l_param.eps*l_param.tau_p*l_param.I;
            cc = -l_param.beta*l_param.tau_p*l_param.I;
            laser_secondary_param.S_ss=(-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa); % ok
            laser_secondary_param.N_ss= l_param.tau_n/Const.q*(l_param.I*(1 + l_param.eps*laser_secondary_param.S_ss) + l_param.g*l_param.N0*Const.q*laser_secondary_param.S_ss)/(1 + l_param.eps*laser_secondary_param.S_ss + l_param.g*laser_secondary_param.S_ss*l_param.tau_n); % ok
            laser_secondary_param.Ps = Const.h*Const.c/(l_param.lambda*l_param.tau_p)*laser_secondary_param.S_ss;
            laser_secondary_param.I_th = Const.q/l_param.tau_n*(l_param.N0 + 1/(l_param.g*l_param.tau_p));
            
            % detector secondary parameters
            % POWER AND SNR
            detector_secondary_param.sigma_s = 2*Const.q*d_param.eta_d*(laser_secondary_param.Ps + d_param.P_LO)*d_param.df; % shot noise variance
            detector_secondary_param.sigma_th = (4*Const.kB*d_param.T)/d_param.R_L*d_param.Fn*d_param.df; % thermal noise variance
            detector_secondary_param.R = diag(repmat(detector_secondary_param.sigma_th+detector_secondary_param.sigma_s,2,1)); %measurement noise covariance matrix
            detector_secondary_param.SNR = 10*log10(2*(d_param.eta_d*Const.q/(Const.h*l_param.nu))^2*(laser_secondary_param.Ps*d_param.P_LO)/(detector_secondary_param.sigma_th+detector_secondary_param.sigma_s));
            detector_secondary_param.omega = 2*pi*(l_param.nu - d_param.nu_LO);
            
            f1 = fieldnames(l_param);
            for i = 1:length(f1)
                theta.(f1{i}) = l_param.(f1{i});
                theta_primary.(f1{i}) = l_param.(f1{i});
            end
            f2 = fieldnames(d_param);
            for i = 1:length(f2)
                theta.(f2{i}) = d_param.(f2{i});
                theta_primary.(f2{i}) = d_param.(f2{i});
            end
            f3 = fieldnames(laser_secondary_param);
            for i = 1:length(f3)
                theta.(f3{i}) = laser_secondary_param.(f3{i});
                theta_secondary.(f3{i}) = laser_secondary_param.(f3{i});
            end
            f4 = fieldnames(detector_secondary_param);
            for i = 1:length(f4)
                theta.(f4{i}) = detector_secondary_param.(f4{i});
                theta_secondary.(f4{i}) = detector_secondary_param.(f4{i});
            end
            
        end
        function [PSD,f] = PSD(x,Ts) % power spectral density
            % PSD power spectral density of a signal
            % ----- INPUT ----.
            % x: input signal
            % Ts: sampling time
            % ----- OUTPUT -----
            % PSD: power spectral density
            % f: frequency
            if nargin <2
                Ts=1;
            end
            N = length(x);
            f = (0:N-1).'/N/Ts;
            PSD = abs(fft(x)*Ts).^2/N/Ts;
        end
        function [DSPSD,f] = DSPSD(x,Ts) %Double side power spectral density (rescaled)
            % DSPSD double side power spectral density of a signal
            % ----- INPUT ----.
            % x: input signal
            % Ts: sampling time
            % ----- OUTPUT -----
            % DSPSD: double side power spectral density
            % f: frequency
            if nargin <2
                Ts=1;
            end
            if (size(x,1)<size(x,2)) % if this is a row vector, turn into a column one
                x = x.';
            end
            [N,~] = size(x);
            PSD = abs(fft(x)*Ts).^2/N/Ts;
            DSPSD = fftshift(PSD);
            if mod(N,2)==0
                f=(-N/2:N/2-1).'./Ts/N; % N even
            else
                f=(-(N-1)/2:(N-1)/2).'./Ts/N; % N odd
            end
        end
        function [SSPSD,f_ss] = SSPSD(x,Ts,W,s) % Single side power spectral density
            if nargin <2
                Ts=1;
            end
            if (size(x,1)<size(x,2)) % if this is a row vector, turn into a column one
                x = x.';
            end
            [Nt,~] = size(x);
            
                       
            %             f_ss = (0:(1/Ts)/N:(1/Ts)/2).';
            
            if nargin < 3
                X = fft(x.*hann(length(x)));
                f = (0:Nt-1).'/Nt/Ts; %frequency vector
                SSPSD = abs(X(1:floor(Nt/2)+1,:)*Ts).^2/Nt/Ts;
            else
                if W > Nt
                    W = Nt;
                end
                K = floor((Nt-W)/s)+1;
                f = (0:W-1).'/W/Ts; %frequency vector
                SSPSD = zeros(floor(W/2)+1,1);
                for k = 1:K
                    x_k = x((k-1)*s + (1:W));
                    X_k = fft(x_k.*hann(length(x_k)));
                    SSPSD = (k-1)/k.*SSPSD + abs(X_k(1:floor(W/2)+1,:)*Ts).^2/W/Ts./k;
%                     disp(['Single Side PSD calculation: ',num2str(k/K*100),' % completed.'])
                end
            end
            f_ss = f(f <= 1/2/Ts);
            SSPSD(2:end-1,:) = 2*SSPSD(2:end-1,:);
        end
        % SSPSD single side power spectral density of a signal
        % ----- INPUT ----.
        % x: input signal
        % Ts: sampling time
        % ----- OUTPUT -----
        % SSPSD: single side power spectral density
        % f: frequency
        function [LW] = linewidth_from_FN(y,x,method)
            % LINEWIDTH_FROM_FN calculate the linewidth from the phase noise time
            %series, using the approximaiton formula from Di Domenico et al (2010).
            % ----- INPUT ----
            % y: signal (in time or frequency)
            % x: time or frequency
            % method:
            %   - 'time' y is the time trace of the phase, x is the time index of the samples (default)
            %   - 'frequency' y is the frequency noise spectra, x is the frequency index
            % ----- OUTPUT -----
            % LW: linewidth in hertz of FWHM
            if nargin < 2
                method = 'time';
            end
            [row,col] = size(y);
            if row < col
                y = y.';
            end
            if strcmp(method,'time') % operate in time domain - ideally, the sampling rate is constant
                y = detrend(y); % remove eventual constants - linear and costant trends that unaffect the linewidth
                [FN,f] = Laser.FN(y,x(2)); % get tFhe frequency noise
                FN = movmean(FN,length(FN)/30);% smooth out the FN
                beta_line = f.*(8*log(2)/pi^2);
                A = trapz(f(2:end),FN(2:end).*(FN(2:end) >= beta_line(2:end)));
                LW = sqrt(A.*8*log(2));
            elseif strcmp(method,'frequency') % operate in frequency domain, meaning that the input it is already the spectrum, which has been already averaged
                beta_line = x.*(8*log(2)/pi^2);
                A = trapz(x(1:end),y(1:end).*(y(1:end) >= beta_line(1:end)));
                LW = sqrt(A.*8*log(2));
            else
                error('Wrong mehtod. Choose if "time" or "frequency".')
            end
        end
        function [LW] = linewidth_from_PSD(y,res,method)
            % LINEWIDTH_FROM_PSD calculate the linewidth from the power spectral density of the signal.
            % it evaluates the peak, and then the width of the FWHM
            % ----- INPUT ----
            % y: signal (in time or frequency)
            % resolution: resolution in time or frequency
            % method:
            %   - 'time' y is the time trace of the signal, res is the sampling time
            %   - 'frequency' y is the PSD of the signal, res is the lowest representable frequency
            % ----- OUTPUT -----
            % LW: linewidth in hertz
            if strcmp(method,'time')
                N = length(y);
                f_res = 1/res/N;
                psd_y = Laser.SSPSD(y,res);
                maximum = max(psd_y);
                half_maximum = maximum/2;
                LW = sum(psd_y>=half_maximum)*(f_res);
            elseif strcmp(method,'frequency')
                maximum = max(y);
                half_maximum = maximum/2;
                LW = sum(y>=half_maximum)*(res);
            else
                error('Wrong mehtod. Choose if "time" or "frequency".')
            end
            
        end
        function [av] = allan_variance(x,Ts)
            % ALLAN_VARIANCE calcualte the allan variance of a signal
            % ----- INPUT ----
            % x: signal
            % Ts: sampling time
            % ----- OUTPUT -----
            % av: allan variance
            av = 0.5*mean(diff(x).^2);
        end
        function[f,RIN,FN] = analytic_spectrums_2(theta,f) % (from small signal analysis)
            % ANALYTIC_SPECTRUMS_2 return the analytic spectrums from small signal analysis of the rate equations. This theoretical spectrums are then an approximation of the spectra calculated on the generated internal states
            % this static method is faster than the other one
            % ----- INPUT ----.
            % theta: rate equation parameters
            % f: desired frequency vector
            % ----- OUTPUT -----
            % f: frequency vector
            % RIN: relative intensity noise
            % FN: frequency noise
            omega = 2*pi*f; % angular frequency
            DNN = 2*theta.N_ss*(1+theta.beta*theta.S_ss)/theta.tau_n;
            DSS = 2*theta.beta*theta.N_ss*theta.S_ss/theta.tau_n;
            DSN = -DSS;
            DPP = theta.beta*theta.N_ss/(2*theta.tau_n*theta.S_ss);
            Aa = 1i.*omega + 1/theta.tau_n + theta.g*theta.S_ss/(1+theta.eps*theta.S_ss);
            Ab = theta.g*(theta.N_ss - theta.N0)/((1+theta.eps*theta.S_ss)^2);
            Ac = -theta.g*theta.S_ss/(1+theta.eps*theta.S_ss)- theta.beta/theta.tau_n;
            Ad = 1i.*omega + 1/theta.tau_p - theta.g*(theta.N_ss - theta.N0)/((1 + theta.eps*theta.S_ss)^2);
            detA = Aa.*Ad - Ab.*Ac;
            Bw_11 = (1i.*omega + 1/theta.tau_p - theta.g*(theta.N_ss - theta.N0)/((1 + theta.eps*theta.S_ss)^2))./detA;
            Bw_21 = (theta.g*theta.S_ss/(1+theta.eps*theta.S_ss) + theta.beta/theta.tau_n)./detA;
            Bw_12 =  - theta.g*(theta.N_ss - theta.N0)/((1+theta.eps*theta.S_ss)^2)./detA;
            Bw_22 = (1i.*omega + 1/theta.tau_n + theta.g*theta.S_ss/(1+theta.eps*theta.S_ss))./detA;
            FN = 1/(4*pi^2)*((theta.alpha*theta.g/2)^2*(Bw_11.*conj(Bw_11)*DNN + Bw_12.*conj(Bw_12)*DSS + DSN*(Bw_11.*conj(Bw_12) + Bw_12.*conj(Bw_11))) + DPP);
            RIN = 10*log10((Bw_21.*conj(Bw_21)*DNN + Bw_22.*conj(Bw_22)*DSS + DSN*(Bw_21.*conj(Bw_22) + Bw_22.*conj(Bw_21)))/(theta.S_ss^2));
        end
        function[SS_1,SS_2] = get_steady_states(theta)
            % GET_STEADY_STATES return the steady states from the rate equation parameters
            % ----- INPUT ----.
            % theta: rate equation parameters
            % ----- OUTPUT -----
            % SS_1: positive steady states (feasible)
            % SS_2: negative steady states (unfeasible, only a mathematical solution)
            aa = Const.q*(theta.eps + theta.g*theta.tau_n);
            bb = Const.q + theta.g*theta.tau_p*((1-theta.beta)*theta.N0*Const.q-theta.tau_n*theta.I) - theta.beta*theta.eps*theta.tau_p*theta.I;
            cc = -theta.beta*theta.tau_p*theta.I;
            % positive (and valid) steady states
            SS_1(1)=(-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa); % ok
            SS_1(2)= theta.tau_n/Const.q*(theta.I*(1 + theta.eps*SS_1(1)) + theta.g*theta.N0*Const.q*SS_1(1))/(1 + theta.eps*SS_1(1) + theta.g*SS_1(1)*theta.tau_n); % ok
            SS_1(3) = 0;
            % negative steady states
            SS_2(1)=(-bb - sqrt(bb^2 - 4*aa*cc))/(2*aa); % ok
            SS_2(2)= theta.tau_n/Const.q*(theta.I*(1 + theta.eps*SS_2(1)) + theta.g*theta.N0*Const.q*SS_2(1))/(1 + theta.eps*SS_2(1) + theta.g*SS_2(1)*theta.tau_n); % ok
            SS_2(3) = 0;
            
        end
        function [tau] = get_time_constants(theta)
            % GET_TIME_CONSTANTS return the time constants from the rate equations. It maybe have some connection with numerical solution stability.
            % ----- INPUT -----
            % theta: rate equation parameters
            % ----- OUTPUT -----
            % tau: time constants
            tau(1,1) = (theta.I/Const.q/theta.N_ss)^-1;
            tau(1,2) = theta.tau_n;
            tau(1,3) = (theta.g*(theta.N_ss-theta.N0)*theta.S_ss/(1+theta.eps*theta.S_ss)/theta.N_ss)^-1;
            tau(1,4) = (sqrt(theta.N_ss*(1+theta.beta*theta.S_ss)/theta.tau_n)/theta.N_ss)^-1;
            tau(1,5) = (sqrt(2*theta.beta^2*theta.N_ss*theta.S_ss^2/(theta.tau_n*(1+theta.beta*theta.S_ss)))/theta.N_ss)^-1;
            
            tau(2,1) = (theta.g*(theta.N_ss-theta.N0)*theta.S_ss/(1+theta.eps*theta.S_ss)/theta.S_ss)^-1;
            tau(2,2) = theta.tau_p;
            tau(2,3) = (theta.N_ss*theta.beta/theta.tau_n/theta.S_ss)^-1;
            tau(2,4) = (sqrt(2*theta.beta^2*theta.N_ss*theta.S_ss^2/(theta.tau_n*(1+theta.beta*theta.S_ss)))/theta.S_ss)^-1;
            tau(2,5) = (sqrt(2*theta.beta*theta.N_ss*theta.S_ss/(theta.tau_n*(1+theta.beta*theta.S_ss)))/theta.S_ss)^-1;
            
            
        end
        function [C,Ck] = corr_time_series(x,W,s)
            % CORR_TIME_SERIES return the correlation matrix of a multivariate time series, for a given length of time samples
            % ----- INPUT -----
            % x: multidimensional time series
            % W: window in samples. The matrix is calculated on W samples in time
            % s: step. Advance the window over the time series, to smooth the results. move every s steps
            % ----- OUTPUT -----
            % C: correlation matrix
            [row,col] = size(x);
            if row < col
                x = x.';
                N = col;
                Nx = row;
            else
                N = row;
                Nx = col;
            end
            if nargin == 1
                W = N;
                s = 1;
            elseif nargin == 2
                s = 1;
            end
            s = max(1,floor(s));
            C = zeros(Nx,Nx);
            if W > N
                W = N;
            end
            K = floor((N-W)/s)+1;
            
                if nargout >1
            Ck = zeros(Nx,Nx,K);
                end
            for k = 1:K
                if nargout < 2
                C = (k-1)/k.*C + corr(x((k-1)*s + (1:W),:))./k;
                else
                Ck(:,:,k) = corr(x((k-1)*s + (1:W),:));
                C = (k-1)/k.*C + Ck(:,:,k)./k;
                end
            end
            
        end
        function [C,Ck] = cov_time_series(x,W,s)
            % COV_TIME_SERIES return the covariance matrix of a multivariate time series, for a given length of time samples
            % ----- INPUT -----
            % x: multidimensional time series
            % W: window in samples. The matrix is calculated on W samples in time
            % s: step. Advance the window over the time series, to smooth the results. move every s steps
            % ----- OUTPUT -----
            % C: covariance matrix
            [row,col] = size(x);
            if row < col
                x = x.';
                N = col;
                Nx = row;
            else
                N = row;
                Nx = col;
            end
            if nargin == 1
                W = N;
                s = 1;
            elseif nargin == 2
                s = 1;
            end
            s = max(1,floor(s));
            C = zeros(Nx,Nx);
            if W > N
                W = N;
            end
            K = floor((N-W)/s)+1;
            if nargout >1
            Ck = zeros(Nx,Nx,K);
            end
            for k = 1:K
                if nargout < 2
                C = (k-1)/k.*C + cov(x((k-1)*s + (1:W),:))./k;
                else
                Ck(:,:,k) = cov(x((k-1)*s + (1:W),:));
                C = (k-1)/k.*C + Ck(:,:,k)./k;
                end
            end
            
        end
        function [xc,xck] = xcorr_time_series(x,y,W,s)
            % XCORR_TIME_SERIES return the cross correlation estimation of two 1-dimensional time series, for a given length of time samples
            % ----- INPUT -----
            % x,y : time series
            % W: window in samples. The matrix is calculated on W samples in time. if missing, the whole sequences are used, i.e. W = min(lenght(x),lenght(y))
            % s: step. Advance the window over the time series, to smooth the results. move every s steps. If missing
            % ----- OUTPUT -----
            % xc: cross correlation between x and y of length W samples
            if ~iscolumn(x)
                x = x.';
            end
            if ~iscolumn(y)
                y = y.';
            end
            N = max(length(x),length(y));
             if nargin == 2
                W = N;
                s = W;
            elseif nargin == 3
                s = W;
             end
            if W > N
                W = N;
            end
            s = max(1,floor(s));
            K = floor((N-W)/s)+1;
            xc = 0;
            if nargout >1
                xck = zeros(1,K);
            end
            for k = 1:K
                if nargout < 2
                xc = (k-1)/k.*xc + corr(x((k-1)*s + (1:W),:),y((k-1)*s + (1:W),:))./k;
                else
                xck(k) = corr(x((k-1)*s + (1:W),:),y((k-1)*s + (1:W),:));
                xc = (k-1)/k.*xc + xck(k)./k;
                end
            end
        end
        function [xC,xCk] = xcov_time_series(x,y,W,s)
            % XCOV_TIME_SERIES return the cross covariance estimation of two 1-dimensional time series, for a given length of time samples
            % ----- INPUT -----
            % x,y : time series
            % W: window in samples. The matrix is calculated on W samples in time. if missing, the whole sequences are used, i.e. W = min(lenght(x),lenght(y))
            % s: step. Advance the window over the time series, to smooth the results. move every s steps. If missing
            % ----- OUTPUT -----
            % xC: cross covariance between x and y of length W samples
            if ~iscolumn(x)
                x = x.';
            end
            if ~iscolumn(y)
                y = y.';
            end
            N = max(length(x),length(y));
             if nargin == 2
                W = N;
                s = W;
            elseif nargin == 3
                s = W;
             end
            if W > N
                W = N;
            end
            s = max(1,floor(s));
            K = floor((N-W)/s)+1;
            xC = 0;
            if nargout >1
                xCk = zeros(1,K);
            end
            for k = 1:K
                cov_temp = cov(x((k-1)*s + (1:W),:),y((k-1)*s + (1:W),:));
                if nargout < 2
                   xC = (k-1)/k.*xC + cov_temp(1,2)./k; 
                else
                   xCk(k) = cov_temp(1,2);
                   xC = (k-1)/k.*xC + cov_temp(1,2)./k; 
                end
                
            end
        end
        
        function [y,Ts] = generate_process(PSD,f_PSD,varargin)
            % GENERATE_PROCESS % generate a random process (in time) from a given PSD. The generated process will have the desired PSD
            % the desired psd can have few points, the rest are going to be interpolated using a spline fit
            % ----- INPUT -----
            % PSD: desired psd
            % f_PSD: frequency corresponding to PSD
            % varargin: variable number of input arguments. a bad combination can affect the quality of the output
            % 'Samples' : desired number of samples
            % 'SamplingTime' : desired sampling time
            % setting only one of this will automatically adjust the other one
            % ----- OUTPUT -----
            % y: generated process
            % Ts: effective sampling time
            
            % check for vertical/horizontal and return the same shape
            if sum(diff(f_PSD) < 0) > 0
                error('Frequency values must be strictly non descreasing')
            end
            if (size(PSD,1)<size(PSD,2)) % if this is a row vector, turn into a column one
                PSD = PSD.';
            end
            if (size(f_PSD,1)<size(f_PSD,2)) % if this is a row vector, turn into a column one
                f_PSD = f_PSD.';
            end
            f_nyquist_ideal = f_PSD(end); % Since we have to represent this frequency
            f_min_ideal = f_PSD(1); % lowest representable frequency
            Ts_ideal = 1/2/f_nyquist_ideal;
            N_ideal = 1/Ts_ideal/f_min_ideal;
            
            
            
            if isempty(varargin)% N and Ts are determined by the frequency vector
                f_nyquist = f_nyquist_ideal; % Since we have to represent this frequency
                f_min = f_min_ideal; % lowest representable frequency
                Ts =Ts_ideal;
                N = N_ideal;
            elseif length(varargin) == 2
                if strcmp(varargin{1},'Samples') % Ts is determined by the frequency vector, but N is choosen
                    N = varargin{2}; % sacrificing the lowest frequency part
                    Ts =Ts_ideal;
                    f_nyquist = f_nyquist_ideal;
                    f_min = 1/Ts/N;
                    if N < N_ideal
                        warning(['Raising the lowest frequency resolution. Lowest representable frequency : ', num2str(f_min), ' Hz';...
                            ' It is advised to use ',num2str(N_ideal), ' samples to represent ',num2str(f_min_ideal), ' Hz'])
                    end
                elseif strcmp(varargin{1},'SamplingTime')% N is determined by the frequency vector, but Ts is choosen
                    Ts = varargin{2};
                    N = N_ideal;
                    f_nyquist = 1/Ts/N;
                    f_min = f_min_ideal;
                    if Ts > Ts_ideal
                        warning(['Lowering the highest frequency resolution. Highest representable frequency : ', num2str(f_nyquist), ' Hz';...
                            ' It is advised to use ',num2str(Ts_ideal), ' s to represent ',num2str(f_nyquist_ideal), ' Hz'])
                    end
                else
                    error('Invalid parameter.')
                end
            elseif length(varargin) == 4 % custom N and Ts
                if strcmp(varargin{1},'Samples') && strcmp(varargin{3},'SamplingTime')
                    Ts = varargin{4};
                    N = varargin{2};
                elseif strcmp(varargin{3},'Samples') && strcmp(varargin{1},'SamplingTime')
                    Ts = varargin{2};
                    N = varargin{4};
                    
                else
                    error('Invalid parameters.')
                end
                if N < N_ideal
                    warning(['Raising the lowest frequency resolution. Lowest representable frequency : ', num2str(f_min), ' Hz';...
                        ' It is advised to use ',num2str(N_ideal), ' samples to represent ',num2str(f_min_ideal), ' Hz'])
                end
                if Ts > Ts_ideal
                    warning(['Lowering the highest frequency resolution. Highest representable frequency : ', num2str(f_nyquist), ' Hz';...
                        ' It is advised to use ',num2str(Ts_ideal), ' s to represent ',num2str(f_nyquist_ideal), ' Hz'])
                end
                f_nyquist = 1/2/Ts;
                f_min = 1/Ts/N;
            else
                error('Invalid parameters.')
            end
            
            
            
            
            f_half =  (0:f_min:f_nyquist).';
            
            
            % fit a spiline to get the whole PSD
            SplineFit = fit([0;f_PSD;f_nyquist], log([PSD(1); PSD; PSD(end)]), 'pchip');
            %               SplineFit = fit(f_PSD, (PSD), 'pchip');
            % evaluate such curve on frequency query points
            log_PSD_y_SS = feval(SplineFit,f_half);
            PSD_y_SS = exp(log_PSD_y_SS);
            %               PSD_y_SS = feval(SplineFit,f_half);
            
            
            % suppose that PSD_y is single side PSD
            PSD_y_DS = [0.5*PSD_y_SS;flipud(PSD_y_SS(2:end))*0.5];
            Nf = length(PSD_y_DS);
            %               sigma_x = sqrt(mean(PSD_y_DS));
            sigma_x = sqrt(1/Ts);
            if Nf > N
                x = sigma_x*randn(Nf,1);
                X = fft(x);
                y = ifft(X.*sqrt(PSD_y_DS));
                y = y(1:N);
            else
                x = sigma_x*randn(N,1);
                y = zeros(size(x));
                s = floor(Nf/2); % improve computational speed. But It can be changed if we want
                % we apply window filtering
                K = floor((N-Nf)/s)+1;
                for k = 1:K
                    y((k-1)*s + (1:Nf)) = ifft(fft(x((k-1)*s + (1:Nf))).*sqrt(PSD_y_DS));
                    disp([num2str(k/K*100),' % completed.'])
                end
            end
            
        end
        function R = get_correlation(SIGMA)
            D = (diag(SIGMA));
            invD = diag(D.^-0.5);
            R = invD*SIGMA*invD;
        end
    end
end

