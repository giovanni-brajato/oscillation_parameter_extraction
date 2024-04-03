classdef correlation_analysis
    %CORRELATION_ANALYSIS This class contains methods for calculate the
    %correlation among time series traces and perform the multidimensional PCA
    %for different time durations
    
    
    properties
        Fs; % sampling frequency of the time seriers
        win_step;% window step advancement in data samples
        win_length;% window lengths in data samples (can be more than one)
        win_type;% windowing function for the data (optional)
        mode;% choose what operation to perform, 'corr' as correlation or 'cov' as covariance
    end
    
    methods
        function obj = correlation_analysis(param)
            %CORRELATION_ANALYSIS Construct an instance of this class
            obj.Fs = param.Fs; % 
            if length(param.win_step) < length(param.win_length)
                if length(param.win_step) == 1
                    obj.win_step = param.win_step.*ones(length(param.win_length),1); 
                    obj.win_length = param.win_length;
                else
                    error('Window lenght and step are not conistent')
                end
            elseif length(param.win_step) == length(param.win_length)
                obj.win_step = param.win_step; 
                obj.win_length = param.win_length;
            else
                error('Window lenght and step are not conistent')
            end
            obj.win_type = param.win_type; 
            obj.mode = param.mode; 
        end
        function [C,tau,n_avg] = covariance_time_series(obj,x)
            % CORR_TIME_SERIES return the correlation (or covariance) matrix of a multivariate time series, for a given length of time samples
            % ----- INPUT -----
            % x: multidimensional time series
            % ----- OUTPUT -----
            % C: correlation (or covariance) matrix for every defined
            % observation time
            % tau: observation times
            % n_avg: number of time each C is estimated for a given window
            % length
            [row,col] = size(x);
            if row < col
                x = x.';
                Nt = col; % time index
                Nx = row; % channel index
            else
                Nt = row;
                Nx = col;
            end
            
            Nw = length(obj.win_length);
            C = zeros(Nx,Nx,Nw);
            n_avg = zeros(Nw,1);
            tau = zeros(Nw,1);
            for w = 1:Nw
                s = obj.win_step(w);
                wl = obj.win_length(w);
                tau(w) = wl/obj.Fs;
                K = floor((Nt-wl)/s)+1;
                n_avg(w) = K;
                for k = 1:K
                        switch obj.mode
                            case 'cov'
                                C(:,:,w) = (k-1)/k.*C(:,:,w) + cov(obj.win_type(size(x((k-1)*s + (1:wl),:))).*x((k-1)*s + (1:wl),:))./k;
                            case 'corr'
                                C(:,:,w) = (k-1)/k.*C(:,:,w) + corr(obj.win_type(size(x((k-1)*s + (1:wl),:))).*x((k-1)*s + (1:wl),:))./k;
                        end
                        
                end
            end
            
            
            
        end
        function [eigVa,eigVe,C,tau] = eig_analysis(obj,x)
            % EIG_ANALYSIS return the evolution of the eigenvalues and
            % eigenvectors over time
            % ----- INPUT -----
            % x: multidimensional time series
            % ----- OUTPUT -----
            % eigVa: eigenvalues over observation times. dimension 1 is the
            % space dimension,dimension 2 is the time dimension
            % eigVe: eigenvectors over observation times. dimension 1 is the
            % space dimension,dimension 2 is the time dimension, dimension
            % 3 is the eigenvector number
            % C: the covariance matrix,3 is the observation time dimension
            % tau: observation times
            [row,col] = size(x);
            if row < col
                x = x.';
                Nt = col; % time index
                Nx = row; % channel index
            else
                Nt = row;
                Nx = col;
            end
            [C,tau] = obj.covariance_time_series(x);
            Nw = size(C,3);
            eigVa = zeros(Nx,Nw);
            eigVe = zeros(Nx,Nw,Nx);
            for w = 1:Nw
                [V,D] = eig(C(:,:,w));
                [eigVa(:,w),indsort] = sort(diag(D));
                for d = 1:Nx
                eigVe(:,w,d) = V(:,indsort(d));
                end
            end
        end
    end  
    methods(Static)
        
    end
end

