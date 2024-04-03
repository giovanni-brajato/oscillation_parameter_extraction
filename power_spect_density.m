classdef power_spect_density
    %POWER_SPECT_DENSITY an object of this class has methods for different
    %power spectral density conputations
    
    properties
        Fs; % sampling frequency
        QL_LW; % quantum limited linewidth
        Noise_floor; % detector noise floor
    end
    
    methods
        function obj = power_spect_density(param) % define the object by assigning values to its properties
            obj.Fs = param.Fs;
            obj.QL_LW = param.QL_LW;
            obj.Noise_floor = param.Noise_floor;
        end
        
        function [freq,PSD_phase_noise] = phase_noise(obj,X) 
            % calculate the power spectral density of a signal X (in this
            % case assumed to be the phase noise). The frequency vector is
            % only positive, and with the same length as the input signal X
            if ~iscolumn(X)
                X = X';
            end
            [L,freq] = freq_axis(obj,X);
            T_end = (1/obj.Fs)*L;
            Ts = 1/obj.Fs;
            L = length(X);
            W = hann(L); % Hanning window
            if  ~iscolumn(W)
                W = W';
            end
            PSD_pn = abs(fft(X.*W)*Ts).^2;
            PSD_phase_noise = 1/T_end*PSD_pn;
        end
        
        function [freq, CPSD] = CPSD(obj,X,Y)
        % calculate the cross power spectral density between the signals X and Y
        if numel(X) ~= numel(Y)
           error('X and Y must have the same size.') 
        end
         if ~iscolumn(X)
                X = X';
         end
         if ~iscolumn(Y)
                Y = Y';
         end
           [L,freq] = freq_axis(obj,X);
            T_end = (1/obj.Fs)*L;
            Ts = 1/obj.Fs;
            L = length(X);
            W = hann(L); % Hanning window
            CPSD = 1/T_end*(fft(X.*W)*Ts).*conj(fft(Y.*W)*Ts);
        end
        
        
        function [freq, RIN] = rel_intensity_noise(obj,X)
            % Calculate the relative intensity noise, starting from the
            % signal power X. X has to be the power of the signal over
            % time, not the amplitude
            [L,freq] = freq_axis(obj,X);
            T_end = (1/obj.Fs)*L;
            Ts = 1/obj.Fs;
            RIN = (1/(mean(X).^2)).*(1/(T_end).*abs(fft( (X - mean(X)).*Ts )).^2);
        end
        
        function [L,freq] = freq_axis(obj,X)
            % return the fft frequency representation from the signal
            % length L and the sampling frequency obj.Fs
            L = length(X);
            freq = (0 : L-1)*obj.Fs/L;
        end
        
        function [PSD_Lorentz] = Lorentzian(obj,freq)
            % return the Lorentzian phase noise profile given the frequency
            % vector freq and the quantum limited linewidth obj.QL_LW
            PSD_Lorentz = (obj.QL_LW./(2*pi*(freq.^2 + (obj.QL_LW/2)^2)))';
        end
        
        function PSD_noise_floor = Noise(obj,freq)
            % return the photodetector noise floor
            PSD_noise_floor = (obj.Noise_floor*ones(1,length(freq)))';
        end
    end
end

