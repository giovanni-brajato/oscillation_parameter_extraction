function [y,Ts] = color_noise_generator(PSD,f_PSD,varargin)
% COLOR_NOISE_GENERATOR % generate a random process (in time) from a given PSD. The generated process will have the desired PSD. It works by "coloring" white noise accordingly to the given power spectral density
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

