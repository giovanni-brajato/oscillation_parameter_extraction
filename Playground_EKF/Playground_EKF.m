clear all
close all

format short
% addpath /zhome/72/8/113687/Desktop/Common_files
%% initialization
Fs = 50e9; % sampling frequency
f_nyquist = Fs/2;
Ts = 1/Fs; % sampling time
N =  1e5; % number of data points (maximum is 62.5e6)
t = (0:N-1).'.*Ts; % time vector
f = (0:N-1).'./Ts/N; % frequency vector



% for spectra visualizations
movmean_factor = 1; % for smoothing out spectra


% EM parameter initialization
R0 = 2.2e-6;
Q0p = 1.1e-8;
Q0a = 4.2e-13;


% noise correlation parameters
N_max = N;
N_min = 1e2;

T_slow = N_max*Ts;
T_fast = N_min*Ts;

f_high = 1/T_fast;
f_low = 1/T_slow;

N_points = 100;
spanlengths = round(logspace(log10(N_min),log10(N_max),N_points));
observation_times = spanlengths.*Ts;
observation_frequencies = 1./observation_times;
step_factor = N; % window length /  step advancement

saving_data = 0;
saving_figures = 1;
calculate_correlations = 1;

trace =3; % from 3 to 80 records
trace_number = (num2str(trace.','%03u'));
% addpath '/work3/gibra/laser_measurements/'
addpath 'O:\Ftnk\Communication-MLPS\UCSB MEASUREMENTS'

RIN_correction_factor_dB = 0;
PN_correction_factor_dB = 0;

%% data import

% channel 1
filename = ['Morton_K_',trace_number,'_Ch1.h5'];
data_info = h5info(filename);
S_i = h5read(filename,'/Data');
% channel 2
%     filename = ['Morton_K_',current_trace_number,'_Ch2.h5'];
%     data_info = h5info(filename);
%     S_q = h5read(filename,'/Data');
%     Ts = h5read(filename,'/Spacing');

%     S = S_i+1i*S_q;
S = S_i;
clear S_i

%% import laser data for comparison
pn_traces = importdata('MORTON_391MA_PHASE');
freq_pn_ref = pn_traces.data(:,1);
PSD_pn_ref = pn_traces.data(:,2)+3;

rin_traces = importdata('MORTON_391MA_RIN');
freq_rin_ref = rin_traces.data(:,1);
PSD_rin_ref = rin_traces.data(:,2);

%% data preparation
% locate the beat frequency
pwr = abs(fft(S));
f_ext = (0:length(S)-1).'/length(S)/Ts;
[~, index] = max(pwr(2:end));
beat_f = f_ext(index);
max_p = pwr(index);
clear pwr f_ext




%% EM init
theta_EM.R = R0;
theta_EM.Q = [Q0a,0;0,Q0p];
theta_EM.dt = Ts;
theta_EM.A = 1;
theta_EM.H = 1;
theta_EM.omega = 2*pi*beat_f;
theta_EM.a = double(mean(abs(hilbert(S))));
EM_param_track = [];
EM_param_track.monitor = 0;
m0p = 0;
m0a = 0;
m0 = [m0a,m0p];
P0 = theta_EM.Q ;
prev_E = + inf;
DE = 0;
saved_param = theta_EM;

%% set up figure
% Create UIFigure
UIFigure = uifigure('Name','EKF parameter analysis v0.1','NumberTitle','off','Position',[100 100 866 698],'Visible','off');
% Create GridLayout
GridLayout = uigridlayout(UIFigure);
GridLayout.ColumnWidth = {'1x', '1x', '1x'};
GridLayout.RowHeight = {'1x', '1x', '1x', '1x'};

% Create Panel
Panel = uipanel(GridLayout);
Panel.Layout.Row = [1 2];
Panel.Layout.Column = 3;

% Create R_min
R_min = uieditfield(Panel, 'numeric');
R_min.Position = [8 360 41 22];

% Create R_max
R_max = uieditfield(Panel, 'numeric');
R_max.Position = [167 362 41 22];

% Create Min_Label
Min_Label = uilabel(Panel);
Min_Label.Position = [8 381 44 22];
Min_Label.Text = 'Min.';

% Create Max_Label
Max_Label = uilabel(Panel);
Max_Label.Position = [167 384 51 22];
Max_Label.Text = 'Max.';

% Create Curr_Label
Curr_Label = uilabel(Panel);
Curr_Label.Position = [238 383 46 22];
Curr_Label.Text = 'Current';

% Create Qp_min
Qp_min = uieditfield(Panel, 'numeric');
Qp_min.Position = [8 334 41 22];

% Create Qp_max
Qp_max = uieditfield(Panel, 'numeric');
Qp_max.Position = [167 335 41 22];

% Create Qa_min
Qa_min = uieditfield(Panel, 'numeric');
Qa_min.Position = [8 308 41 22];

% Create Qa_max
Qa_max = uieditfield(Panel, 'numeric');
Qa_max.Position = [167 307 41 22];

% Create Qap_max
Qap_max = uieditfield(Panel, 'numeric');
Qap_max.Position = [167 279 41 22];

% Create Qap_min
Qap_min = uieditfield(Panel, 'numeric');
Qap_min.Position = [8 280 41 22];

% Create omega_max
omega_max = uieditfield(Panel, 'numeric');
omega_max.Position = [167 251 41 22];

% Create omega_min
omega_min = uieditfield(Panel, 'numeric');
omega_min.Position = [8 252 41 22];

% Create A_min
A_min = uieditfield(Panel, 'numeric');
A_min.Position = [8 224 41 22];

% Create A_max
A_max = uieditfield(Panel, 'numeric');
A_max.Position = [167 223 41 22];

% Create dt_min
dt_min = uieditfield(Panel, 'numeric');
dt_min.Position = [8 196 41 22];

% Create dt_max
dt_max = uieditfield(Panel, 'numeric');
dt_max.Position = [167 193 41 22];

% Create RSliderLabel
RSliderLabel = uilabel(Panel);
RSliderLabel.HorizontalAlignment = 'right';
RSliderLabel.Position = [143 364 13 22];
RSliderLabel.Text = 'R';

% Create R_Slider
R_Slider = uislider(Panel);
R_Slider.MajorTicks = [];
R_Slider.MinorTicks = [];
R_Slider.Position = [60 373 75 3];
R_Slider.Value = theta_EM.R;
R_Slider.Limits = [0 1]; 

% Create QpSliderLabel
QpSliderLabel = uilabel(Panel);
QpSliderLabel.HorizontalAlignment = 'right';
QpSliderLabel.Position = [135 335 25 22];
QpSliderLabel.Text = 'Qp';

% Create Qp_Slider
Qp_Slider = uislider(Panel);
Qp_Slider.MajorTicks = [];
Qp_Slider.MinorTicks = [];
Qp_Slider.Position = [61 344 75 3];
Qp_Slider.Value = theta_EM.Q(2,2);
Qp_Slider.Limits = [0 1]; 

% Create QaSliderLabel
QaSliderLabel = uilabel(Panel);
QaSliderLabel.HorizontalAlignment = 'right';
QaSliderLabel.Position = [135 309 25 22];
QaSliderLabel.Text = 'Qa';

% Create Qa_Slider
Qa_Slider = uislider(Panel);
Qa_Slider.MajorTicks = [];
Qa_Slider.MinorTicks = [];
Qa_Slider.Position = [61 318 75 3];
Qa_Slider.Value = theta_EM.Q(1,1);
Qa_Slider.Limits = [0 1]; 

% Create QapSliderLabel
QapSliderLabel = uilabel(Panel);
QapSliderLabel.HorizontalAlignment = 'right';
QapSliderLabel.Position = [138 280 28 22];
QapSliderLabel.Text = 'Qap';

% Create Qap_Slider
Qap_Slider = uislider(Panel);
Qap_Slider.MajorTicks = [];
Qap_Slider.MinorTicks = [];
Qap_Slider.Position = [60 289 75 3];
Qap_Slider.Value = theta_EM.Q(1,2);
Qap_Slider.Limits = [0 1]; 

% Create omegaSliderLabel
omegaSliderLabel = uilabel(Panel);
omegaSliderLabel.HorizontalAlignment = 'right';
omegaSliderLabel.Position = [124 252 42 22];
omegaSliderLabel.Text = 'omega';

% Create omega_Slider
omega_Slider = uislider(Panel);
omega_Slider.MajorTicks = [];
omega_Slider.MinorTicks = [];
omega_Slider.Position = [60 261 65 3];
omega_Slider.Limits = [0 1e20]; 
omega_Slider.Value = theta_EM.omega;

% Create ASliderLabel
ASliderLabel = uilabel(Panel);
ASliderLabel.HorizontalAlignment = 'right';
ASliderLabel.Position = [132 224 25 22];
ASliderLabel.Text = 'A';

% Create A_Slider
A_Slider = uislider(Panel);
A_Slider.MajorTicks = [];
A_Slider.MinorTicks = [];
A_Slider.Position = [61 233 75 3];

A_Slider.Value = theta_EM.a;

% Create dtSliderLabel
dtSliderLabel = uilabel(Panel);
dtSliderLabel.HorizontalAlignment = 'right';
dtSliderLabel.Position = [133 197 25 22];
dtSliderLabel.Text = 'dt';

% Create dt_Slider
dt_Slider = uislider(Panel);
dt_Slider.MajorTicks = [];
dt_Slider.MinorTicks = [];
dt_Slider.Position = [62 206 75 3];
dt_Slider.Value = theta_EM.dt;

% Create R_curr
R_curr = uispinner(Panel);
R_curr.Step = 0.001;
R_curr.Position = [217 362 100 22];
R_curr.Value = theta_EM.R;

% Create Qp_curr
Qp_curr = uispinner(Panel);
Qp_curr.Step = 0.001;
Qp_curr.Position = [217 336 100 22];
Qp_curr.Value = theta_EM.Q(2,2);

% Create Qa_curr
Qa_curr = uispinner(Panel);
Qa_curr.Step = 0.001;
Qa_curr.Position = [217 308 100 22];
Qa_curr.Value = theta_EM.Q(1,1);

% Create Qap_curr
Qap_curr = uispinner(Panel);
Qap_curr.Step = 0.001;
Qap_curr.Position = [217 280 100 22];
Qap_curr.Value = theta_EM.Q(1,2);

% Create omega_curr
omega_curr = uispinner(Panel);
omega_curr.Step = 0.001;
omega_curr.Position = [217 251 100 22];
omega_curr.Value = theta_EM.omega;

% Create A_curr
A_curr = uispinner(Panel);
A_curr.Step = 0.001;
A_curr.Position = [217 223 100 22];
A_curr.Value = theta_EM.a;

% Create dt_curr
dt_curr = uispinner(Panel);
dt_curr.Step = 0.001;
dt_curr.Position = [217 196 100 22];
dt_curr.Value = theta_EM.dt;


% Create StartStop_Button
Exit_Button = uibutton(Panel, 'state');
Exit_Button.Text = 'Exit';
Exit_Button.Position = [253 162 54 22];


% Create FilteringLabel
FilteringLabel = uilabel(Panel);
FilteringLabel.HorizontalAlignment = 'center';
FilteringLabel.FontColor = [0 0 1];
FilteringLabel.Position = [166.5 88 48 22];
FilteringLabel.Text = 'Filtering';

% Create Track_Switch
Track_Switch = uiswitch(Panel, 'slider');
Track_Switch.FontColor = [0 0 1];
Track_Switch.Position = [90 89 45 20];
Track_Switch.Enable = 0;

% Create Reset_Button
Reset_Button = uibutton(Panel, 'state');
Reset_Button.Text = 'Reset';
Reset_Button.Position = [183 162 55 22];


% Create SmoothingLabel
SmoothingLabel = uilabel(Panel);
SmoothingLabel.HorizontalAlignment = 'center';
SmoothingLabel.FontColor = [1 0 0];
SmoothingLabel.Position = [162 53 63 22];
SmoothingLabel.Text = 'Smoothing';

% Create Track_Switch_2
Smooth_Switch = uiswitch(Panel, 'slider');
Smooth_Switch.FontColor = [1 0 0];
Smooth_Switch.Position = [87 54 45 20];
Smooth_Switch.Enable = 0;


% Create SimulateSwitchLabel
SimulateSwitchLabel = uilabel(Panel);
SimulateSwitchLabel.HorizontalAlignment = 'center';
SimulateSwitchLabel.FontColor = [1 0 1];
SimulateSwitchLabel.Position = [164 14 52 22];
SimulateSwitchLabel.Text = 'Simulate';

% Create Simulate_Switch
Simulate_Switch = uiswitch(Panel, 'slider');
Simulate_Switch.FontColor = [1 0 1];
Simulate_Switch.Position = [89 15 45 20];

% Create Usesliders_CheckBox
Usesliders_CheckBox = uiswitch(Panel, 'toggle');
Usesliders_CheckBox.Items = {'Spinners', 'Sliders'};
Usesliders_CheckBox.Orientation = 'horizontal';
Usesliders_CheckBox.Position = [68 163 45 20];
Usesliders_CheckBox.Value = 'Spinners';

% Create VVtrackingLabel
VVtrackingLabel = uilabel(Panel);
VVtrackingLabel.HorizontalAlignment = 'center';
VVtrackingLabel.Position = [161 125 67 22];
VVtrackingLabel.Text = 'VV tracking';

% Create Conv_Switch
Conv_Switch = uiswitch(Panel, 'slider');
Conv_Switch.Position = [90 125 45 20];

% Create Tracking_Samples_EditField_2
N_Samples_EditField = uieditfield(Panel, 'numeric');
N_Samples_EditField.Position = [237 124 79 22];
N_Samples_EditField.Value = N;

% Create ResultsPanel
ResultsPanel = uipanel(GridLayout);
ResultsPanel.Title = 'Results';
ResultsPanel.Layout.Row = 3;
ResultsPanel.Layout.Column = 3;

% Create E_Label
E_Label = uilabel(ResultsPanel);
E_Label.Position = [31 98 25 22];
E_Label.Text = '-LL';

% Create DE_Label
DE_Label = uilabel(ResultsPanel);
DE_Label.Position = [31 67 25 22];
DE_Label.Text = 'dLL';

% Create E
E = uilabel(ResultsPanel);
E.Position = [81 98 100 22];
E.Text = num2str(nan);
% Create DE
DE = uilabel(ResultsPanel);
DE.Position = [81 67 100 22];
DE.Text = num2str(nan);
% Create LW_Label
LW_Label = uilabel(ResultsPanel);
LW_Label.Position = [31 38 25 22];
LW_Label.Text = 'LW';

% Create LW
LW = uilabel(ResultsPanel);
LW.Position = [81 38 100 22];
LW.Text = num2str(nan);
% Create Theoreticalspectra_CheckBox
Theoreticalspectra_CheckBox = uicheckbox(ResultsPanel);
Theoreticalspectra_CheckBox.Text = {'Theoretical spectra +',' measurement noise level'};
Theoreticalspectra_CheckBox.FontColor = [0 1 1];
Theoreticalspectra_CheckBox.Position = [160 38 124 22];

% Create Hz_Label
Hz_Label = uilabel(ResultsPanel);
Hz_Label.Position = [140 38 20 22];
Hz_Label.Text = 'Hz';

% Create ParameterlearningPanel
ParameterlearningPanel = uipanel(GridLayout);
ParameterlearningPanel.Title = 'Parameter learning ';
ParameterlearningPanel.Layout.Row = 4;
ParameterlearningPanel.Layout.Column = 3;

% Create Gradientdescent_CheckBox
Gradientdescent_CheckBox = uicheckbox(ParameterlearningPanel);
Gradientdescent_CheckBox.Text = 'Gradient descent';
Gradientdescent_CheckBox.Position = [30 139 114 22];
Gradientdescent_CheckBox.Enable = 0;
% Create Expectationmaximization_CheckBox
Expectationmaximization_CheckBox = uicheckbox(ParameterlearningPanel);
Expectationmaximization_CheckBox.Text = 'Expectation maximization';
Expectationmaximization_CheckBox.Position = [30 108 158 22];
Expectationmaximization_CheckBox.Enable = 0;

% Create Refineomega_CheckBox
Refineomega_CheckBox = uicheckbox(ParameterlearningPanel);
Refineomega_CheckBox.Text = 'Refine omega';
Refineomega_CheckBox.Position = [30 66 97 22];
Refineomega_CheckBox.Enable = 0;
% Create RefineA_CheckBox
RefineA_CheckBox = uicheckbox(ParameterlearningPanel);
RefineA_CheckBox.Text = 'Refine A';
RefineA_CheckBox.Position = [127 66 68 22];
RefineA_CheckBox.Enable = 0;

% Create Signal_time
Signal_time = uiaxes(GridLayout);
title(Signal_time, 'Signal')
xlabel(Signal_time, 'time [s]')
ylabel(Signal_time, '[a.u.]')
Signal_time.Layout.Row = 1;
Signal_time.Layout.Column = 1;

% Create Signal_frequency
Signal_frequency = uiaxes(GridLayout);
title(Signal_frequency, 'Signal')
xlabel(Signal_frequency, 'frequency [Hz]')
ylabel(Signal_frequency, '[dB/Hz]')
Signal_frequency.Layout.Row = 1;
Signal_frequency.Layout.Column = 2;

% Create Amplitude_time
Amplitude_time = uiaxes(GridLayout);
title(Amplitude_time, 'Amplitude noise')
xlabel(Amplitude_time, 'time [s]')
ylabel(Amplitude_time, '[a.u.]')
Amplitude_time.Layout.Row = 2;
Amplitude_time.Layout.Column = 1;

% Create RIN
RIN = uiaxes(GridLayout);
title(RIN, 'RIN')
xlabel(RIN, 'frequency [Hz]')
ylabel(RIN, '[dB/Hz]')
RIN.XMinorTick = 'on';
RIN.XGrid = 'on';
RIN.XMinorGrid = 'on';
RIN.YGrid = 'on';
RIN.YMinorGrid = 'on';
RIN.XScale = 'log';
RIN.Layout.Row = 3;
RIN.Layout.Column = 1;

% Create PN
PN = uiaxes(GridLayout);
title(PN, 'PN')
xlabel(PN, 'frequency [Hz]')
ylabel(PN, '[dB rad^2/Hz]')
PN.XMinorTick = 'on';
PN.XGrid = 'on';
PN.XMinorGrid = 'on';
PN.YGrid = 'on';
PN.YMinorGrid = 'on';
PN.XScale = 'log';
PN.Layout.Row = 3;
PN.Layout.Column = 2;

% Create amp_autocorr
amp_autocorr = uiaxes(GridLayout);
title(amp_autocorr, 'Amp. noise autocorr')
xlabel(amp_autocorr, 'lag [s]')
ylabel(amp_autocorr, '')
amp_autocorr.Layout.Row = 4;
amp_autocorr.Layout.Column = 1;

% Create phase_autocorr
phase_autocorr = uiaxes(GridLayout);
title(phase_autocorr, 'Phase noise autocorr.')
xlabel(phase_autocorr, 'lag [s]')
ylabel(phase_autocorr, '')
phase_autocorr.Layout.Row = 4;
phase_autocorr.Layout.Column = 2;

% Create Phase_time
Phase_time = uiaxes(GridLayout);
title(Phase_time, 'Phase noise')
xlabel(Phase_time, 'time [s]')
ylabel(Phase_time, '[rad]')
Phase_time.Layout.Row = 2;
Phase_time.Layout.Column = 2;


%%
UIFigure.Visible = 'on';
warning off
exit = 0;
    y = double(S(1:N));
    y = y - mean(y);
    [sspsd_y,fss] = Laser.SSPSD(y,Ts);
in_lag = N/10;
        out_lag = N;
        window_size = N;
        step_length = N;
while ~exit
    % assigning some values
    if N_Samples_EditField.Value ~= N
        y = double(S(1:N));
            y = y - mean(y);
    [sspsd_y,fss] = Laser.SSPSD(y,Ts);
        
       N = N_Samples_EditField.Value;
       in_lag = N/10;
        out_lag = N;
        window_size = N;
        step_length = N;

       
    end
    
        % for spectra calculations
    
    
    if strcmp(Conv_Switch.Value,'On')
        Track_Switch.Enable = 1;
        Theoreticalspectra_CheckBox.Enable = 1;
        
        % conventional phase and amplitude extraction
        ya= hilbert(y);
        a = sign(randn(N,1)); % random sequence 
        ya_mod = a.*ya; % modulated sequence
        
        
        conv_total_phase = unwrap(angle(ya));
        conv_amplitude = abs(ya);
        mean_conv_amplitude = mean(conv_amplitude);
        conv_phase_noise = conv_total_phase - 2*pi*beat_f*t;
        conv_amplitude_noise= conv_amplitude - mean_conv_amplitude;
        
        y_tone_conv = conv_amplitude.*exp(1i.*(conv_phase_noise));
        
        % conv BER calculation
        y_demod_conv = ya_mod.*exp(-1i.*(conv_total_phase));
        a_conv = sign(real(y_demod_conv));
        errors_conv = sum(a_conv ~= a);
        ber_conv= errors_conv/N;
        disp(['Conventional BER : ',num2str(ber_conv)])
        
        % RIN, FN and PN calculation
        [psd_conv_phase_noise,f_pn] = Laser.SSPSD(conv_phase_noise(in_lag:out_lag),Ts,window_size,step_length);
        [rin_conv, f_rin] = Laser.SSPSD(detrend(conv_amplitude(in_lag:out_lag).^2),Ts,window_size,step_length);
        rin_conv = rin_conv./(mean(conv_amplitude(in_lag:out_lag).^2));
        
        
        [conv_autocorr_amp,lags] = xcorr(conv_amplitude_noise(in_lag:out_lag),conv_amplitude_noise(in_lag:out_lag),'normalized');
        conv_autocorr_phase = xcorr(conv_phase_noise(in_lag:out_lag),conv_phase_noise(in_lag:out_lag),'normalized');
        
    end
    
    if strcmp(Track_Switch.Value,'On')
        Smooth_Switch.Enable = 1;
        if strcmp(Smooth_Switch.Value,'On') % tracking + smoothing
            EM_param_track.filtering_only = 0;
            Expectationmaximization_CheckBox.Enable = 1;
            [m0s,P0s,Qs,Rs,~,~,m,~,yf_ap,ET,ms,~,ys_ap] = Filter.EM_OFC_ap(m0,P0,y(1:N),theta_EM,EM_param_track);
            if Expectationmaximization_CheckBox.Value
                    theta_EM.Q = Qs;
                    theta_EM.R = Rs;
                    m0 = m0s;
                    P0 = P0s;
                 if strcmp(Usesliders_CheckBox.Value,'Sliders')
                    R_Slider.Value = theta_EM.R;
                    Qa_Slider.Value = theta_EM.Q(1,1);
                    Qp_Slider.Value = theta_EM.Q(2,2);
                    Qap_Slider.Value = theta_EM.Q(1,2);
                 else
                    R_curr.Value = theta_EM.R;
                    Qa_curr.Value = theta_EM.Q(1,1);
                    Qp_curr.Value = theta_EM.Q(2,2);
                    Qap_curr.Value = theta_EM.Q(1,2);   
                 end
                 
            end
            
            if RefineA_CheckBox.Value
                theta_EM.a = theta_EM.a + mean(ms(:,1));
            end
            
            if Refineomega_CheckBox.Value
                PP = polyfit(t,ms(:,2),1);
                theta_EM.omega = theta_EM.omega + PP(1);
            end
            smooth_amplitude = theta_EM.a + ms(:,1);
            smooth_phase_noise_ap = ms(:,2);
            smooth_amplitude_noise = smooth_amplitude - mean(smooth_amplitude);
            
            % calculate smoothed spectra
            sspsd_ys_ap = Laser.SSPSD(ys_ap,Ts,window_size,step_length);
            [psd_smooth_phase_noise_ap,f_pn] = Laser.SSPSD(smooth_phase_noise_ap(in_lag:out_lag),Ts,window_size,step_length);
            [rin_smooth, f_rin]  = Laser.SSPSD(detrend(smooth_amplitude(in_lag:out_lag).^2),Ts,window_size,step_length);
            rin_smooth = rin_smooth./mean(smooth_amplitude(in_lag:out_lag).^2);
            
            %calculate autocorrelations
            [smooth_autocorr_amp,lags] = xcorr(smooth_amplitude_noise(in_lag:out_lag),smooth_amplitude_noise(in_lag:out_lag),'normalized');
            smooth_autocorr_phase_ap = xcorr(smooth_phase_noise_ap(in_lag:out_lag),smooth_phase_noise_ap(in_lag:out_lag),'normalized');
            
            % BER estimation
             y_demod_smooth_ap = ya_mod.*exp(-1i.*(2*pi*beat_f.*t + smooth_phase_noise_ap)); % correct for phase noise term;
                a_smooth_ap = sign(real(y_demod_smooth_ap));
                errors_smooth_ap = sum(a_smooth_ap ~= a);
                ber_smooth_ap = errors_smooth_ap/N;  
                disp(['Smoothing BER : ',num2str(ber_smooth_ap)])
            
        else % tracking only
            Expectationmaximization_CheckBox.Enable = 0;
            EM_param_track.filtering_only = 1;
            [~,~,~,~,~,~,m,~,yf_ap,ET,~,~,~] = Filter.EM_OFC_ap(m0,P0,y(1:N),theta_EM,EM_param_track);
            
        end
        
        filt_amplitude = theta_EM.a + m(:,1);
        filt_phase_noise_ap = m(:,2);
        filt_amplitude_noise = filt_amplitude - mean(filt_amplitude);
        
        % calculate filtered spectra
        sspsd_y_ap = Laser.SSPSD(yf_ap,Ts,window_size,step_length);
        [psd_filt_phase_noise_ap,f_pn] = Laser.SSPSD(filt_phase_noise_ap(in_lag:out_lag),Ts,window_size,step_length);
        [rin_filt, f_rin]  = Laser.SSPSD(detrend(filt_amplitude(in_lag:out_lag).^2),Ts,window_size,step_length);
        rin_filt = rin_filt./mean(filt_amplitude(in_lag:out_lag).^2);
        
        %calculate autocorrelations
        [filt_autocorr_amp, lags] = xcorr(filt_amplitude_noise(in_lag:out_lag),filt_amplitude_noise(in_lag:out_lag),'normalized');
        filt_autocorr_phase_ap = xcorr(filt_phase_noise_ap(in_lag:out_lag),filt_phase_noise_ap(in_lag:out_lag),'normalized');
        
        % BER comparison - parameter validation
        y_demod_filt_ap = ya_mod.*exp(-1i.*(2*pi*beat_f.*t + filt_phase_noise_ap)); % correct for phase noise term;
        a_filt_ap = sign(real(y_demod_filt_ap));
        errors_filt_ap = sum(a_filt_ap ~= a);
        ber_filt_ap = errors_filt_ap/N;  
        disp(['Filtering BER : ',num2str(ber_filt_ap)])
        E.Text = num2str(sum(ET));
        dE = sum(ET)-prev_E;
        DE.Text = num2str(dE);
        prev_E = sum(ET);
    else
        Smooth_Switch.Enable = 0;
    end
    
    % run EKF (always)
    
    
    
    if strcmp(Simulate_Switch.Value,'On')
        
        q = mvnrnd(zeros(N,2),theta_EM.Q);
        sim_phase_noise = cumsum(q(:,2));
        sim_amplitude_noise = cumsum(q(:,1));
        sim_amplitude = sim_amplitude_noise + theta_EM.a;
        
        y_sim = sim_amplitude.*cos(theta_EM.omega.*t + sim_phase_noise);
        sspsd_y_sim = Laser.SSPSD(y_sim,Ts,window_size,step_length);
        
        % calculate spectra
        [psd_sim_phase_noise_ap,f_pn] = Laser.SSPSD(sim_phase_noise(in_lag:out_lag),Ts,window_size,step_length);
        [rin_sim, f_rin]  = Laser.SSPSD(detrend(sim_amplitude(in_lag:out_lag).^2),Ts,window_size,step_length);
        rin_sim = rin_sim./mean(sim_amplitude(in_lag:out_lag).^2);
        
        %calculate autocorrelations
        [sim_autocorr_amp,lags] = xcorr(sim_amplitude_noise(in_lag:out_lag),sim_amplitude_noise(in_lag:out_lag),'normalized');
        sim_autocorr_phase_ap = xcorr(sim_phase_noise(in_lag:out_lag),sim_phase_noise(in_lag:out_lag),'normalized');
        
    end
    
    if Theoreticalspectra_CheckBox.Value
       LW_est = theta_EM.Q(2,2)/Ts/2/pi;
       LW.Text = num2str(LW_est);
       psd_lorentian_phase_noise = LW_est/pi./(f_pn.^2 ); 
       
       lorentzian_profile = (0.5*LW_est)./((fss - theta_EM.omega/2/pi).^2 + (0.5*LW_est)^2)/pi/Fs;
    end
    
    % slider priority
    if strcmp(Usesliders_CheckBox.Value,'Sliders')
           
        R_Slider.Enable = 1;
        Qp_Slider.Enable = 1;
        Qa_Slider.Enable = 1;
        Qap_Slider.Enable = 1;
        omega_Slider.Enable = 1;
        A_Slider.Enable = 1;
        dt_Slider.Enable = 1;
        
        R_min.Enable = 1;
        R_max.Enable = 1;
        Qp_min.Enable = 1;
        Qp_max.Enable = 1;
        Qa_min.Enable = 1;
        Qa_max.Enable = 1;
        Qap_min.Enable = 1;
        Qap_max.Enable = 1;
        omega_min.Enable = 1;
        omega_max.Enable = 1;
        A_min.Enable = 1;
        A_max.Enable = 1;
        dt_min.Enable = 1;
        dt_max.Enable = 1;
        
        R_curr.Enable = 0;
        Qp_curr.Enable = 0;
        Qa_curr.Enable = 0;
        Qap_curr.Enable = 0;
        omega_curr.Enable = 0;
        A_curr.Enable = 0;
        dt_curr.Enable = 0;
        
        R_Slider.Limits(1) = R_min.Value;
        R_Slider.Limits(2) = R_max.Value;
        R_curr.Value = R_Slider.Value;
        R_curr.Step = (R_max.Value - R_min.Value)/1e3;
        theta_EM.R = R_Slider.Value;
        
        Qp_Slider.Limits(1) = Qp_min.Value;
        Qp_Slider.Limits(2) = Qp_max.Value;
        Qp_curr.Value = Qp_Slider.Value;
        Qp_curr.Step = (Qp_max.Value - Qp_min.Value)/1e3;
        theta_EM.Q(2,2) = Qp_Slider.Value;
        
        Qa_Slider.Limits(1) = Qa_min.Value;
        Qa_Slider.Limits(2) = Qa_max.Value;
        Qa_curr.Value = Qa_Slider.Value;
        Qa_curr.Step = (Qa_max.Value - Qa_min.Value)/1e3;
        theta_EM.Q(1,1) = Qa_Slider.Value;
        
        Qap_Slider.Limits(1) = Qap_min.Value;
        Qap_Slider.Limits(2) = Qap_max.Value;
        Qap_curr.Value = Qap_Slider.Value;
        Qap_curr.Step = (Qap_max.Value - Qap_min.Value)/1e3;
        theta_EM.Q(1,2) = Qap_Slider.Value;
        theta_EM.Q(2,1) = Qap_Slider.Value;
        
        omega_Slider.Limits(1) = omega_min.Value;
        omega_Slider.Limits(2) = omega_max.Value;
        omega_curr.Value = omega_Slider.Value;
        omega_curr.Step = (omega_max.Value - omega_min.Value)/1e3;
        theta_EM.omega = omega_Slider.Value;
        
        A_Slider.Limits(1) = A_min.Value;
        A_Slider.Limits(2) = A_max.Value;
        A_curr.Value = A_Slider.Value;
        A_curr.Step = (A_max.Value - A_min.Value)/1e3;
        theta_EM.a = A_Slider.Value;
        
        dt_Slider.Limits(1) = dt_min.Value;
        dt_Slider.Limits(2) = dt_max.Value;
        dt_curr.Value = dt_Slider.Value;
        dt_curr.Step = (dt_max.Value - dt_min.Value)/1e3;
        theta_EM.dt = dt_Slider.Value;
        
        
    else
        
        R_Slider.Enable = 0;
        Qp_Slider.Enable = 0;
        Qa_Slider.Enable = 0;
        Qap_Slider.Enable = 0;
        omega_Slider.Enable = 0;
        A_Slider.Enable = 0;
        dt_Slider.Enable = 0;
        
        R_min.Enable = 0;
        R_max.Enable = 0;
        Qp_min.Enable = 0;
        Qp_max.Enable = 0;
        Qa_min.Enable = 0;
        Qa_max.Enable = 0;
        Qap_min.Enable = 0;
        Qap_max.Enable = 0;
        omega_min.Enable = 0;
        omega_max.Enable = 0;
        A_min.Enable = 0;
        A_max.Enable = 0;
        dt_min.Enable = 0;
        dt_max.Enable = 0;
        
        R_curr.Enable = 1;
        Qp_curr.Enable = 1;
        Qa_curr.Enable = 1;
        Qap_curr.Enable = 1;
        omega_curr.Enable = 1;
        A_curr.Enable = 1;
        dt_curr.Enable = 1;
        
        if  R_curr.Value >= theta_EM.R % increase
            
            theta_EM.R = R_curr.Value;
            R_Slider.Limits(2) = R_curr.Value*2;
            R_Slider.Value = R_curr.Value;
            R_Slider.Limits(1) = R_curr.Value*0.1;
            R_max.Value = R_Slider.Limits(2);
            R_min.Value = R_Slider.Limits(1);
            
        else % decrease
            
            theta_EM.R = R_curr.Value;
            R_Slider.Limits(1) = R_curr.Value*0.1;
            R_Slider.Value = R_curr.Value;
            R_Slider.Limits(2) = R_curr.Value*2;
            R_max.Value = R_Slider.Limits(2);
            R_min.Value = R_Slider.Limits(1);
            
        end
        R_curr.Step = (R_max.Value - R_min.Value)/1e1;
        
        if  Qp_curr.Value >= theta_EM.Q(2,2) % increase
            
            theta_EM.Q(2,2) = Qp_curr.Value;
            Qp_Slider.Limits(2) = Qp_curr.Value*2;
            Qp_Slider.Value = Qp_curr.Value;
            Qp_Slider.Limits(1) = Qp_curr.Value*0.1;
            Qp_max.Value = Qp_Slider.Limits(2);
            Qp_min.Value = Qp_Slider.Limits(1);
            
        else % decrease
            theta_EM.Q(2,2) = Qp_curr.Value;
            Qp_Slider.Limits(1) = Qp_curr.Value*0.1;
            Qp_Slider.Value = Qp_curr.Value;
            Qp_Slider.Limits(2) = Qp_curr.Value*2;
            Qp_max.Value = Qp_Slider.Limits(2);
            Qp_min.Value = Qp_Slider.Limits(1);
            
        end
        Qp_curr.Step = (Qp_max.Value - Qp_min.Value)/1e1;
        
        
        if  Qa_curr.Value >= theta_EM.Q(1,1) % increase
            theta_EM.Q(1,1) = Qa_curr.Value;
            Qa_Slider.Limits(1) = Qa_Slider.Value*0.1;
            Qa_Slider.Limits(2) = Qa_Slider.Value*2;
            Qa_Slider.Value = Qa_curr.Value;
            Qa_max.Value = Qa_Slider.Limits(2);
            Qa_min.Value = Qa_Slider.Limits(1);
         else % decrease
             theta_EM.Q(1,1) = Qa_curr.Value;
            Qa_Slider.Limits(2) = Qa_Slider.Value*2; 
            Qa_Slider.Limits(1) = Qa_Slider.Value*0.1;
            Qa_Slider.Value = Qa_curr.Value;
            Qa_max.Value = Qa_Slider.Limits(2);
            Qa_min.Value = Qa_Slider.Limits(1);
        end     
        Qa_curr.Step = (Qa_max.Value - Qa_min.Value)/1e1;
        
        
         if  Qap_curr.Value >= theta_EM.Q(1,2) % increase
            theta_EM.Q(1,2) = Qap_curr.Value;
            theta_EM.Q(2,1) = Qap_curr.Value;
            Qap_Slider.Limits(1) = -abs(Qap_curr.Value)*2;
            Qap_Slider.Limits(2) = abs(Qap_curr.Value)*2+eps;
            Qap_Slider.Value = Qap_curr.Value;
            Qap_max.Value = Qap_Slider.Limits(2);
            Qap_min.Value = Qap_Slider.Limits(1);
         else % decrease
            theta_EM.Q(1,2) = Qap_curr.Value;
            theta_EM.Q(2,1) = Qap_curr.Value;
            Qap_Slider.Limits(2) = abs(Qap_curr.Value)*2+eps;
            Qap_Slider.Limits(1) = -abs(Qap_curr.Value)*2;
            Qap_Slider.Value = Qap_curr.Value;  
            Qap_max.Value = Qap_Slider.Limits(2);
            Qap_min.Value = Qap_Slider.Limits(1);
             
        end     
        Qap_curr.Step = (Qap_max.Value - Qap_min.Value)/1e1;
        
        if  omega_curr.Value >= theta_EM.omega % increase
             theta_EM.omega = omega_curr.Value;
            omega_Slider.Limits(2) = omega_curr.Value*2;
            omega_Slider.Limits(1) = omega_curr.Value*0.1;
            omega_Slider.Value = omega_curr.Value;
            omega_max.Value = omega_Slider.Limits(2);
            omega_min.Value = omega_Slider.Limits(1);
        else % decrease
             theta_EM.omega = omega_curr.Value;
            omega_Slider.Limits(1) = omega_curr.Value*0.1;
            omega_Slider.Limits(2) = omega_curr.Value*2;
            omega_Slider.Value = omega_curr.Value;
            omega_max.Value = omega_Slider.Limits(2);
            omega_min.Value = omega_Slider.Limits(1);
            
        end 
        omega_curr.Step = (omega_max.Value - omega_min.Value)/1e1;
        
        if  A_curr.Value >= theta_EM.a % increase
            theta_EM.a = A_curr.Value;
            A_Slider.Limits(1) = A_Slider.Value*0.1;
            A_Slider.Limits(2) = A_Slider.Value*2; 
            A_Slider.Value = A_curr.Value;
            A_max.Value = A_Slider.Limits(2);
            A_min.Value = A_Slider.Limits(1);
        else % decrease
            theta_EM.a = A_curr.Value;
            A_Slider.Limits(2) = A_Slider.Value*2; 
            A_Slider.Limits(1) = A_Slider.Value*0.1;
            A_Slider.Value = A_curr.Value;
            A_max.Value = A_Slider.Limits(2);
            A_min.Value = A_Slider.Limits(1);
        end     
        A_curr.Step = (A_max.Value - A_min.Value)/1e1;
        
        
        if  dt_curr.Value >= theta_EM.dt % increase
            theta_EM.dt = dt_curr.Value;
            dt_Slider.Limits(1) = dt_Slider.Value*0.1;
            dt_Slider.Limits(2) = dt_Slider.Value*2;
            dt_Slider.Value = dt_curr.Value;
            dt_max.Value = dt_Slider.Limits(2);
            dt_min.Value = dt_Slider.Limits(1);
        else % decrease
            theta_EM.dt = dt_curr.Value;
            dt_Slider.Limits(2) = dt_Slider.Value*2;
            dt_Slider.Limits(1) = dt_Slider.Value*0.1;
            dt_Slider.Value = dt_curr.Value;
            dt_max.Value = dt_Slider.Limits(2);
            dt_min.Value = dt_Slider.Limits(1);
        end  
        dt_curr.Step = (dt_max.Value - dt_min.Value)/1e1;
    end
    
    if Reset_Button.Value
       Usesliders_CheckBox.Value = 'Spinners';
       theta_EM = saved_param;
       R_curr.Value = theta_EM.R;
      Qp_curr.Value =  theta_EM.Q(2,2);
       Qa_curr.Value =  theta_EM.Q(1,1);
        Qap_curr.Value = theta_EM.Q(1,2);
         omega_curr.Value = theta_EM.omega;
        A_curr.Value = theta_EM.a ;
        dt_curr.Value = theta_EM.dt;
        
        Reset_Button.Value = 0;
    end
    
    
    
    
    % plotting
    cla(Signal_time)
    hold(Signal_time,'on')
    if strcmp(Conv_Switch.Value,'On')
    plot(Signal_time,t,y,'k');
    end
    if strcmp(Track_Switch.Value,'On')
    plot(Signal_time,t,yf_ap,'b');    
    end
    if strcmp(Smooth_Switch.Value,'On')
    plot(Signal_time,t,ys_ap,'r');    
    end
    if strcmp(Simulate_Switch.Value,'On')
    plot(Signal_time,t,y_sim,'m');
    end
    hold(Signal_time,'off')
    
    
    cla(Signal_frequency)
    hold(Signal_frequency,'on')
    if strcmp(Conv_Switch.Value,'On')
    plot(Signal_frequency,fss,10*log10(sspsd_y),'k');
    end
    if strcmp(Track_Switch.Value,'On')
    plot(Signal_frequency,fss,10*log10(sspsd_y_ap),'b');
    end
    if strcmp(Smooth_Switch.Value,'On')
    plot(Signal_frequency,fss,10*log10(sspsd_ys_ap),'r');
    end
    if strcmp(Simulate_Switch.Value,'On')
    plot(Signal_frequency,fss,10*log10(sspsd_y_sim),'m');
    end
    if Theoreticalspectra_CheckBox.Value
    plot(Signal_frequency,fss,10*log10(ones(length(fss),1).*2*Ts*theta_EM.R),'c'); 
    end
    hold(Signal_frequency,'off')
    
    cla(Amplitude_time)
    hold(Amplitude_time,'on')
    if strcmp(Conv_Switch.Value,'On')
    plot(Amplitude_time,t,conv_amplitude_noise,'k');
    end
    if strcmp(Track_Switch.Value,'On')
    plot(Amplitude_time,t,filt_amplitude_noise,'b');
    end
    if strcmp(Smooth_Switch.Value,'On')
    plot(Amplitude_time,t,smooth_amplitude_noise,'r'); 
    end
    hold(Amplitude_time,'off')
    
    
    cla(Phase_time)
    hold(Phase_time,'on')
    if strcmp(Conv_Switch.Value,'On')
    plot(Phase_time,t,conv_phase_noise,'k');
    end
    if strcmp(Track_Switch.Value,'On')
    plot(Phase_time,t,filt_phase_noise_ap,'b');
    end
    if strcmp(Smooth_Switch.Value,'On')
    plot(Phase_time,t,smooth_phase_noise_ap,'r');
    end
    hold(Phase_time,'off')
    
    
    cla(RIN)
    hold(RIN,'on')
    if strcmp(Conv_Switch.Value,'On')
    semilogx(RIN,f_rin,10*log10(rin_conv),'k');
    end
    if strcmp(Track_Switch.Value,'On')
    semilogx(RIN,f_rin,10*log10(rin_filt),'b');
    end
    if strcmp(Smooth_Switch.Value,'On')
    semilogx(RIN,f_rin,10*log10(rin_smooth),'r');
    end
    if strcmp(Simulate_Switch.Value,'On')
    semilogx(RIN,f_rin,10*log10(rin_sim),'m');
    end
    semilogx(RIN,freq_rin_ref,PSD_rin_ref,'g');
    hold(RIN,'off')
    

    cla(PN)
    hold(PN,'on')
    if strcmp(Conv_Switch.Value,'On')
    semilogx(PN,f_pn,10*log10(psd_conv_phase_noise),'k');
    end
    if strcmp(Track_Switch.Value,'On')
    semilogx(PN,f_pn,10*log10(psd_filt_phase_noise_ap),'b');
    end
    if strcmp(Smooth_Switch.Value,'On')
    semilogx(PN,f_pn,10*log10(psd_smooth_phase_noise_ap),'r');
    end
    if strcmp(Simulate_Switch.Value,'On')
    semilogx(PN,f_pn,10*log10(psd_sim_phase_noise_ap),'m');
    end
    if Theoreticalspectra_CheckBox.Value
       semilogx(PN,f_pn,10*log10(psd_lorentian_phase_noise),'c');
    end
    semilogx(PN,freq_pn_ref,PSD_pn_ref,'g');
    hold(PN,'off')
    
    
    cla(amp_autocorr)
    hold(amp_autocorr,'on')
    if strcmp(Conv_Switch.Value,'On')
    plot(amp_autocorr,lags.*Ts,conv_autocorr_amp,'k');
    end
    if strcmp(Track_Switch.Value,'On')
    plot(amp_autocorr,lags.*Ts,filt_autocorr_amp,'b');
    end
    if strcmp(Smooth_Switch.Value,'On')
    plot(amp_autocorr,lags.*Ts,smooth_autocorr_amp,'r');
    end
    if strcmp(Simulate_Switch.Value,'On')
    plot(amp_autocorr,lags.*Ts,sim_autocorr_amp,'m');
    end
    hold(amp_autocorr,'off')
  
    
    cla(phase_autocorr)
    hold(phase_autocorr,'on')
    if strcmp(Conv_Switch.Value,'On')
    plot(phase_autocorr,lags.*Ts,conv_autocorr_phase,'k');
    end
    if strcmp(Track_Switch.Value,'On')
    plot(phase_autocorr,lags.*Ts,filt_autocorr_phase_ap,'b');
    end
    if strcmp(Smooth_Switch.Value,'On')
    plot(phase_autocorr,lags.*Ts,filt_autocorr_phase_ap,'r');
    end
    if strcmp(Simulate_Switch.Value,'On')
    plot(phase_autocorr,lags.*Ts,sim_autocorr_phase_ap,'m');
    end
    hold(phase_autocorr,'off')
    

    drawnow
    disp('Refresh')
    
    exit = Exit_Button.Value;
    
end
save('starting_parameters','theta_EM')
