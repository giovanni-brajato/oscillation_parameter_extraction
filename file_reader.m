function [X] = file_reader(param)
%FILE_READER read the file containing the trace
%   param.file_type: define if the file has a .mat extension or .h5
%   extension
%   param.file_name: name of the trace file (including the full path)
%   param.L_trace: number of samples we take from the trace

%  X: the trace signal, without the DC component
switch param.file_type
    case 'mat'
        load(param.file_name);
        X = signal;
    case 'h5'
        X = h5read(param.file_name,'/signal_trace');
end
X = X(1:param.L_trace);
X = X - mean(X);
end

