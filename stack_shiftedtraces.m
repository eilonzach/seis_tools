function stack = stack_shiftedtraces(data,Nsamp_shift)
% stack = stack_shiftedtraces(data,Nsamp_shift)
% function to take a matrix of multiple data series and shift them each by
% an some amount of samples (e.g. from a xcorr) before stacking over them
% and producing a final time series with the same length as the original,
% padding with zeros where necessary.
%  INPUTS:
%   data         NxM matrix of M time series each length N samples
%   Nsamp_shift  Mx1 vector with number of samples by which to shift each
%                vector of data. Positive values shift their traces
%                backwards in time - i.e. if a phase is LATE, we apply a
%                POSITIVE shift bringing it earlier to line up w/ the stack
% 
%  OUTPUTS:
%   stack        Nx1 vector of the final, stacked data;

[N,M] = size(data);