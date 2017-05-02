function [dof,tim]=scdofcalc(xx,option)
% [dof,tim]=scdofcalc(xx,option)
% scdofcalc.m  calculate degrees of freedom a la Silver and Chan for a time
% series xx or from a noise series XY filename
% optional arguments: 
% [] = xx is the time series. 
% 'filenam' = xx is the name of the file with the time series
% 
%  Z. Eilon, July 2012

if nargin==2
    if strcmp(option,'filenam')==1
infile=xx;
timeseries=load(infile);  
    end
else
timeseries=xx;
end

% clf;
f=fft(timeseries(:,1)); %NB was f=fft(timeseries(:,2));
f2=f.*conj(f);
nyq=floor(length(f2)/2);
f2 = f2(1:nyq);
f4 = f2 .* f2;
ee = sum(f2) - 0.5*(f2(1) + f2(nyq));
e4 = sum(f4) - 0.5*(f4(1) + f4(nyq));

dof =  2 * (2*ee*ee/e4 - 1)	;	% SC equation A12
tim=(timeseries(length(timeseries),1) - timeseries(1,1))/dof; 
% fprintf(' %f degrees of freedom for %f seconds per free parameter\n',dof,tim);
end