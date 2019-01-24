function windata = flat_hanning_win(taxis,data,winbgt,winendt,tapertime)
% windata = flat_hanning_win(taxis,data,winbgt,winendt,tapertime)
% function to apply a box-car window with a hanning taper end.
% 
% Written by Ge Jin
% Adapted by Z. Eilon 04/2016

if size(data,2)>size(data,1)
    data = data';
end

windata = data;

outind = (taxis < winbgt | taxis > winendt);
windata(outind,:) = 0;

inind = find(taxis >= winbgt & taxis <= winendt);

dt = taxis(2)-taxis(1);
taperlength = floor(tapertime./dt);

tapr = hanning(taperlength*2);
windata(inind(1:taperlength),:) = diag(tapr(1:taperlength))*windata(inind(1:taperlength),:);
windata(inind(end-taperlength+1:end),:) = diag(tapr(taperlength+1:end))*windata(inind(end-taperlength+1:end),:);

end
