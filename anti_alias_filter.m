function outdata = anti_alias_filter(indata,indelta,outdelta)

if indelta == outdelta
    outdata = indata; return
end

	N = length(indata);
	dt = indelta;
	
	corner_f = 1/2/outdelta;
	fN = 1./2/dt;
	w = corner_f/fN;

	[b a] = butter(2,w,'low');
	outdata = filtfilt(b,a,indata);

end