function simil = waveform_similarity(tr1,tr2,trref)
% simil = waveform_similarity(tr1,tr2,trref)
%
% Function to calculate similarity between two waveforms by reference to a
% reference waveform. This is drawn from Ritsema et al., 2002 (section 2.1)
% who sought to compare P and PP (or S/SS) waveform amplitudes with
% reference to a synthetic trace. 
%
% We assume for this analysis that the traces are already optimally aligned
% in time and are of the same length; 
%
% ZE 01/2026

if nargin > 2 &&  ~isequal(tr2(:),trref(:)) 
    % case intended in the paper, where two traces are being compared to
    % each other via a reference


    % first calculate all correlation functions
    % autocor
    psi11 = dot(tr1,tr1);
    psi22 = dot(tr2,tr2);
    psirr = dot(trref,trref);
    % xcorr
    psi1r = dot(tr1,trref);
    psi2r = dot(tr2,trref);
    % psi12 = dot(tr1,tr2); never used
    
    % Now calculate the amplitude ratios 
    A1_1 = psi1r./psirr;
    A2_1 = psi11./psi1r;
    
    A1_2 = psi2r./psirr;
    A2_2 = psi22./psi2r;
    
    % take ratios of ratios:
    R1 = min([A1_2,A2_2])./max([A1_1,A2_1]);
    R2 = max([A1_2,A2_2])./min([A1_1,A2_1]);
    
    simil = abs(R2 - R1)./(R1+R2);

else  
    % modified to just compare 1 trace to a reference, basically equate tr2
    % with trref

    % first calculate all correlation functions
    % autocorr
    psi11 = dot(tr1,tr1);
    psirr = dot(trref,trref);
    % xcorr
    psi1r = dot(tr1,trref);
    
    % Now calculate the amplitude ratios 
    A1_1 = psi1r./psirr;
    A2_1 = psi11./psi1r;

    % take ratios of ratios:
    R1 = 1./max([A1_1,A2_1]);
    R2 = 1./min([A1_1,A2_1]);
    
    simil = abs(R2 - R1)./(R1+R2);

end