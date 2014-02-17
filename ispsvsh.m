function [ psvsh ] = ispsvsh( phase_name )
%[ psvsh ] = ispsvsh( phase_name )
%  quick and dirty function to parse name of phase and work out if the
%  radiation pattern that's relevant is:
% p  = 1  -- for phases beginning with P
% sv = 2  -- for phases beginning with S that then undergo transition to P
% sh = 3  -- for phases with only S(or ScS, Sdiff, etc)

if strncmpi(phase_name,'P',1)
psvsh = 1;
else
    if isempty(strfind(phase_name,'P')) && isempty(strfind(phase_name,'K'))
        psvsh = 3;
    else
        psvsh = 2;
    end
end

end

