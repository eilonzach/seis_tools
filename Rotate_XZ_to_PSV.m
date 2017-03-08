function [P,SV] = Rotate_XZ_to_PSV(X,Z,VpSurf,VsSurf,RP)
    % [P,SV] = Rotate_XZ_to_PSV(X,Z,VpSurf,VsSurf,RP)
    % 
    % here RP is the flat Earth ray parameter, in s/km
    %         
    % Decompose vertical/radial into P/SV
    % (i.e., Kennett [1991] and Bostock [1998])
    %
    % !! NOTE !! Z in this equation is positive down !!
    % We simply make Z_comp_sac (above) negative

    a   = VpSurf;
    b   = VsSurf;

    qa  = sqrt((a^-2)-(RP^2));
    qb  = sqrt((b^-2)-(RP^2));
    A   = RP*(b^2)/a;
    B   = ((b^2)*(RP^2)-0.5)/(a*qa);
    C   = (0.5-(b^2)*(RP^2))/(b*qb);
    D   = RP*b;

    P  = A*X+B*Z;
    SV = C*X+D*Z;

end