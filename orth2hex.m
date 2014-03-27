function [ CC,A,C,F,L,N,eta ] = orth2hex( CC,axis )
%[ CC,A,C,F,L,N,eta ] = orth2hex( CC )
% convert a 6x6 stiffness tensor CC with orthorhombic symmetry to its
% hexagonal equivalent, by spinning up (averaging) about "axis"

% if axis==1
% A = CC(1,1);
% C = (3/8)*( CC(2,2) + CC(3,3) ) + CC(2,3)/4 + CC(4,4)/2;
% F = ( CC(1,2) + CC(1,3) )/2;
% L = ( CC(5,5) + CC(6,6) )/2;
% N = ( CC(2,2) + CC(3,3) )/8     - CC(2,3)/4 + CC(4,4)/2;
% 
% CC = [  C    F     F   0 0 0;
%         F    C   A-2*N 0 0 0;
%         F  A-2*N   A   0 0 0;
%         0    0     0   N 0 0;
%         0    0     0   0 N 0;
%         0    0     0   0 0 L];
% 
% elseif axis==3

if axis==1
    CC = [CC(3,3) CC(2,3) CC(1,3) 0 0 0;
          CC(2,3) CC(2,2) CC(1,2) 0 0 0;
          CC(1,3) CC(1,2) CC(1,1) 0 0 0;
          0 0 0             CC(6,6) 0 0;
          0 0 0             0 CC(5,5) 0;
          0 0 0             0 0 CC(4,4)];
end
    
A = (3/8)*( CC(1,1) + CC(2,2) ) + CC(1,2)/4 + CC(6,6)/2;
C = CC(3,3);
F = ( CC(1,3) + CC(2,3) )/2;
L = ( CC(4,4) + CC(5,5) )/2;
N = ( CC(1,1) + CC(2,2) )/8     - CC(1,2)/4 + CC(6,6)/2;

CC = [  A    A-2*N  F   0 0 0;
      A-2*N    A    F   0 0 0;
        F      F    C   0 0 0;
        0      0    0   L 0 0;
        0      0    0   0 L 0;
        0      0    0   0 0 N];

eta = F/(A-2*L);



end

