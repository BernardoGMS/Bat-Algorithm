function [A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,h,E]=DRBF_VIB_BAT

propriedades=1;

% Construcao de K:

if propriedades==1
    
[A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,h,E]=sandwichFGM_BAT;

elseif propriedades==2
    
% [A,B,D,I0,I1,I2,h,E]=FGM_Mix;
[A,B,D,I0,I1,I2,h,E]=FGM_Mori;

end


end