clc
clear all
close all
warning('off','all')

delete DistM
delete DiffMX
delete DiffMY

    
    %Funcao---> 1. SANDWICH e FGM; 3. VIBRACOES LIVRES
   
    funcao=6; save funcao
    
    [A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,k,E,h,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,N,A,B,D]=Points_BAT;
            
    [Q,bestvar,bestvarG,BESTvar,fitness,sol,minimo,maximo,Ngen,minG,maxG,optim,BEST,ng,A,r]=BAT_algorithm(A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,k,E,h,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,N,A,B,D);
            
            if size(bestvarG,2)==1
                if funcao==1 || funcao==2
                MASSAF=MASSA0(BESTvar(end,1));
                fprintf('MASSA: %6.4f\n',MASSAF);
                end
            elseif size(bestvarG,2)==2
            MASSAF=MASSA0(BESTvar(end,1),BESTvar(end,2));
            fprintf('MASSA: %6.4f\n',MASSAF);
            end           
           
           fprintf('Valor optimo: %E\n',BEST);
           fprintf('Encontrado na iteracao nr: %g\n',ng);

