function [Q,bestvar,bestvarG,fitness,sol,minimo,maximo,minG,maxG]=...
    inicializa(n,Ngen,alpha,gama,optim,Qmin,Qmax,epsilon,d,Fperi,Fperf,adaptive,hibrid,Lb,Ub,A,r,A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,EE,hh,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN,MA,MB,MD)

Q=zeros(1,n);
bestvar=zeros(Ngen,d);bestvarG=zeros(Ngen,d);
fitness=zeros(1,n); sol=zeros(n,d);
SWobj=zeros(n,1);

for i=1:n
    for j=1:d
        aleat=rand;
        
        sol(i,j)=Lb(1,j)+(Ub(1,j)-Lb(1,j))*aleat;
        
    end
    
    load funcao
    
    MASSAeval=MASSA0(sol(i,1),sol(i,2));
%     MASSAeval=MASSA0(sol(i,1));
    while MASSAeval>inf
        
        for j=1:d
        sol(i,j)=Lb(1,j)+(Ub(1,j)-Lb(1,j))*rand;
        end
       MASSAeval=MASSA0(sol(i,1),sol(i,2));
%        MASSAeval=MASSA0(sol(i,1));
%        disp('PASSOU')
    end

    if funcao==1
    
%         MINIMIZAR DEFORMADA FGM
      SW=FGM_ss(A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,kkk,sol(i,:),AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN);
    
     SWobj(i,1)=SW;
        
    
     elseif funcao==3
        
            %MAXIMIZAR 1a FRQ SANDWICH FGM
     MODO=FGM_ss_Vib(A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,sol(i,:),AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN);
            
     SWobj(i,1)=MODO;
     
    else
            
        % DIVERSAS FUNCOES
%         SWobj(i,1)=(sol(i,1)^3-sol(i,2))^2+2*(sol(i,2)-sol(i,1))^4;
        SWobj(i,1)=2*sol(i,1)^2+sol(i,1)*sol(i,2)+sol(i,2)^2+sol(i,2)*sol(i,3)+sol(i,3)^2-6*sol(i,1)-7*sol(i,2)-8*sol(i,3)+9;
%         SWobj(i,1)=abs(- 0.01299686456980716395663890921133/(((70000000.0*sol(i,1) - 2800000.0)*(625.0*sol(i,1)^2*sol(i,2)^2 + 1875.0*sol(i,1)^2*sol(i,2) + 1250.0*sol(i,1)^2 + ...
%             50.0*sol(i,1)*sol(i,2) + 50.0*sol(i,1) + 2.0))/(sol(i,2)^3 + 6.0*sol(i,2)^2 + 11.0*sol(i,2) + 6.0) - 14583333333.3333339691162109375*sol(i,1)^3 + 1435897.4358974359929561614990234));
        
    end
    
end

for i= 1:n
    fitness(1,i)= fitness(1,i)+SWobj(i,1);
end

minimo=zeros(1,Ngen);maximo=zeros(1,Ngen);minG=zeros(1,Ngen);maxG=zeros(1,Ngen);
minimo(1,1)=fitness(1,1);maximo(1,1)=fitness(1,1);

for i=2:n
    if strcmp(optim,'m') && fitness(1,i)<=fitness(1,i-1)
        minimo(1,1)=fitness(1,i);
        bestvar(1,1:d)=sol(i,1:d);
    elseif strcmp(optim,'M') && fitness(1,i)>=fitness(1,i-1)
        maximo(1,1)=fitness(1,i);
        bestvar(1,1:d)=sol(i,1:d);
    end
end

minG(1,1)=minimo(1,1);
maxG(1,1)=maximo(1,1);
bestvarG(1,1:d)=bestvar(1,1:d);

end