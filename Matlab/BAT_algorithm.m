function [Q,bestvar,bestvarG,BESTvar,fitness,sol,minimo,maximo,Ngen,minG,maxG,optim,BEST,ng,A,r]=BAT_algorithm(A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,EE,hh,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN,MA,MB,MD)
       


% No caso de haver BAT hidrido DE

Fperi=1; %continua com este valor se adaptive==0
Fperf=.3;
adaptive=0;
hibrid=1;

% A entrada tambem pode ser interactiva

n=60; %dimensao da populacao
Ngen=60; %numero de geracoes
A0=.9; %intensidade (loudness, constante ou crescente)
r0=.4; %taxa de impulsos (pulse rate, constante ou decrescente)
alpha=.9;
gama=.9;

optim='m'; % para minimizar ou maximizar

% Gama de freq define o escalamento, pelo que a gama podera requerer
% ajustamento

Qmin=0;Qmax=3; %fre min e max
epsilon=1/1000; %factor de escala
d=3; %dimensao das variaveis de projecto

% Limites inferior ou superior dos vectores das variaveis de projecto

Lb=zeros(1,d);Ub=zeros(1,d);

%1 ---> SANDWICH ou outra funcao (colocar LB e Ub manualmente); 2 ---> CAMADA FGM; 3 ---> FUNCAO EXPLICITA 
tipo_prop=1;

if tipo_prop==1
    
     Lb(1,1)=-5;
     Ub(1,1)=3;
     
     if d==2 || d==3
    
     Lb(1,2)=-1;
     Ub(1,2)=10;
     
     end
     
     if d==3
     
     Lb(1,3)=0;
     Ub(1,3)=5;
     
     end
    
elseif tipo_prop==2
    
    Lb(1,1)=0;
    Ub(1,1)=7;
    
else
    
        for i=1:d
    
        Lb(1,i)=0;
        Ub(1,i)=2;

        end

end


% Inicializacao de arrays

A=zeros(Ngen,n);r=zeros(Ngen,n);

for i=1: n
    
    A(1,i)=A0+rand*.1;
    r(1,i)=r0+rand*.1;

end


[Q,bestvar,~,fitness,sol,minimo,maximo,~,~]=...
    inicializa(n,Ngen,alpha,gama,optim,Qmin,Qmax,epsilon,d,Fperi,Fperf,adaptive,hibrid,Lb,Ub,A,r,A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,EE,hh,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN,MA,MB,MD);

[Ngen,minG,maxG,optim,BEST,bestvarG,BESTvar,ng,A,r]=...
    BAT(n,Ngen,alpha,gama,optim,Qmin,Qmax,epsilon,d,Fperi,Fperf,adaptive,hibrid,Lb,Ub,A,r,A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,EE,hh,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN,MA,MB,MD);


end