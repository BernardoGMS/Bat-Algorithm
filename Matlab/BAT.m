function [Ngen,minG,maxG,optim,BEST,bestvarG,BESTvar,ng,A,r]=...
    BAT(n,Ngen,alpha,gama,optim,Qmin,Qmax,epsilon,d,Fperi,Fperf,adaptive,hibrid,Lb,Ub,A,r,A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,EE,hh,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN,MA,MB,MD)
    
v=zeros(n,d);

[Q,bestvar,bestvarG,fitness,sol,minimo,maximo,minG,maxG]=...
    inicializa(n,Ngen,alpha,gama,optim,Qmin,Qmax,epsilon,d,Fperi,Fperf,adaptive,hibrid,Lb,Ub,A,r,A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,EE,hh,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN,MA,MB,MD);


% Ciclo as geracoes

S=zeros(n,d);vec=zeros(1,d); 

A(2,1:n)=A(1,1:n);r(2,1:n)=r(1,1:n);

for t=2:Ngen
    
    % Se BAT hibrido com DE
    if hibrid==1
        if adaptive==1
            
            deltaF=(Fperf-Fperi)/Ngen;
            Fper=Fperi+deltaF*(t-2);
        else
            
            Fper=Fperi;
            
        end
    end
    
    fitnessnew=zeros(1,n);
    
    for i=1:n
        aleat=rand;
        Q(1,i)=Qmin+(Qmax-Qmin)*aleat;
        for j=1:d
           
        v(i,j)=v(i,j)+(sol(i,j)-bestvarG(t,j))*Q(1,i);
        S(i,j)=S(i,j)+v(i,j);
            
        end
                
    end
    
    disp('Variaveis de projecto (Equacoes v e p)')
    disp(S)
    
    if hibrid==0
        
        for i=1:n
        raleat=rand;
        bat=randi(n);
        while bat==i
          bat=randi(n); 
        end

        if raleat>r(t,bat)
            
            for a1=1:1
                for a2=1:d
                    vec(a1,a2)=rand;
                end
            end
            
            S(i,1:d)=bestvarG(t-1,1:d)+vec(1,1:d)*epsilon;
            
            
        end
        end
        
        disp('Variaveis de projecto (Sem hibrido)')
        disp(S)
        
        for i=1:n
        for j=1:d
            
            if S(i,j)>Ub(1,j)
                S(i,j)=Ub(1,j);
            elseif S(i,j)<Lb(1,j)
                S(i,j)=Lb(1,j);
            end
            
        end
        end
        
%         disp(S(bat,1:d))
        
        
    elseif hibrid==1
        
        trial=zeros(n,d);
%         auxil=zeros(n,d);
        
        for i1=1:n
            for j1=1:d
            ind2=randi(n);if ind2==i1; ind2=randi(n);end
            ind3=randi(n);if ind3==i1; ind3=randi(n);end
            ind4=randi(n);if ind4==i1; ind4=randi(n);end
            ind5=randi(n);if ind5==i1; ind5=randi(n);end
            
            if n<=6
                trial(i1,j1)=bestvarG(t-1,j1)+Fper*(S(ind2,j1)-S(ind3,j1));
            else
                trial(i1,j1)=bestvarG(t-1,j1)+Fper*(S(ind2,j1)-S(ind3,j1)+S(ind4,j1)-S(ind5,j1));
            end
            
            aleat=rand;
            S(i1,j1)=aleat*bestvarG(t-1,j1)+(1-aleat)*trial(i1,j1);
            

                        
            end
        end
        
        for i=1:n
        for j=1:d
            
            if S(i,j)>Ub(1,j)
                S(i,j)=Ub(1,j);
            elseif S(i,j)<Lb(1,j)
                S(i,j)=Lb(1,j);
            end
            
        end
        end
    
    end

    disp('Variaveis de projecto (Apos devidas tecnicas e verificacao de limites)')
    disp(S)
    
    % Avaliacao da funcao para novas populacoes
    SWobj=zeros(n,1);
    for i=1:n
%         for j=1:d-1

            load funcao
            
            if funcao==1 || funcao==3
            valmassa=inf;
            else
            valmassa=inf;
            end

             MASSAeval=MASSA0(S(i,1),S(i,2));
%            MASSAeval=MASSA0(S(i,1));
        while MASSAeval>valmassa
            for j=1:d
            S(i,j)=Lb(1,j)+(Ub(1,j)-Lb(1,j))*rand;
            end
            MASSAeval=MASSA0(S(i,1),S(i,2));
%           MASSAeval=MASSA0(S(i,1));
        end
        

    
    if funcao==1
    
%         MINIMIZAR DEFORMADA FGM
      SW=FGM_ss(A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,kkk,S(i,:),AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN);
            
    
     SWobj(i,1)=SW;
     
 
     elseif funcao==3
         
         
         %MAXIMIZAR 1a FRQ SANDWICH FGM
     MODO=FGM_ss_Vib(A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,kkk,S(i,:),AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,NNN);
        
     SWobj(i,1)=MODO;
     
    else
        
        %DIVERSAS FUNCOES
%        SWobj(i,1)=(S(i,1)^3-S(i,2))^2+2*(S(i,2)-S(i,1))^4;
         SWobj(i,1)=2*S(i,1)^2+S(i,1)*S(i,2)+S(i,2)^2+S(i,2)*S(i,3)+S(i,3)^2-6*S(i,1)-7*S(i,2)-8*S(i,3)+9;
%        SWobj(i,1)=abs(- 0.01299686456980716395663890921133/(((70000000.0*S(i,1) - 2800000.0)*(625.0*S(i,1)^2*S(i,2)^2 + 1875.0*S(i,1)^2*S(i,2) + 1250.0*S(i,1)^2 + ...
%             50.0*S(i,1)*S(i,2) + 50.0*S(i,1) + 2.0))/(S(i,2)^3 + 6.0*S(i,2)^2 + 11.0*S(i,2) + 6.0) - 14583333333.3333339691162109375*S(i,1)^3 + 1435897.4358974359929561614990234));
        
    end

            fitnessnew(1,i)=fitnessnew(1,i)+SWobj(i,1);
    end
    
    disp('VALOR DA FUNCAO EM TODOS O MORCEGOS')
    disp(fitnessnew)

    % Actualizacao da solucao
    if strcmp(optim,'m')
    minG(1,t)=fitnessnew(1,1);for j=1:d;bestvarG(t,j)=S(1,j);end;
    if t==2;BEST=minG(1,t);ng=t;for j=1:d;BESTvar(1,j)=S(1,j);end;end
    elseif  strcmp(optim,'M')
    maxG(1,t)=fitnessnew(1,1);for j=1:d;bestvarG(t,j)=S(1,j);end;
    if t==2;BEST=maxG(1,t);ng=t;for j=1:d;BESTvar(1,j)=S(1,j);end;end
    end
    fprintf('Geracao:  %g - %g\n',t)
    for i=2:n
        for j=1:d
            

            fprintf('OPTIMO GLOBAL DA GERACAO: %E\n',minG(1,t))
            if strcmp(optim,'m') && fitnessnew(1,i)<=minG(1,t)
                bestvar(t,j)=S(i,j);
                minimo(1,t)=fitnessnew(1,i);
                sol(i,j)=S(i,j);
                fitness(1,i)=fitnessnew(1,i);
                bestvarG(t,j)=S(i,j);
                
                if minimo(1,t)<BEST

                    minG(1,t)=minimo(1,t);
                    bestvarG(t,j)=S(i,j);ng=t;
                    BEST=minG(1,t);BESTvar(1,1)=S(i,1);if d==2 || d==3;BESTvar(1,2)=S(i,2);end; if d==3;BESTvar(1,3)=S(i,3);end
                    fprintf('CONVERGIU: %6.4f\n',BEST)
                    fprintf('SS: %6.4f - %6.4f\n',S(i,1),S(i,2))
                    % Actualizacao dos parametros loudness e pulse rate
                    
                    if t<Ngen
                        A(t+1,i)=alpha*A(t,i);
                        r(t+1,i)=r(1,i)*(1-exp(-gama*t));
                    end
       
                else
                        
%                     minG(1,t)=fitnessnew(1,i);%minG(1,t-1);
%                     bestvarG(t,j)=bestvarG(t-1,j);
%                     BEST=minG(1,t);
                    fprintf('NAO CONVERGIU: %6.4f\n',BEST)

                    if t<Ngen
                        A(t+1,1:n)=A(t,1:n);
                        r(t+1,1:n)=r(t,1:n);
                    end
                    
                end
                
            elseif strcmp(optim,'M') && fitnessnew(1,i)>=maxG(1,t)

                bestvar(t,j)=S(i,j);
                maximo(1,t)=fitnessnew(1,i);
                sol(i,j)=S(i,j);
                fitness(1,i)=fitnessnew(1,i);
                bestvarG(t,j)=S(i,j);
                
                if maximo(1,t)>BEST
        
                    maxG(1,t)=maximo(1,t);
                    bestvarG(t,j)=S(i,j);ng=t;
                    BEST=maxG(1,t);BESTvar(1,1)=S(i,1);if d==2 || d==3;BESTvar(1,2)=S(i,2);end; if d==3;BESTvar(1,3)=S(i,3);end
                    fprintf('CONVERGIU: %6.4f\n',BEST)
                    fprintf('SS: %6.4f - %6.4f\n',S(i,1),S(i,2))
                    
                    % Actualizacao dos parametros loudness e pulse rate
                    
                    if t<=Ngen
                        A(t+1,i)=alpha*A(t,i);
                        r(t+1,i)=r(1,i)*(1-exp(-gama*t));
                    end
                        
                else
                    
%                     maxG(1,t)=maxG(1,t-1);
%                     bestvarG(t,j)=bestvarG(t-1,j);
%                     BEST=maxG(1,t-1);
                    fprintf('NAO CONVERGIU: %6.4f\n',BEST)

                    if t<Ngen
                        A(t+1,1:n)=A(t,1:n);
                        r(t+1,1:n)=r(t,1:n);
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
end


end