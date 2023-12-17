function [A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,k,E,h,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,N,A,B,D]=Points_BAT


%%-------------------------------------------------------------------------
%               PLACAS - MODELACAO POR FUNCOES DE BASE RADIAL
%%-------------------------------------------------------------------------

%% Parametros dimensionais: 

% Pontos de fronteira:

LXF=1;
LYF=1;

XNDF=10;
YNDF=10;

XSTF=LXF/XNDF;
YSTF=LYF/YNDF;
nnF=0;

% Pontos internos:

LX=LXF-2*XSTF;
LY=LYF-2*YSTF;

% Criacao dos pontos de fronteira da malha:

nXY=2;
NFRT=2*(XNDF-1)+2*(YNDF-1)+4;
CF=zeros(nXY,NFRT);

for j=1:YNDF+1;nnF=nnF+1;CF(1,nnF)=-LXF./2;CF(2,nnF)=-(LYF./2)+(j-1)*YSTF;end
for j=1:YNDF+1;nnF=nnF+1;CF(1,nnF)=LXF./2;CF(2,nnF)=-(LYF./2)+(j-1)*YSTF;end
for i=1:XNDF-1;nnF=nnF+1;CF(1,nnF)=-(LXF./2)+i*XSTF;CF(2,nnF)=-LYF./2;end
for i=1:XNDF-1;nnF=nnF+1;CF(1,nnF)=-(LXF./2)+i*XSTF;CF(2,nnF)=LYF./2;end;NF=nnF;

% Criacao dos pontos interiores da malha:

xx=linspace(-LX/2,LX/2,XNDF-1);yy=linspace(-LY/2,LY/2,YNDF-1);
[xt,yt]=meshgrid(xx,yy);
CI=[xt(:) yt(:)]';N=size(CI,2);NT=N+NF;

CXY=zeros(nXY,NT);

%% Criacao das coordenadas finais de todos os pontos:

for i=1:NF;CXY(1,i)=CF(1,i);CXY(2,i)=CF(2,i);end
for j=1:N;CXY(1,NF+j)=CI(1,j);CXY(2,NF+j)=CI(2,j);end

DistM=DistanceMatrix(CXY',CXY');
DiffMX=DifferenceMatrix(CXY(1,:)',CXY(1,:)');
DiffMY=DifferenceMatrix(CXY(2,:)',CXY(2,:)');


%% Parametros de forma

q=0.5;
ep=1;

%% Matrizes de Colocacao

rbf=@(c,r)((ep.*r).^2+c^2).^q;
rbfdx=@(c,r,dx)2*dx.*(ep^2).*q.*((ep.*r).^2+c^2).^(q-1);
rbfdy=@(c,r,dy)2*dy.*(ep^2).*q.*((ep.*r).^2+c^2).^(q-1);
rbf2dx=@(c,r,dx)2.*(ep^2).*q.*(((ep.*r).^2+c^2).^(q-1))+4.*(ep^4).*(dx.^2).*(q.*(q-1)).*(((ep.*r).^2+c^2).^(q-2));
rbf2dy=@(c,r,dy)2.*(ep^2).*q.*(((ep.*r).^2+c^2).^(q-1))+4.*(ep^4).*(dy.^2).*(q.*(q-1)).*(((ep.*r).^2+c^2).^(q-2));
rbfdxy=@(c,r,dx,dy)4*(ep^4).*(dx.*dy).*(q.*(q-1)).*(((ep.*r).^2+c^2).^(q-2));

TEORIA=1;

if TEORIA==1
    
[A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,k,E,h,A,B,D]=DRBF_FDST_BAT;

pf=2/sqrt(XNDF+1);      
AA=rbf(pf,DistM); Ax=rbfdx(pf,DistM,DiffMX); Ay=rbfdy(pf,DistM,DiffMY);
Axx=rbf2dx(pf,DistM,DiffMX); Ayy=rbf2dy(pf,DistM,DiffMY); Axy=rbfdxy(pf,DistM,DiffMX,DiffMY);

elseif TEORIA==2

[A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,k,E,h,A,B,D]=DRBF_FDST_BAT;
     
AA=rbf; Ax=rbfdx; Ay=rbfdy;
Axx=rbf2dx; Ayy=rbf2dy; Axy=rbfdxy;

save DistM
save DiffMX
save DiffMY

elseif TEORIA==3
    
    [A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,h,E]=DRBF_VIB_BAT;
            
   pf=2/sqrt(XNDF+1);      
   AA=rbf(pf,DistM); Ax=rbfdx(pf,DistM,DiffMX); Ay=rbfdy(pf,DistM,DiffMY);
   Axx=rbf2dx(pf,DistM,DiffMX); Ayy=rbf2dy(pf,DistM,DiffMY); Axy=rbfdxy(pf,DistM,DiffMX,DiffMY);

  A=[];B=[];D=[];k=5/6;

end

    
end


function DM = DistanceMatrix(dsites,ctrs)
[M,s] = size(dsites); [N,s] = size(ctrs);

DM = zeros(M,N);
for d=1:s
[dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
DM = DM + (dr-cc).^2;
end
DM = sqrt(DM);

end

function DM = DifferenceMatrix(dsites,ctrs)

[dr,cc] = ndgrid(dsites,ctrs);
DM = dr-cc;

end