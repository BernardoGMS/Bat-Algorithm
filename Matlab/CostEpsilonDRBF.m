function ceps=CostEpsilonDRBF(vp,k,XNDF,NT,NF,rbf,rbfdx,rbfdy,rbf2dx,rbf2dy,rbfdxy,A,B,D)

load DistM
load DiffMX
load DiffMY

AA=rbf(vp,DistM); Ax=rbfdx(vp,DistM,DiffMX); Ay=rbfdy(vp,DistM,DiffMY);
Axx=rbf2dx(vp,DistM,DiffMX); Ayy=rbf2dy(vp,DistM,DiffMY); Axy=rbfdxy(vp,DistM,DiffMX,DiffMY);

% propriedades=0;

% Construcao de K:

% if propriedades==0
% 
% [A,B,D,~,~,~,~,h,E,k]=layerwisetheory1;%Laminados
% 
% k=5/6;
% 
% elseif propriedades==1
% 
% %FGM tem MortiTanaka para sandwich (rever)
% [A,B,D,~,~,~,h,E]=sandwichFGM;%Sandwich FGM - Sem MoriTanaka
% k=5/6;
% 
% A=double(A);B=double(B);D=double(D);
% 
% elseif propriedades==2
%     
% [A,B,D,E,h]=FGMproperties;%Camada FGM (Misturas + MoriTanaka) - Sem Massas
% k=5/6;
% 
% end

% 1� equacao (grau de liberdade u0):
K11=A(1,1)*Axx+2*A(1,3)*Axy+A(3,3)*Ayy;K12=A(1,2)*Axy+A(1,3)*Axx+A(2,3)*Ayy+A(3,3)*Axy; K13=zeros(NT,NT); K14=B(1,1)*Axx+2*B(1,3)*Axy+B(3,3)*Ayy; K15=B(1,2)*Axy+B(1,3)*Axx+B(2,3)*Ayy+B(3,3)*Axy;
% 2� equacao (grau de liberdade v0):
K21=A(1,3)*Axx+A(3,3)*Axy+A(1,2)*Axy+A(2,3)*Ayy; K22=A(3,3)*Axx+A(2,2)*Ayy+2*A(2,3)*Axy; K23=zeros(NT,NT); K24=B(1,3)*Axx+B(3,3)*Axy+B(1,2)*Axy+B(2,3)*Ayy; K25=2*B(2,3)*Axy+B(3,3)*Axx+B(2,2)*Ayy;
% 3� equacao (grau de liberdade w0):
K31=zeros(NT,NT); K32=zeros(NT,NT); K33=k*(2*A(4,5)*Axy+A(5,5)*Axx+A(4,4)*Ayy); K34=k*(A(5,5)*Ax+A(4,5)*Ay); K35=k*(A(4,5)*Ax+A(4,4)*Ay);
% 4� equacao (grau de liberdade phi_x):
K41=B(1,1)*Axx+2*B(1,3)*Axy+B(3,3)*Ayy; K42=B(1,3)*Axx+(B(1,2)+B(3,3))*Axy+B(2,3)*Ayy; K43=-k*(A(4,5)*Ay+A(5,5)*Ax); K44=D(1,1)*Axx+2*D(1,3)*Axy+D(3,3)*Ayy-k*A(5,5)*AA; K45=D(1,3)*Axx+(D(1,2)+D(3,3))*Axy+D(2,3)*Ayy-k*A(4,5)*AA;
% 5� equacao (grau de liberdade phi_y):
K51=B(1,3)*Axx+(B(1,2)+B(1,3))*Axy+B(2,3)*Ayy; K52=B(3,3)*Axx+2*B(2,3)*Axy+B(2,2)*Ayy; K53=-k*(A(4,4)*Ay+A(4,5)*Ax); K54=D(1,3)*Axx+(D(1,2)+D(3,3))*Axy+D(2,3)*Ayy-k*A(4,5)*AA; K55=D(3,3)*Axx+2*D(2,3)*Axy+D(2,2)*Ayy-k*A(4,4)*AA;

K=[K11 K12 K13 K14 K15;...
    K21 K22 K23 K24 K25;...
    K31 K32 K33 K34 K35;...
    K41 K42 K43 K44 K45;...
    K51 K52 K53 K54 K55];

boundary=1;

if boundary==1 %Simplesmente apoiado
    
%u0=0:
K(2*(XNDF+1)+1:NF,1:NT)=AA(2*(XNDF+1)+1:NF,:);
K(2*(XNDF+1)+1:NF,NT+1:2*NT)=0; 
K(2*(XNDF+1)+1:NF,2*NT+1:3*NT)=0;
K(2*(XNDF+1)+1:NF,3*NT+1:4*NT)=0;
K(2*(XNDF+1)+1:NF,4*NT+1:5*NT)=0;

%u0=0 nos pontos dos cantos
K(1,1:NT)=AA(1,:);K(1,NT+1:2*NT)=0;K(1,2*NT+1:3*NT)=0;K(1,3*NT+1:4*NT)=0;K(1,4*NT+1:5*NT)=0;
K(XNDF+1,1:NT)=AA(XNDF+1,:);K(XNDF+1,NT+1:2*NT)=0;K(XNDF+1,2*NT+1:3*NT)=0;K(XNDF+1,3*NT+1:4*NT)=0;K(XNDF+1,4*NT+1:5*NT)=0;
K(XNDF+2,1:NT)=AA(XNDF+2,:);K(XNDF+2,NT+1:2*NT)=0;K(XNDF+2,2*NT+1:3*NT)=0;K(XNDF+2,3*NT+1:4*NT)=0;K(XNDF+2,4*NT+1:5*NT)=0;
K(2*(XNDF+1),1:NT)=AA(2*(XNDF+1),:);K(2*(XNDF+1),NT+1:2*NT)=0;K(2*(XNDF+1),2*NT+1:3*NT)=0;K(2*(XNDF+1),3*NT+1:4*NT)=0;K(2*(XNDF+1),4*NT+1:5*NT)=0;

%% 2-----------------------1
%Nxx=0
% K(2:2*(XNDF+1)-1,1:NT)=A(1,1)*Ax(2:2*(XNDF+1)-1,:)+A(1,3)*Ay(2:2*(XNDF+1)-1,:);
% K(2:2*(XNDF+1)-1,NT+1:2*NT)=A(1,2)*Ay(2:2*(XNDF+1)-1,:)+A(1,3)*Ax(2:2*(XNDF+1)-1,:); 
% K(2:2*(XNDF+1)-1,2*NT+1:3*NT)=0;
% K(2:2*(XNDF+1)-1,3*NT+1:4*NT)=B(1,1)*Ax(2:2*(XNDF+1)-1,:)+B(1,3)*Ay(2:2*(XNDF+1)-1,:);
% K(2:2*(XNDF+1)-1,4*NT+1:5*NT)=B(1,2)*Ay(2:2*(XNDF+1)-1,:)+B(1,3)*Ax(2:2*(XNDF+1)-1,:);
%Nxx=0 nos cantos inclusive
K(1:2*(XNDF+1),1:NT)=A(1,1)*Ax(1:2*(XNDF+1),:)+A(1,3)*Ay(1:2*(XNDF+1),:);
K(1:2*(XNDF+1),NT+1:2*NT)=A(1,2)*Ay(1:2*(XNDF+1),:)+A(1,3)*Ax(1:2*(XNDF+1),:); 
K(1:2*(XNDF+1),2*NT+1:3*NT)=0;
K(1:2*(XNDF+1),3*NT+1:4*NT)=B(1,1)*Ax(1:2*(XNDF+1),:)+B(1,3)*Ay(1:2*(XNDF+1),:);
K(1:2*(XNDF+1),4*NT+1:5*NT)=B(1,2)*Ay(1:2*(XNDF+1),:)+B(1,3)*Ax(1:2*(XNDF+1),:);

%v0=0:
% K(NT+2:NT+2*(XNDF+1)-1,1:NT)=0;
% K(NT+2:NT+2*(XNDF+1)-1,NT+1:2*NT)=AA(2:2*(XNDF+1)-1,:);
% K(NT+2:NT+2*(XNDF+1)-1,2*NT+1:3*NT)=0;
% K(NT+2:NT+2*(XNDF+1)-1,3*NT+1:4*NT)=0; 
% K(NT+2:NT+2*(XNDF+1)-1,4*NT+1:5*NT)=0;
%v0=0 nos pontos dos cantos inclusive
K(NT+1:NT+2*(XNDF+1),1:NT)=0;
K(NT+1:NT+2*(XNDF+1),NT+1:2*NT)=AA(1:2*(XNDF+1),:);
K(NT+1:NT+2*(XNDF+1),2*NT+1:3*NT)=0;
K(NT+1:NT+2*(XNDF+1),3*NT+1:4*NT)=0; 
K(NT+1:NT+2*(XNDF+1),4*NT+1:5*NT)=0;

%Nyy=0
K(NT+2*(XNDF+1)+1:NT+NF,1:NT)=A(1,2)*Ax(2*(XNDF+1)+1:NF,:)+A(2,3)*Ay(2*(XNDF+1)+1:NF,:);
K(NT+2*(XNDF+1)+1:NT+NF,NT+1:2*NT)=A(2,2)*Ay(2*(XNDF+1)+1:NF,:)+A(2,3)*Ax(2*(XNDF+1)+1:NF,:);
K(NT+2*(XNDF+1)+1:NT+NF,2*NT+1:3*NT)=0;
K(NT+2*(XNDF+1)+1:NT+NF,3*NT+1:4*NT)=B(1,2)*Ax(2*(XNDF+1)+1:NF,:)+B(2,3)*Ay(2*(XNDF+1)+1:NF,:); 
K(NT+2*(XNDF+1)+1:NT+NF,4*NT+1:5*NT)=B(2,2)*Ay(2*(XNDF+1)+1:NF,:)+B(2,3)*Ax(2*(XNDF+1)+1:NF,:);

%% !!!!!!!!!!!!!!!!!!!!!!!!!!
% %Nyy=0 nos pontos dos cantos
K(NT+1,1:NT)=A(1,2)*Ax(1,:)+A(2,3)*Ay(1,:);
K(NT+1,NT+1:2*NT)=A(2,2)*Ay(1,:)+A(2,3)*Ax(1,:);
K(NT+1,2*NT+1:3*NT)=0;
K(NT+1,3*NT+1:4*NT)=B(1,2)*Ax(1,:)+B(2,3)*Ay(1,:); 
K(NT+1,4*NT+1:5*NT)=B(2,2)*Ay(1,:)+B(2,3)*Ax(1,:);

K(NT+XNDF+1,1:NT)=A(1,2)*Ax(XNDF+1,:)+A(2,3)*Ay(XNDF+1,:);
K(NT+XNDF+1,NT+1:2*NT)=A(2,2)*Ay(XNDF+1,:)+A(2,3)*Ax(XNDF+1,:);
K(NT+XNDF+1,2*NT+1:3*NT)=0;
K(NT+XNDF+1,3*NT+1:4*NT)=B(1,2)*Ax(XNDF+1,:)+B(2,3)*Ay(XNDF+1,:);
K(NT+XNDF+1,4*NT+1:5*NT)=B(2,2)*Ay(XNDF+1,:)+B(2,3)*Ax(XNDF+1,:);

K(NT+XNDF+2,1:NT)=A(1,2)*Ax(XNDF+2,:)+A(2,3)*Ay(XNDF+2,:);
K(NT+XNDF+2,NT+1:2*NT)=A(2,2)*Ay(XNDF+2,:)+A(2,3)*Ax(XNDF+2,:);
K(NT+XNDF+2,2*NT+1:3*NT)=0;
K(NT+XNDF+2,3*NT+1:4*NT)=B(1,2)*Ax(XNDF+2,:)+B(2,3)*Ay(XNDF+2,:);
K(NT+XNDF+2,4*NT+1:5*NT)=B(2,2)*Ay(XNDF+2,:)+B(2,3)*Ax(XNDF+2,:);

K(NT+2*(XNDF+1),1:NT)=A(1,2)*Ax(2*(XNDF+1),:)+A(2,3)*Ay(2*(XNDF+1),:);
K(NT+2*(XNDF+1),NT+1:2*NT)=A(2,2)*Ay(2*(XNDF+1),:)+A(2,3)*Ax(2*(XNDF+1),:);
K(NT+2*(XNDF+1),2*NT+1:3*NT)=0;
K(NT+2*(XNDF+1),3*NT+1:4*NT)=B(1,2)*Ax(2*(XNDF+1),:)+B(2,3)*Ay(2*(XNDF+1),:);
K(NT+2*(XNDF+1),4*NT+1:5*NT)=B(2,2)*Ay(2*(XNDF+1),:)+B(2,3)*Ax(2*(XNDF+1),:);

%w0=0 nos cantos inclusive
K(2*NT+1:2*NT+NF,1:NT)=0;
K(2*NT+1:2*NT+NF,NT+1:2*NT)=0;
K(2*NT+1:2*NT+NF,2*NT+1:3*NT)=AA(1:NF,:);
K(2*NT+1:2*NT+NF,3*NT+1:4*NT)=0; 
K(2*NT+1:2*NT+NF,4*NT+1:5*NT)=0;

%w0=0:
% K(2*NT+2:2*NT+2*(XNDF+1)-1,1:NT)=0;
% K(2*NT+2:2*NT+2*(XNDF+1)-1,NT+1:2*NT)=0;
% K(2*NT+2:2*NT+2*(XNDF+1)-1,2*NT+1:3*NT)=AA(2:2*(XNDF+1)-1,:);
% K(2*NT+2:2*NT+2*(XNDF+1)-1,3*NT+1:4*NT)=0; 
% K(2*NT+2:2*NT+2*(XNDF+1)-1,4*NT+1:5*NT)=0;
% 
% K(2*NT+2*(XNDF+1)+1:2*NT+NF,1:NT)=0;
% K(2*NT+2*(XNDF+1)+1:2*NT+NF,NT+1:2*NT)=0;
% K(2*NT+2*(XNDF+1)+1:2*NT+NF,2*NT+1:3*NT)=AA(2*(XNDF+1)+1:NF,:);
% K(2*NT+2*(XNDF+1)+1:2*NT+NF,3*NT+1:4*NT)=0; 
% K(2*NT+2*(XNDF+1)+1:2*NT+NF,4*NT+1:5*NT)=0;

%phi_x=0:
K(3*NT+2*(XNDF+1)+1:3*NT+NF,1:NT)=0;
K(3*NT+2*(XNDF+1)+1:3*NT+NF,NT+1:2*NT)=0;
K(3*NT+2*(XNDF+1)+1:3*NT+NF,2*NT+1:3*NT)=0;
K(3*NT+2*(XNDF+1)+1:3*NT+NF,3*NT+1:4*NT)=AA(2*(XNDF+1)+1:NF,:);
K(3*NT+2*(XNDF+1)+1:3*NT+NF,4*NT+1:5*NT)=0;

%phi_x=0 nos pontos dos cantos
K(3*NT+1,1:NT)=0;K(3*NT+1,NT+1:2*NT)=0;K(3*NT+1,2*NT+1:3*NT)=0;K(3*NT+1,3*NT+1:4*NT)=AA(1,:);K(3*NT+1,4*NT+1:5*NT)=0;
K(3*NT+XNDF+1,1:NT)=0;K(3*NT+XNDF+1,NT+1:2*NT)=0;K(3*NT+XNDF+1,2*NT+1:3*NT)=0;K(3*NT+XNDF+1,3*NT+1:4*NT)=AA(XNDF+1,:);K(3*NT+XNDF+1,4*NT+1:5*NT)=0;
K(3*NT+XNDF+2,1:NT)=0;K(3*NT+XNDF+2,NT+1:2*NT)=0;K(3*NT+XNDF+2,2*NT+1:3*NT)=0;K(3*NT+XNDF+2,3*NT+1:4*NT)=AA(XNDF+2,:);K(3*NT+XNDF+2,4*NT+1:5*NT)=0;
K(3*NT+2*(XNDF+1),1:NT)=0;K(3*NT+2*(XNDF+1),NT+1:2*NT)=0;K(3*NT+2*(XNDF+1),2*NT+1:3*NT)=0;K(3*NT+2*(XNDF+1),3*NT+1:4*NT)=AA(2*(XNDF+1),:);K(3*NT+2*(XNDF+1),4*NT+1:5*NT)=0;

%% !!!!!!!!!!!!!!!!!!!!!!!!!!
% %Mxx=0
% K(3*NT+1:3*NT+2*(XNDF+1),1:NT)=B(1,1)*Ax(1:2*(XNDF+1),:)+B(1,3)*Ay(1:2*(XNDF+1),:);
% K(3*NT+1:3*NT+2*(XNDF+1),NT+1:2*NT)=B(1,2)*Ay(1:2*(XNDF+1),:)+B(1,3)*Ax(1:2*(XNDF+1),:);
% K(3*NT+1:3*NT+2*(XNDF+1),2*NT+1:3*NT)=0;
% K(3*NT+1:3*NT+2*(XNDF+1),3*NT+1:4*NT)=D(1,1)*Ax(1:2*(XNDF+1),:)+D(1,3)*Ay(1:2*(XNDF+1),:);
% K(3*NT+1:3*NT+2*(XNDF+1),4*NT+1:5*NT)=D(1,2)*Ay(1:2*(XNDF+1),:)+D(1,3)*Ax(1:2*(XNDF+1),:);

% 2-----------------------1
%Mxx=0
% K(3*NT+2:3*NT+2*(XNDF+1)-1,1:NT)=B(1,1)*Ax(2:2*(XNDF+1)-1,:)+B(1,3)*Ay(2:2*(XNDF+1)-1,:);
% K(3*NT+2:3*NT+2*(XNDF+1)-1,NT+1:2*NT)=B(1,2)*Ay(2:2*(XNDF+1)-1,:)+B(1,3)*Ax(2:2*(XNDF+1)-1,:);
% K(3*NT+2:3*NT+2*(XNDF+1)-1,2*NT+1:3*NT)=0;
% K(3*NT+2:3*NT+2*(XNDF+1)-1,3*NT+1:4*NT)=D(1,1)*Ax(2:2*(XNDF+1)-1,:)+D(1,3)*Ay(2:2*(XNDF+1)-1,:);
% K(3*NT+2:3*NT+2*(XNDF+1)-1,4*NT+1:5*NT)=D(1,2)*Ay(2:2*(XNDF+1)-1,:)+D(1,3)*Ax(2:2*(XNDF+1)-1,:);
% %Mxx=0 nos cantos inclusive
K(3*NT+1:3*NT+2*(XNDF+1),1:NT)=B(1,1)*Ax(1:2*(XNDF+1),:)+B(1,3)*Ay(1:2*(XNDF+1),:);
K(3*NT+1:3*NT+2*(XNDF+1),NT+1:2*NT)=B(1,2)*Ay(1:2*(XNDF+1),:)+B(1,3)*Ax(1:2*(XNDF+1),:);
K(3*NT+1:3*NT+2*(XNDF+1),2*NT+1:3*NT)=0;
K(3*NT+1:3*NT+2*(XNDF+1),3*NT+1:4*NT)=D(1,1)*Ax(1:2*(XNDF+1),:)+D(1,3)*Ay(1:2*(XNDF+1),:);
K(3*NT+1:3*NT+2*(XNDF+1),4*NT+1:5*NT)=D(1,2)*Ay(1:2*(XNDF+1),:)+D(1,3)*Ax(1:2*(XNDF+1),:);

%phi_y=0 nos cantos inclusive
K(4*NT+1:4*NT+2*(XNDF+1),1:NT)=0;
K(4*NT+1:4*NT+2*(XNDF+1),NT+1:2*NT)=0;
K(4*NT+1:4*NT+2*(XNDF+1),2*NT+1:3*NT)=0;
K(4*NT+1:4*NT+2*(XNDF+1),3*NT+1:4*NT)=0;
K(4*NT+1:4*NT+2*(XNDF+1),4*NT+1:5*NT)=AA(1:2*(XNDF+1),:);
%phi_y=0:
% K(4*NT+2:4*NT+2*(XNDF+1)-1,1:NT)=0;
% K(4*NT+2:4*NT+2*(XNDF+1)-1,NT+1:2*NT)=0;
% K(4*NT+2:4*NT+2*(XNDF+1)-1,2*NT+1:3*NT)=0;
% K(4*NT+2:4*NT+2*(XNDF+1)-1,3*NT+1:4*NT)=0;
% K(4*NT+2:4*NT+2*(XNDF+1)-1,4*NT+1:5*NT)=AA(2:2*(XNDF+1)-1,:);

%Myy=0
K(4*NT+2*(XNDF+1)+1:4*NT+NF,1:NT)=B(1,2)*Ax(2*(XNDF+1)+1:NF,:)+B(2,3)*Ay(2*(XNDF+1)+1:NF,:);
K(4*NT+2*(XNDF+1)+1:4*NT+NF,NT+1:2*NT)=B(2,2)*Ay(2*(XNDF+1)+1:NF,:)+B(2,3)*Ax(2*(XNDF+1)+1:NF,:);
K(4*NT+2*(XNDF+1)+1:4*NT+NF,2*NT+1:3*NT)=0;
K(4*NT+2*(XNDF+1)+1:4*NT+NF,3*NT+1:4*NT)=D(1,2)*Ax(2*(XNDF+1)+1:NF,:)+D(2,3)*Ay(2*(XNDF+1)+1:NF,:);
K(4*NT+2*(XNDF+1)+1:4*NT+NF,4*NT+1:5*NT)=D(2,2)*Ay(2*(XNDF+1)+1:NF,:)+D(2,3)*Ax(2*(XNDF+1)+1:NF,:);

%% !!!!!!!!!!!!!!!!!!!!!!!!!!
%Myy=0 nos pontos dos cantos
K(4*NT+1,1:NT)=B(1,2)*Ax(1,:)+B(2,3)*Ay(1,:);
K(4*NT+1,NT+1:2*NT)=B(2,2)*Ay(1,:)+B(2,3)*Ax(1,:);
K(4*NT+1,2*NT+1:3*NT)=0;
K(4*NT+1,3*NT+1:4*NT)=D(1,2)*Ax(1,:)+D(2,3)*Ay(1,:);
K(4*NT+1,4*NT+1:5*NT)=D(2,2)*Ay(1,:)+D(2,3)*Ax(1,:);

K(4*NT+XNDF+1,1:NT)=B(1,2)*Ax(XNDF+1,:)+B(2,3)*Ay(XNDF+1,:);
K(4*NT+XNDF+1,NT+1:2*NT)=B(2,2)*Ay(XNDF+1,:)+B(2,3)*Ax(XNDF+1,:);
K(4*NT+XNDF+1,2*NT+1:3*NT)=0;
K(4*NT+XNDF+1,3*NT+1:4*NT)=D(1,2)*Ax(XNDF+1,:)+D(2,3)*Ay(XNDF+1,:);
K(4*NT+XNDF+1,4*NT+1:5*NT)=D(2,2)*Ay(XNDF+1,:)+D(2,3)*Ax(XNDF+1,:);

K(4*NT+XNDF+2,1:NT)=B(1,2)*Ax(XNDF+2,:)+B(2,3)*Ay(XNDF+2,:);
K(4*NT+XNDF+2,NT+1:2*NT)=B(2,2)*Ay(XNDF+2,:)+B(2,3)*Ax(XNDF+2,:);
K(4*NT+XNDF+2,2*NT+1:3*NT)=0;
K(4*NT+XNDF+2,3*NT+1:4*NT)=D(1,2)*Ax(XNDF+2,:)+D(2,3)*Ay(XNDF+2,:);
K(4*NT+XNDF+2,4*NT+1:5*NT)=D(2,2)*Ay(XNDF+2,:)+D(2,3)*Ax(XNDF+2,:);

K(4*NT+2*(XNDF+1),1:NT)=B(1,2)*Ax(2*(XNDF+1),:)+B(2,3)*Ay(2*(XNDF+1),:);
K(4*NT+2*(XNDF+1),NT+1:2*NT)=B(2,2)*Ay(2*(XNDF+1),:)+B(2,3)*Ax(2*(XNDF+1),:);
K(4*NT+2*(XNDF+1),2*NT+1:3*NT)=0;
K(4*NT+2*(XNDF+1),3*NT+1:4*NT)=D(1,2)*Ax(2*(XNDF+1),:)+D(2,3)*Ay(2*(XNDF+1),:);
K(4*NT+2*(XNDF+1),4*NT+1:5*NT)=D(2,2)*Ay(2*(XNDF+1),:)+D(2,3)*Ax(2*(XNDF+1),:);
    
elseif boundary==2 %Encastramento
  
%uz=0:
K(1:NF,1:NT)=A(1:NF,:);
K(1:NF,NT+1:2*NT)=0; 
K(1:NF,2*NT+1:3*NT)=0;

%phi_x=0:
K(NT+1:NT+NF,1:NT)=0;
K(NT+1:NT+NF,NT+1:2*NT)=A(1:NF,:);
K(NT+1:NT+NF,2*NT+1:3*NT)=0; 

%phi_y=0:
K(2*NT+1:2*NT+NF,1:NT)=0;
K(2*NT+1:2*NT+NF,NT+1:2*NT)=0;
K(2*NT+1:2*NT+NF,2*NT+1:3*NT)=A(1:NF,:);
    
end

%% Construcao do vector de cargas

force=2;

p0=1;

if force==1;

    ff=zeros(5*NT,1);ff(NF+((NT-NF-1)/2+1),1)=NT*p0;
        
elseif force==2
    
    ff=zeros(5*NT,1); 
    ff(2*NT+NF+1:3*NT,1)=p0;
   
end

lambda=K\ff;
invA=pinv(K);EF=zeros(NT,1);
for i=1:NT
EF(i,1)=lambda(i,1)/invA(i,i);
end

% EF=(invA*ff)./repmat(diag(invA),1,3*NT);
ceps=norm(EF(:));
disp(vp)
disp(EF')

end
