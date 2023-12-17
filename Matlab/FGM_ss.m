function SW=FGM_ss(A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,k,vp,AA,Ax,Ay,Axx,Ayy,Axy,NT,NF,XNDF,N)

%1 ---> SANDWICH; 2 ---> CAMADA FGM; 3 ---> FUNCAO EXPLICITA 
tipo_prop=1;

if tipo_prop==1
    
A=zeros(5,5);B=zeros(3,3);D=zeros(3,3);       
A(1,1)=A11(vp(1),vp(2));A(1,2)=A12(vp(1),vp(2));A(1,3)=0;A(1,4)=0;A(1,5)=0;
A(2,1)=A21(vp(1),vp(2));A(2,2)=A22(vp(1),vp(2));A(2,3)=0;A(2,4)=0;A(2,5)=0;
A(3,1)=0;A(3,2)=0;A(3,3)=A33(vp(1),vp(2));A(3,4)=0;A(3,5)=0;
A(4,1)=0;A(4,2)=0;A(4,3)=0;A(4,4)=A44(vp(1),vp(2));A(4,5)=0;
A(5,1)=0;A(5,2)=0;A(5,3)=0;A(5,4)=0;A(5,5)=A55(vp(1),vp(2));

% B(1,1)=B11(vp(1),vp(2));B(1,2)=B12(vp(1),vp(2));B(1,3)=0;
% B(2,1)=B21(vp(1),vp(2));B(2,2)=B22(vp(1),vp(2));B(2,3)=0;
% B(3,1)=0;B(3,2)=0;B(3,3)=B33(vp(1),vp(2));

B(1,1)=0;B(1,2)=0;B(1,3)=0;
B(2,1)=0;B(2,2)=0;B(2,3)=0;
B(3,1)=0;B(3,2)=0;B(3,3)=0;

D(1,1)=D11(vp(1),vp(2));D(1,2)=D12(vp(1),vp(2));D(1,3)=0;
D(2,1)=D21(vp(1),vp(2));D(2,2)=D22(vp(1),vp(2));D(2,3)=0;
D(3,1)=0;D(3,2)=0;D(3,3)=D33(vp(1),vp(2));
    
else
            
A=zeros(5,5);B=zeros(3,3);D=zeros(3,3);       
A(1,1)=A11(vp);A(1,2)=A12(vp);A(1,3)=0;A(1,4)=0;A(1,5)=0;
A(2,1)=A21(vp);A(2,2)=A22(vp);A(2,3)=0;A(2,4)=0;A(2,5)=0;
A(3,1)=0;A(3,2)=0;A(3,3)=A33(vp);A(3,4)=0;A(3,5)=0;
A(4,1)=0;A(4,2)=0;A(4,3)=0;A(4,4)=A44(vp);A(4,5)=0;
A(5,1)=0;A(5,2)=0;A(5,3)=0;A(5,4)=0;A(5,5)=A55(vp);

B(1,1)=B11(vp);B(1,2)=B12(vp);B(1,3)=0;
B(2,1)=B21(vp);B(2,2)=B22(vp);B(2,3)=0;
B(3,1)=0;B(3,2)=0;B(3,3)=B33(vp);

D(1,1)=D11(vp);D(1,2)=D12(vp);D(1,3)=0;
D(2,1)=D21(vp);D(2,2)=D22(vp);D(2,3)=0;
D(3,1)=0;D(3,2)=0;D(3,3)=D33(vp);
end

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

p0=1;

ff=zeros(5*NT,1); 
ff(2*NT+NF+1:3*NT,1)=p0;

lambda=K\ff;

SWT=AA*lambda(2*NT+1:3*NT,1);

SW=max(abs(SWT));

% indice=find(-SW==SWT);
% 
% fprintf('INDICE\n')
% disp(indice)


end