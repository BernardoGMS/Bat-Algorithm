%fgm PROPERTIES
function [A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
    ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
    D11,D12,D13,D21,D22,D23,D31,D32,D33,Ez,MEf0,MEf1,MEf2,h] = FGMproperties_BAT% op�ao de regra de misturas ou torikama)
colordef white

Et=70e9;
Eb=200e9;
nut=.3;
nub=.3;
a=1;
rah=5;

%Gemometrical Properties
h=a/rah;
opcao=1;
%Mechanical Properties
EE_m=Et;
EE_c=Eb;

rhom=2702;
rhoc=5700;

%Fiber volume
syms z real
pp=sym('pp','positive');
Vc=((z+h/2)/h)^pp;

%ezplot(Vc,[-h/2,h/2]);

%Properties in function of z
if opcao==1 %misturas
Ez=Vc*EE_c+(1-Vc)*EE_m;
%ezplot(Ez,[-h/2,h/2]);
nuz=Vc*nub+(1-Vc)*nut;

A=sym(zeros(5,5));

Q=(Ez/(1-nuz^2))*[1,nuz,0;nuz,1,0;0,0,(1-nuz)/2];
G=(Ez/(2*(1+nuz)))*[1,0;0,1];

    A(1:3,1:3)=Q;
    A(4:5,4:5)=G;

%integra��o matriz A B D As
 A=int(A,z,-h/2,h/2);
 B=int(Q*z,z,-h/2,h/2);    
 D=int(Q*z^2,z,-h/2,h/2);
%  As=int(G,z,-h/2,h/2);
 
%  A=matlabFunction(A);
%  B=matlabFunction(B);
%  D=matlabFunction(D);
%  As=matlabFunction(As);

A11=A(1,1);A12=A(1,2);A13=A(1,3)*pp;A14=A(1,4)*pp;A15=A(1,5)*pp;
A21=A(2,1);A22=A(2,2);A23=A(2,3)*pp;A24=A(2,4)*pp;A25=A(2,5)*pp;
A31=A(3,1)*pp;A32=A(3,2)*pp;A33=A(3,3);A34=A(3,4)*pp;A35=A(3,5)*pp;
A41=A(4,1)*pp;A42=A(4,2)*pp;A43=A(4,3)*pp;A44=A(4,4);A45=A(4,5);
A51=A(5,1)*pp;A52=A(5,2)*pp;A53=A(5,3)*pp;A54=A(5,4);A55=A(5,5);
 
B11=B(1,1);B12=B(1,2);B13=B(1,3);
B21=B(2,1);B22=B(2,2);B23=B(2,3);
B31=B(3,1);B32=B(3,2);B33=B(3,3);

D11=D(1,1);D12=D(1,2);D13=D(1,3);
D21=D(2,1);D22=D(2,2);D23=D(2,3);
D31=D(3,1);D32=D(3,2);D33=D(3,3);

A11=matlabFunction(A11);A12=matlabFunction(A12);A13=matlabFunction(A13);A14=matlabFunction(A14);A15=matlabFunction(A15);
A21=matlabFunction(A21);A22=matlabFunction(A22);A23=matlabFunction(A23);A24=matlabFunction(A24);A25=matlabFunction(A25);
A31=matlabFunction(A31);A32=matlabFunction(A32);A33=matlabFunction(A33);A34=matlabFunction(A34);A35=matlabFunction(A35);
A41=matlabFunction(A41);A42=matlabFunction(A42);A43=matlabFunction(A43);A44=matlabFunction(A44);A45=matlabFunction(A45);
A51=matlabFunction(A51);A52=matlabFunction(A52);A53=matlabFunction(A53);A54=matlabFunction(A54);A55=matlabFunction(A55);
 
B11=matlabFunction(B11);B12=matlabFunction(B12);B13=matlabFunction(B13)*pp;
B21=matlabFunction(B21);B22=matlabFunction(B22);B23=matlabFunction(B23)*pp;
B31=matlabFunction(B31)*pp;B32=matlabFunction(B32)*pp;B33=matlabFunction(B33);

D11=matlabFunction(D11);D12=matlabFunction(D12);D13=matlabFunction(D13)*pp;
D21=matlabFunction(D21);D22=matlabFunction(D22);D23=matlabFunction(D23)*pp;
D31=matlabFunction(D31)*pp;D32=matlabFunction(D32)*pp;D33=matlabFunction(D33);

 
else  %moritanaka
    Vol2=Vc;
    G1=EE_m/(2*(1+nut)); G2=EE_c/(2*(1+nub));
K1=EE_m/(3*(1-2*nut)); K2=EE_c/(3*(1-2*nub));

f1=(G1*(9*K1+8*G1))/(6*(K1+2*G1));

K11s=(4*G1*K1 + 3*K1*K2 - 4*G1*K1*Vol2 + 4*G1*K2*Vol2)/(4*G1 + 3*K2 + 3*K1*Vol2 - 3*K2*Vol2);
GG11s=(G1*f1 + G1*G2 - G1*Vol2*f1 + G2*Vol2*f1)/(G2 + f1 + G1*Vol2 - G2*Vol2);
   
Ez=vpa((9*K11s*GG11s)/(3*K11s+GG11s));
nuz=vpa((3*K11s-2*GG11s)/(2*(3*K11s+GG11s)));

Q=(Ez/(1-nuz^2))*[1,nuz,0;nuz,1,0;0,0,(1-nuz)/2];
G=(Ez/(2*(1+nuz)))*[1,0;0,1];
% 
    Qb=Q*z;
    Qd=Q*z^2;

    g11 = matlabFunction(Q(1,1));
    g12 = matlabFunction(Q(1,2));
    g13 = matlabFunction(Q(1,3));
    
    g22 = matlabFunction(Q(2,2));
    g23 = matlabFunction(Q(2,3));
    
    g33 = matlabFunction(Q(3,3));
    
    gb11 = matlabFunction(Qb(1,1));
    gb12 = matlabFunction(Qb(1,2));
    gb13 = matlabFunction(Qb(1,3));
    
    gb22 = matlabFunction(Qb(2,2));
    gb23 = matlabFunction(Qb(2,3));
    
    gb33 = matlabFunction(Qb(3,3));
    
    gd11 = matlabFunction(Qd(1,1));
    gd12 = matlabFunction(Qd(1,2));
    gd13 = matlabFunction(Qd(1,3));
    
    gd22 = matlabFunction(Qd(2,2));
    gd23 = matlabFunction(Qd(2,3));
    
    gd33 = matlabFunction(Qd(3,3));
  
    
    gas11 = matlabFunction(G(1,1));
    gas22 = matlabFunction(G(2,2));
   
    
  
    %integra��o matriz A B D As
    
 A(1,1)=quadv(g11,-h/2,h/2);
 A(1,2)=quadv(g12,-h/2,h/2);
 A(1,3)=0;%integral(g13,-h/2,h/2);
 

 A(2,1)= A(1,2);
 A(2,2)=quadv(g22,-h/2,h/2);
 A(2,3)=0;%integral(g23,-h/2,h/2);
 
 A(3,1)= A(1,3);
 A(3,2)=A(2,3);
 A(3,3)=quadv(g33,-h/2,h/2);
 
 B(1,1)=quadv(gb11,-h/2,h/2);
 B(1,2)=quadv(gb12,-h/2,h/2);
 B(1,3)=0;%integral(gb13,-h/2,h/2);
 

 B(2,1)= B(1,2);
 B(2,2)=quadv(gb22,-h/2,h/2);
 B(2,3)=0;%integral(gb23,-h/2,h/2);
 
 B(3,1)= B(1,3);
 B(3,2)=B(2,3);
 B(3,3)=quadv(gb33,-h/2,h/2);
 
 D(1,1)=quadv(gd11,-h/2,h/2);
 D(1,2)=quadv(gd12,-h/2,h/2);
 D(1,3)=0;%integral(gd13,-h/2,h/2);
 

 D(2,1)= D(1,2);
 D(2,2)=quadv(gd22,-h/2,h/2);
 D(2,3)=0;%integral(gd23,-h/2,h/2);
 
 D(3,1)= D(1,3);
 D(3,2)=D(2,3);
 D(3,3)=quadv(gd33,-h/2,h/2);
 
 As(1,1)=quadv(gas11,-h/2,h/2);
 As(1,2)=0;%integral(gas12,-h/2,h/2);

 As(2,1)= As(1,2);
 As(2,2)=quadv(gas22,-h/2,h/2);

 A=round(A);
 B=round(B);
 D=round(D);
 As=round(As);
 A(4:5,4:5)=As;
 
end

MMix=(rhoc-rhom)*Vc+rhom;
MMix=int(MMix,z,-h/2,h/2);
MEf1=int(z*MMix,z,-h/2,h/2);
MEf2=int(z^2*MMix,z,-h/2,h/2);

MEf0 = matlabFunction(MMix);
MEf1 = matlabFunction(MEf1);
MEf2 = matlabFunction(MEf2);

end

%Factor de correc��o ao corte


