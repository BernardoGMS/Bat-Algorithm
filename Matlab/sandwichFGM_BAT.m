function [A11,A12,A13,A14,A15,A21,A22,A23,A24,A25,A31,A32,A33,A34,A35,A41,A42,A43,A44,A45,A51,A52,A53,A54,A55...
                ,B11,B12,B13,B21,B22,B23,B31,B32,B33...
                ,D11,D12,D13,D21,D22,D23,D31,D32,D33,MASSA0,MASSA1,MASSA2,et,Young_core]=sandwichFGM_BAT

a=1;%p=9.9556665725;
razao_a_h=5;
    
et=a/razao_a_h;
% ec=et/2;

Young_core=70e9;
Young_inf=200e9;
Young_sup=Young_inf;
poisson_core=0.3;
poisson_inf=poisson_core;
poisson_sup=poisson_core;
rho_sup=5700;
rho_inf=rho_sup;
rho_core=2702;

G_core=Young_core/(2*(1+poisson_core));
G_inf=Young_inf/(2*(1+poisson_inf));
G_sup=G_inf;

syms hx x1
x2=sym('x2','positive');

VolMix1=vpa(((hx-(-et/2))/(-x1/2-(-et/2)))^x2); VolMix3=vpa(((hx-et/2)/(x1/2-et/2))^x2);

EEMix1=(Young_core-Young_inf)*VolMix1+Young_inf; EEMix2=(Young_core-Young_sup)*VolMix3+Young_sup;
NUMix1=(poisson_core-poisson_inf)*VolMix1+poisson_inf; NUMix2=(poisson_core-poisson_sup)*VolMix3+poisson_sup;
GMix1=(G_core-G_inf)*VolMix1+G_inf; GMix2=(G_core-G_sup)*VolMix3+G_sup;
MMix1=(rho_core-rho_inf)*VolMix1+rho_inf; MMix2=(rho_core-rho_sup)*VolMix3+rho_sup;

FG=int(GMix1,hx,-et/2,-x1/2)+int(G_core,hx,-x1/2,x1/2)+int(GMix2,hx,x1/2,et/2);
BFG=int(hx*GMix1,hx,-et/2,-x1/2)+int(hx*G_core,hx,-x1/2,x1/2)+int(hx*GMix2,hx,x1/2,et/2);
DFG=int(hx^2*GMix1,hx,-et/2,-x1/2)+int(hx^2*G_core,hx,-x1/2,x1/2)+int(hx^2*GMix2,hx,x1/2,et/2);

NUEf=int(NUMix1,hx,-et/2,-x1/2)+int(poisson_core,hx,-x1/2,x1/2)+int(NUMix2,hx,x1/2,et/2);
v12=NUEf;v21=v12;

AEf=int(EEMix1,hx,-et/2,-x1/2)+int(Young_core,hx,-x1/2,x1/2)+int(EEMix2,hx,x1/2,et/2);
% AEf=subs(AEf,[x1,x2],[ec,p]);
AEf=AEf/(1-v12*v21);
% GEf=subs(FG,[x1,x2],[ec,p]);

BEf=int(hx*EEMix1,hx,-et/2,-x1/2)+int(hx*Young_core,hx,-x1/2,x1/2)+int(hx*EEMix2,hx,x1/2,et/2);
% BEf=subs(BEf,[x1,x2],[ec,p]);
BEf=BEf/(1-v12*v21);
% BGEf=subs(BFG,[x1,x2],[ec,p]);

DEf=int(hx^2*EEMix1,hx,-et/2,-x1/2)+int(hx^2*Young_core,hx,-x1/2,x1/2)+int(hx^2*EEMix2,hx,x1/2,et/2);
% DEf=subs(DEf,[x1,x2],[ec,p]);
DEf=DEf/(1-v12*v21);
% DGEf=subs(DFG,[x1,x2],[ec,p]);

MEf0=int(MMix1,hx,-et/2,-x1/2)+int(rho_core,hx,-x1/2,x1/2)+int(MMix2,hx,x1/2,et/2);%MEf0=subs(MEf0,[x1,x2],[ec,p]);
MEf1=int(hx*MMix1,hx,-et/2,-x1/2)+int(hx*rho_core,hx,-x1/2,x1/2)+int(hx*MMix2,hx,x1/2,et/2);%MEf1=subs(MEf1,[x1,x2],[ec,p]);
MEf2=int(hx^2*MMix1,hx,-et/2,-x1/2)+int(hx^2*rho_core,hx,-x1/2,x1/2)+int(hx^2*MMix2,hx,x1/2,et/2);%MEf2=subs(MEf2,[x1,x2],[ec,p]);

A=sym(zeros(5,5));B=sym(zeros(3,3));D=sym(zeros(5,5));

A(1,1)=AEf;A(1,2)=v12*A(1,1);
A(2,1)=A(1,2);A(2,2)=A(1,1);
A(3,3)=FG;A(4,4)=FG;A(5,5)=FG;
% A(3,3)=GEf;A(4,4)=GEf;A(5,5)=GEf;

B(1,1)=BEf;B(1,2)=v12*B(1,1);
B(2,1)=B(1,2);B(2,2)=B(1,1);
B(3,3)=BFG;
% B(3,3)=BGEf;

D(1,1)=DEf;D(1,2)=v12*D(1,1);
D(2,1)=D(1,2);D(2,2)=D(1,1);
D(3,3)=DFG;D(4,4)=DFG;D(5,5)=DFG;
% D(3,3)=DGEf;D(4,4)=DGEf;D(5,5)=DGEf;

% A=double(A);B=double(B);D=double(D);

A11=matlabFunction(A(1,1));A12=matlabFunction(A(1,2));A13=0;A14=0;A15=0;
A21=matlabFunction(A(2,1));A22=matlabFunction(A(2,2));A23=0;A24=0;A25=0;
A31=0;A32=0;A33=matlabFunction(A(3,3));A34=0;A35=0;
A41=0;A42=0;A43=0;A44=matlabFunction(A(4,4));A45=0;
A51=0;A52=0;A53=0;A54=0;A55=matlabFunction(A(5,5));
    
B11=matlabFunction(B(1,1));B12=matlabFunction(B(1,2));B13=0;
B21=matlabFunction(B(2,1));B22=matlabFunction(B(2,2));B23=0;
B31=0;B32=0;B33=matlabFunction(B(3,3));B34=0;B35=0;

D11=matlabFunction(D(1,1));D12=matlabFunction(D(1,2));D13=0;
D21=matlabFunction(D(2,1));D22=matlabFunction(D(2,2));D23=0;
D31=0;D32=0;D33=matlabFunction(D(3,3));D34=0;D35=0;

MASSA0=matlabFunction(MEf0);MASSA1=matlabFunction(MEf1);MASSA2=matlabFunction(MEf2);


end