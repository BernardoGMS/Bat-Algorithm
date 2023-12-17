function [A,B,D,MEf0,MEf1,MEf2,h,Et]=FGM_Mori
 
colordef white

Et=380e9;
Eb=70e9;
nut=.3;
nub=.3;
rho_t=3800;
rho_b=2707;
p=10;
a=1;
rah=20;

%Gemometrical Properties
h=a/rah;
%Mechanical Properties
EE_m=Eb;
EE_c=Et;

%Fiber volume
syms z real
Vc=((z+h/2)/h)^p;
% ezplot(Vc,[-h/2 h/2])

    Vol2=Vc;
    G1=EE_m/(2*(1+nut)); G2=EE_c/(2*(1+nub));
    K1=EE_m/(3*(1-2*nut)); K2=EE_c/(3*(1-2*nub));

f1=(G1*(9*K1+8*G1))/(6*(K1+2*G1));

K11s=(4*G1*K1 + 3*K1*K2 - 4*G1*K1*Vol2 + 4*G1*K2*Vol2)/(4*G1 + 3*K2 + 3*K1*Vol2 - 3*K2*Vol2);
GG11s=(G1*f1 + G1*G2 - G1*Vol2*f1 + G2*Vol2*f1)/(G2 + f1 + G1*Vol2 - G2*Vol2);
   
Ez=vpa((9*K11s*GG11s)/(3*K11s+GG11s));
nuz=vpa((3*K11s-2*GG11s)/(2*(3*K11s+GG11s)));

% ezplot(Ez,[-h/2 h/2])
% ezplot(nuz,[-h/2 h/2])

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
   
    
  
    %integração matriz A B D As
    
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

    syms hx x1
%     x2=sym('x2','positive');
    
    Vol1=((hx+h/2)/h)^p;
 
MMix1=(rho_t-rho_b)*Vol1+rho_b;% ezplot(MMix1,[-h/2 h/2])
MEf0=int(MMix1,hx,-h/2,h/2);%MEf0=subs(MEf0,x2,p);
MEf1=int(hx*MMix1,hx,-h/2,h/2);%MEf1=subs(MEf1,x2,p);
MEf2=int(hx^2*MMix1,hx,-h/2,h/2);%MEf2=subs(MEf2,x2,p);

MEf0=double(MEf0);MEf1=double(MEf1);MEf2=double(MEf2);


end