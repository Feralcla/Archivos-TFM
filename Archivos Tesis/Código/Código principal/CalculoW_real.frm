**Este programa servirá para calcular el elemento de matriz asociada al real del Z/Gamma**

**Se definen los elementos que se van a utilizar**
Symbols d,[2pi],m,Qc, g,Vqq,e,gs,[sqrt2],cl,V1,A1,N, Np,mW,mG,Mlepton,daga,D1,D2,D3,D4,D5,s12,s13,s14,s23,s24,cq,sw, Vud;
Indices a1,a2,j1,mu,nu,b1,b2,i1,i2,i3,i4,l1,l2,rho,,epsilon,epsilon1,epsilon2,phi,beta1,beta2,theta,alpha;
Vectors r,p,k,q,p1,p2,p3,q1,q3;
CFunctions DeltaAdj,epsG,Deltafun;
Functions Ub,U,Vb,V,Q,vrtx,quark,electron,up,down,ghost,Ghost,Gamma,lepton,Lepton,neutrino,W,ZEULT1,ZEULT2;




**Se llama a los diagramas**
#-
#include DiagramasW.inc
#include DiagramasW_cc.inc
.sort
#+



**Definición de las reglas de Feynman utilizadas en la definición de los diagramas**
**Vertices**
id vrtx(1,W(mu?),quark(a1?),quark(a2?))= (i_*[2pi]^4)*i_*g*1/(2*[sqrt2])*Deltafun(a1,a2)*g_(1,mu)*g6_(1)*Vud;
id vrtx(2,W(mu?),lepton,neutrino)=(i_*[2pi]^4)*i_*g*1/(2*[sqrt2])*g_(2,mu)*g6_(2);
id vrtx(W(q?,nu?),W(k?,rho?),Gamma(p?,mu?))=(i_*[2pi]^4)*(-1)*g*sw*(d_(mu,nu)*(p(rho)-q(rho))+d_(nu,rho)*(q(mu)-k(mu))+d_(rho,mu)*(k(nu)-p(nu)));
id vrtx(1,Gamma(alpha?),up(a1?),up(a2?))=(i_*[2pi]^4)*(-2/3)*i_*g*sw*Deltafun(a1,a2)*g_(1,alpha);
id vrtx(1,Gamma(alpha?),down(a1?), down(a2?))=(i_*[2pi]^4)*(-1/3)*i_*g*sw*Deltafun(a1,a2)*g_(1,alpha);
id vrtx(2,Gamma(alpha?),lepton,lepton)=(i_*[2pi]^4)*(-1)*i_*g*sw*g_(2,alpha);
id vrtx(W(mu?),ghost(daga?),Gamma(nu?))=(i_*[2pi]^4)*(-1)*(daga)*i_*g*sw*mW*d_(mu,nu);
id vrtx(2,ghost,lepton(daga?),neutrino)=(i_*[2pi]^4)*i_*g*1/(2*[sqrt2])*(Mlepton/mW)*(daga+g5_(2));



*Propagadores*
id Q(1,i1?,i2?,p?) = (-i_/[2pi]^4)*Deltafun(i1,i2)*(-i_*g_(1,p))/p^2;
id Lepton(2,p?,Mlepton?)=(-i_/[2pi]^4)*(-i_*g_(2,p)+Mlepton)/(p^2+Mlepton^2);
id W(mu?,nu?,p?)=(-i_/[2pi]^4)*d_(mu,nu)/(p^2+mW^2);
id Ghost(k?)=(-i_/[2pi]^4)*1/(k^2+mW^2);
.sort


**Interferencia**
Global [AA*]=(V(1,p1)*[d1(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d1+(2)]*ZEULT2)+
(V(1,p1)*[d2(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d2+(2)]*ZEULT2)+
(V(1,p1)*[d3(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d3+(2)]*ZEULT2)+
(V(1,p1)*[d4(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d4+(2)]*ZEULT2)+
(V(1,p1)*[d5(1)]*[d5+(1)]*ZEULT1*U(2,q3)*[d5(2)]*[d5+(2)]*ZEULT2)+
(V(1,p1)*[d1(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d2+(2)]*ZEULT2)+(V(1,p1)*[d2(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d1+(2)]*ZEULT2)
 +(V(1,p1)*[d1(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d3+(2)]*ZEULT2)+(V(1,p1)*[d3(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d1+(2)]*ZEULT2)
+(V(1,p1)*[d1(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d4+(2)]*ZEULT2)+(V(1,p1)*[d4(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d1+(2)]*ZEULT2)
 +(V(1,p1)*[d1(1)]*[d5+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d5+(2)]*ZEULT2)+(V(1,p1)*[d5(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d5(2)]*[d1+(2)]*ZEULT2)
+(V(1,p1)*[d2(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d3+(2)]*ZEULT2)+(V(1,p1)*[d3(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d2+(2)]*ZEULT2)
 +(V(1,p1)*[d2(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d4+(2)]*ZEULT2)+(V(1,p1)*[d4(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d2+(2)]*ZEULT2)
+(V(1,p1)*[d2(1)]*[d5+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d5+(2)]*ZEULT2)+(V(1,p1)*[d5(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d5(2)]*[d2+(2)]*ZEULT2)
 +(V(1,p1)*[d3(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d4+(2)]*ZEULT2)+(V(1,p1)*[d4(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d3+(2)]*ZEULT2)
 +(V(1,p1)*[d3(1)]*[d5+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d5+(2)]*ZEULT2)+(V(1,p1)*[d5(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d5(2)]*[d3+(2)]*ZEULT2)
 +(V(1,p1)*[d4(1)]*[d5+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d5+(2)]*ZEULT2)+(V(1,p1)*[d5(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d5(2)]*[d4+(2)]*ZEULT2);
.sort



**Esta parte del programa se encargará de tener en cuenta la estructura de color**
repeat ;
id once DeltaAdj(i1?,i2?)*DeltaAdj(i2?,j1?) = DeltaAdj(i1,j1) ;
id once DeltaAdj(i1?,i2?)*DeltaAdj(j1?,i2?) = DeltaAdj(i1,j1) ;
id once DeltaAdj(i2?,i1?)*DeltaAdj(i2?,j1?) = DeltaAdj(i1,j1) ;
endrepeat ;

repeat ;
id once Deltafun(i1?,i2?)*Deltafun(i2?,j1?) = Deltafun(i1,j1) ;
id once Deltafun(i1?,i2?)*Deltafun(j1?,i2?) = Deltafun(i1,j1) ;
id once Deltafun(i2?,i1?)*Deltafun(i2?,j1?) = Deltafun(i1,j1) ;
endrepeat ;

id Deltafun(i1?,i1?) = N ;
.sort






**Se definen los denominadores que aparecen en los cálculos de nuestros diagramas**
id 1/(mW^2+(p1+p2)^2)=1/(mW^2-s12);
id 1/(mW^2+(p1+p2-q1)^2)=1/(mW^2-s12-s14-s24);
id 1/(Mlepton^2 + q1*q1 + 2*q1*q3 + q3*q3)=1/(s12+s13+s14+s23+s24-Mlepton^2);
id 1/(p1+p2-q1)^2=1/(-s12-s14-s24);
id 1/(mW^2 + p1*p1 + 2*p1*p2 - 2*p1*q1 + p2*p2 - 2*p2*q1 + q1*q1)=1/(mW^2-s12-s14-s24);

id 1/(mW^2-s12)=D1;
id 1/(mW^2-s12-s14-s24)=D2;
id 1/(s12+s13+s14+s23+s24-Mlepton^2)=D3;
id 1/(-s12-s14-s24)=D4;
id 1/(mW^2-s12-s14-s24)=D5;

id 1/(p1-q1)^2=-1/s14;
.sort



**Operaciones necesarias antes de realizar la traza de la linea férmionica (2)**
id epsG(p?,mu?)*epsG(p?,nu?) = d_(mu,nu);
id d_(mu?,nu?)*d_(nu?,theta?)=d_(mu,theta);
id U(2,q3)*ZEULT2 = 1;
id U(2,q3)*Ub(2,q3) = -i_*g_(2,q3)+Mlepton;
id V(2,p1?)*Vb(2,p1?) = -i_*g_(2,p1);





**Operaciones necesarias antes de realizar la traza de la linea férmionica (1)**
id V(1,p1)*ZEULT1 = 1;
id V(1,p1)*Vb(1,p1) = -i_*g_(1,p1);
id U(1,p2)*Ub(1,p2) = -i_*g_(1,p2);
id d_(mu?,nu?)*d_(nu?,theta?)=d_(mu,theta);




**Órdenes para que Form cálcule la traza**
Trace4, 2;
Trace4, 1;
.sort
contract;



**Se definen los productos escalares (cinemática del problema)**id p1.p1 = 0;
id p1.p2 = -s12/2;
id p1.q1 = s14/2;
id p1.q3 = (s13-Mlepton^2)/2;
id p2.p2 = 0;
id p2.q1 = s24/2;
id p2.q3 = (s23-Mlepton^2)/2;
id q1.q1 = 0;
id q1.q3 = (s12+s13+s14+s23+s24-Mlepton^2)/2;
id q3.q3 = -Mlepton^2;
.sort

id 1/(p2*p2 - 2*p2*q1 + q1*q1)/(p2*p2 - 2*p2*q1 + q1*q1)=1/s24^2;
id 1/(p2*p2 - 2*p2*q1 + q1*q1)=-1/s24;


**Identidades para simplificar el resultado**
id d=4;
id N=3;
id [sqrt2]^-4=1/4;
id s14 = Qc^2-s12-s24;


**Printeamos el resultado y se termina el programa**
Print [AA*];
.end








