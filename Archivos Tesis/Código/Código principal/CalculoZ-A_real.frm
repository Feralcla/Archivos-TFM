**Este programa servirá para calcular el elemento de matriz asociada al real del Z/Gamma**

**Se definen los elementos que se van a utilizar**
Symbols d,[2pi], m,Qc,g,Vqq,e,gs,[sqrt2],cl,N,,mZ,mG,Mlepton,D1,D2,D3,D4,D5,s12,s13,s14,s23,s24,cq,sw,cw,CV,CA,cq;
Indices a1,a2,j1,mu,nu,m1,m2,b1,b2,i1,i2,i3,i4,l1,l2,i,j,rho,beta,theta,alpha,epsilon;
Vectors r,p,k,q,p1,p2,p3,q1,q3;
CFunctions epsG,Deltafun, DeltaAdj;
Functions Ub,U,Vb,V,A,Q,vrtx,quark,electron,antiup,up,antidown,down,Goldstone,Gold,Gamma,antilepton,lepton,Lepton,Z, ZEULT1,ZEULT2;


**Se llama a los diagramas**
nwrite statistics;
.global
#-
#include DiagramasZ-A.inc
#include DiagramasZ-A_cc.inc
.sort
#+



**Definición de las reglas de Feynman utilizadas en la definición de los diagramas**
**Aniquilación/creación del Z**
id vrtx(1,antiup(a1?),up(a2?),Z(nu?))=    (i_*[2pi]^4)*i_*g*(1)/(4*cw)*g_(1,nu)*(CV + CA*g5_(1))*Deltafun(a1,a2);
id vrtx(2,e,e,Z(nu?))=(i_*[2pi]^4)*i_*g*(-1)/(4*cw)*g_(2,nu)*(-1+4*sw^2-g5_(2));

**Aniquilación/creación del Gamma**
id vrtx(1,antiup(a1?),up(a2?),A(nu?))    =(i_*[2pi]^4)*(cq)*i_*g*sw*Deltafun(a1,a2)*g_(1,nu);
id vrtx(2,e,e,A(nu?))=(i_*[2pi]^4)*(-1)*i_*g*sw*g_(2,nu);

**Emisión del fotón**
id vrtx(1,Gamma(alpha?),up(a1?),up(a2?))=(i_*[2pi]^4)*(cq)*i_*g*sw*Deltafun(a1,a2)*g_(1,alpha);
id vrtx(1,Gamma(alpha?),antiup(a1?),antiup(a2?))=(i_*[2pi]^4)*(-cq)*i_*g*sw*Deltafun(a1,a2)*g_(1,alpha);
id vrtx(2,Gamma(alpha?),lepton,lepton)=(i_*[2pi]^4)*(-1)*i_*g*sw*g_(2,alpha);
id vrtx(2,Gamma(alpha?),antilepton,antilepton)=(i_*[2pi]^4)*(1)*i_*g*sw*g_(2,alpha);

**Propagadores**
id Lepton(2,p?)=(-i_/[2pi]^4)*(-i_*g_(2,p)+Mlepton)/(p^2+Mlepton^2);
id Q(1,i1?,i2?,p?) = (-i_/[2pi]^4)*Deltafun(i1,i2)*(-i_*g_(1,p))/p^2;
id Z(k?,mu?,nu?)=(-i_/[2pi]^4)*d_(mu,nu)/(k^2+mZ^2);
id A(k?,mu?,nu?)=(-i_/[2pi]^4)*d_(mu,nu)/(k^2);
.sort



**Interferencia**

Global[Ups-Leptones]=
(V(1,p1)*[d1(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d1+(2)]*ZEULT2)+
(V(1,p1)*[d2(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d2+(2)]*ZEULT2)+
(V(1,p1)*[d3(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d3+(2)]*ZEULT2)+
(V(1,p1)*[d4(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d4+(2)]*ZEULT2)+
(V(1,p1)*[d13(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d13+(2)]*ZEULT2)+
(V(1,p1)*[d14(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d14+(2)]*ZEULT2)+
(V(1,p1)*[d15(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d15+(2)]*ZEULT2)+
(V(1,p1)*[d16(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d16+(2)]*ZEULT2)-
(V(1,p1)*[d1(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d2+(2)]*ZEULT2)-(V(1,p1)*[d2(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d1+(2)]*ZEULT2)-
(V(1,p1)*[d1(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d3+(2)]*ZEULT2)-(V(1,p1)*[d3(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d1+(2)]*ZEULT2)+
(V(1,p1)*[d1(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d4+(2)]*ZEULT2)+(V(1,p1)*[d4(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d1+(2)]*ZEULT2)-
(V(1,p1)*[d1(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d13+(2)]*ZEULT2)-(V(1,p1)*[d13(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d1+(2)]*ZEULT2)+
(V(1,p1)*[d1(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d14+(2)]*ZEULT2)+(V(1,p1)*[d14(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d1+(2)]*ZEULT2)+
(V(1,p1)*[d1(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d15+(2)]*ZEULT2)+(V(1,p1)*[d15(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d1+(2)]*ZEULT2)-
(V(1,p1)*[d1(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d1(2)]*[d16+(2)]*ZEULT2)-(V(1,p1)*[d16(1)]*[d1+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d1+(2)]*ZEULT2)+
(V(1,p1)*[d2(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d3+(2)]*ZEULT2)+(V(1,p1)*[d3(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d2+(2)]*ZEULT2)-
(V(1,p1)*[d2(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d4+(2)]*ZEULT2)-(V(1,p1)*[d4(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d2+(2)]*ZEULT2)+
(V(1,p1)*[d2(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d13+(2)]*ZEULT2)+(V(1,p1)*[d13(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d2+(2)]*ZEULT2)-
(V(1,p1)*[d2(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d14+(2)]*ZEULT2)-(V(1,p1)*[d14(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d2+(2)]*ZEULT2)-
(V(1,p1)*[d2(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d15+(2)]*ZEULT2)-(V(1,p1)*[d15(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d2+(2)]*ZEULT2)+
(V(1,p1)*[d2(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d2(2)]*[d16+(2)]*ZEULT2)+(V(1,p1)*[d16(1)]*[d2+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d2+(2)]*ZEULT2)-
(V(1,p1)*[d3(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d4+(2)]*ZEULT2)-(V(1,p1)*[d4(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d3+(2)]*ZEULT2)+
(V(1,p1)*[d3(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d13+(2)]*ZEULT2)+(V(1,p1)*[d13(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d3+(2)]*ZEULT2)-
(V(1,p1)*[d3(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d14+(2)]*ZEULT2)-(V(1,p1)*[d14(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d3+(2)]*ZEULT2)-
(V(1,p1)*[d3(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d15+(2)]*ZEULT2)-(V(1,p1)*[d15(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d3+(2)]*ZEULT2)+
(V(1,p1)*[d3(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d3(2)]*[d16+(2)]*ZEULT2)+(V(1,p1)*[d16(1)]*[d3+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d3+(2)]*ZEULT2)-
(V(1,p1)*[d4(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d13+(2)]*ZEULT2)-(V(1,p1)*[d13(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d4+(2)]*ZEULT2)+
(V(1,p1)*[d4(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d14+(2)]*ZEULT2)+(V(1,p1)*[d14(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d4+(2)]*ZEULT2)+
(V(1,p1)*[d4(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d15+(2)]*ZEULT2)+(V(1,p1)*[d15(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d4+(2)]*ZEULT2)-
(V(1,p1)*[d4(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d4(2)]*[d16+(2)]*ZEULT2)-(V(1,p1)*[d16(1)]*[d4+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d4+(2)]*ZEULT2)-
(V(1,p1)*[d13(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d14+(2)]*ZEULT2)-(V(1,p1)*[d14(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d13+(2)]*ZEULT2)-
(V(1,p1)*[d13(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d15+(2)]*ZEULT2)-(V(1,p1)*[d15(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d13+(2)]*ZEULT2)+
(V(1,p1)*[d13(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d13(2)]*[d16+(2)]*ZEULT2)+(V(1,p1)*[d16(1)]*[d13+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d13+(2)]*ZEULT2)+
(V(1,p1)*[d14(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d15+(2)]*ZEULT2)+(V(1,p1)*[d15(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d14+(2)]*ZEULT2)-
(V(1,p1)*[d14(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d14(2)]*[d16+(2)]*ZEULT2)-(V(1,p1)*[d16(1)]*[d14+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d14+(2)]*ZEULT2)-
(V(1,p1)*[d15(1)]*[d16+(1)]*ZEULT1*U(2,q3)*[d15(2)]*[d16+(2)]*ZEULT2)-(V(1,p1)*[d16(1)]*[d15+(1)]*ZEULT1*U(2,q3)*[d16(2)]*[d15+(2)]*ZEULT2);
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
id 1/(mZ^2 + p1*p1 + 2*p1*p2 + p2*p2)=1/(mZ^2-s12);
id 1/(Mlepton^2 + p1*p1 + 2*p1*p2 - 2*p1*q3 + p2*p2 - 2*p2*q3 + q3*q3)=1/ (2*Mlepton^2-s12 - s13 - s23);
id 1/(Mlepton^2 + p1*p1 + 2*p1*p2 - 2*p1*q1 + p2*p2 - 2*p2*q1 + q1*q1)=1/(Mlepton^2 -s12 - s14 - s24);
id 1/(Mlepton^2 + q1*q1 + 2*q1*q3 + q3*q3)=1/(s12+s13+s14+s23+s24-2*Mlepton^2);
id 1/(mZ^2      + p1*p1 + 2*p1*p2 - 2*p1*q1 + p2*p2 - 2*p2*q1 + q1*q1)=1/(mZ^2  -s12 - s14 - s24);
id 1/(p1*p1 + 2*p1*p2 - 2*p1*q1 +p2*p2 - 2*p2*q1 + q1*q1)=1/(-s12 - s14 - s24);

id 1/(p2*p2 - 2*p2*q3 + q3*q3)=-1/s23;
id 1/(p2*p2 - 2*p2*q1 + q1*q1)=-1/s24;
id 1/(p1*p1 - 2*p1*q1 + q1*q1)=-1/s14;
id 1/(p1*p1 + 2*p1*p2 + p2*p2)=-1/s12;

id 1/(mZ^2-s12)=D1;
id 1/ (2*Mlepton^2-s12 - s13 - s23)=D2;
id 1/(s12+s13+s14+s23+s24-2*Mlepton^2)=D3;
id 1/(mZ^2  -s12 - s14 - s24)=D4;
id 1/(-s12 - s14 - s24)=D5;






**Operaciones necesarias antes de realizar la traza de la linea férmionica (2)**
id epsG(p?,mu?)*epsG(p?,nu?) = d_(mu,nu);
id d_(mu?,nu?)*d_(nu?,theta?)=d_(mu,theta);
id U(2,q3)*ZEULT2 = 1;
id U(2,q3)*Ub(2,q3) = -i_*g_(2,q3)+Mlepton;
id V(2,p1?)*Vb(2,p1?) = -i_*g_(2,p1)-Mlepton;



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



**Se definen los productos escalares (cinemática del problema)**
id p1.p1 = 0;
id p1.p2 = -s12/2;
id p1.q1 = s14/2;
id p1.q3 = (s13-Mlepton^2)/2;
id p2.p2 = 0;
id p2.q1 = s24/2;
id p2.q3 = (s23-Mlepton^2)/2;
id q1.q1 = 0;
id q1.q3 = (s12+s13+s14+s23+s24-2*Mlepton^2)/2;
id q3.q3 = -Mlepton^2;
.sort

**Identidades para simplificar el resultado**
id d=4;
id N=3;
id s14 = Qc^2-s12-s24;


**Printeamos el resultado y se termina el programa**
print[Ups-Leptones];
.end


