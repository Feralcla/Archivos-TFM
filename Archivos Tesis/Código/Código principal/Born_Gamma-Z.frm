**Este programa servirá para calcular el elemento de matriz asociada al Born del Z/Gamma**

**Se definen los elementos que se van a utilizar**

Symbols d,[2pi], m, g, e,[sqrt2],N, mZ,Mlepton, D1,D2,D3,D4,D5, s12,s13,s14,s23,s24, sw,cw, CV,CA,cq;
Indices a1,a2,mu,nu,b1,b2,i1,i2,j1,j2,rho,theta,alpha;
Vectors p,r,k,q,p1,p2,p3,q1,q3;
CFunctions Deltafun,DeltaAdj;
Functions vrtx, Ub,U,Vb,V, A,Q,quark,electron,antiup,up,Gamma,antilepton,lepton,Lepton,Z, ZEULT1,ZEULT2;




**Se definen los diagramas**

**Diagrama del Z**
Global [dZ(1)]=	Vb(1,p1)*
	vrtx(1,antiup(b1),up(b2),Z(nu))*
	U(1,p2);

Global[dZ(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,Z(mu))*
	V(2,p1+p2-q3)*
	Z(p1+p2,mu,nu);

Global [dZ+(1)]=
	Ub(1,p2)*
	vrtx(1,antiup(b1),up(b2),Z(rho))*
	V(1,p1);

Global[dZ+(2)]=
	-Vb(2,p1+p2-q3)*
	vrtx(2,e,e,Z(theta))*
	U(2,q3)*
	Z(p1+p2,theta,rho);


**Diagrama del Gamma**
Global [dG(1)]=	
	Vb(1,p1)*
	vrtx(1,antiup(b1),up(b2),A(nu))*
	U(1,p2);

Global[dG(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,A(mu))*
	V(2,p1+p2-q3)*
	A(p1+p2,mu,nu);

Global [dG+(1)]=
	Ub(1,p2)*
	vrtx(1,antiup(b1),up(b2),A(rho))*
	V(1,p1);

Global[dG+(2)]=
	-Vb(2,p1+p2-q3)*
	vrtx(2,e,e,A(theta))*
	U(2,q3)*
	A(p1+p2,theta,rho);



**Definición de las reglas de Feynman utilizadas en la definición de los diagramas**
**Aniquilación/creación del Z**
id vrtx(1,antiup(a1?),up(a2?),Z(nu?))=    (i_*[2pi]^4)*i_*g*(-1)/(4*cw)*g_(1,nu)*(CV+CA*g5_(1))*Deltafun(a1,a2);
id vrtx(2,e,e,Z(nu?))=(i_*[2pi]^4)*i_*g*(-1)/(4*cw)*g_(2,nu)*(-1+4*sw^2-g5_(2));

**Aniquilación/creación del Gamma**
id vrtx(1,antiup(a1?),up(a2?),A(nu?))    =(i_*[2pi]^4)*(cq)*i_*g*sw*Deltafun(a1,a2)*g_(1,nu);
id vrtx(2,e,e,A(nu?))=(i_*[2pi]^4)*(-1)*i_*g*sw*g_(2,nu);

**Propagadores**
id Z(k?,mu?,nu?)=(-i_/[2pi]^4)*d_(mu,nu)/(k^2+mZ^2);
id A(k?,mu?,nu?)=(-i_/[2pi]^4)*d_(mu,nu)/(k^2);
.sort



**Se hace explicito el cálculo que se quiere que Form haga con los diagramas definidos**

Global[BB*] = (V(1,p1)*[dZ(1)]*[dZ+(1)]*ZEULT1*U(2,q3)*[dZ(2)]*[dZ+(2)]*ZEULT2) +
	      (V(1,p1)*[dG(1)]*[dG+(1)]*ZEULT1*U(2,q3)*[dG(2)]*[dG+(2)]*ZEULT2)	+
              (V(1,p1)*[dZ(1)]*[dG+(1)]*ZEULT1*U(2,q3)*[dZ(2)]*[dG+(2)]*ZEULT2) + (V(1,p1)*[dG(1)]*[dZ+(1)]*ZEULT1*U(2,q3)*[dG(2)]*[dZ+(2)]*ZEULT2);		
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
id 1/(Mlepton^2 + q1*q1 + 2*q1*q3 + q3*q3)=1/(s12+s13+s14+s23+s24-2*Mlepton^2);
id 1/(mZ^2      + p1*p1 + 2*p1*p2 - 2*p1*q1 + p2*p2 - 2*p2*q1 + q1*q1)=1/(mZ^2  -s12 - s14 - s24);
id 1/(p1*p1 + 2*p1*p2 - 2*p1*q1 +p2*p2 - 2*p2*q1 + q1*q1)=1/(-s12 - s14 - s24);

id 1/(mZ^2-s12)=D1;
id 1/ (2*Mlepton^2-s12 - s13 - s23)=D2;
id 1/(s12+s13+s14+s23+s24-2*Mlepton^2)=D3;
id 1/(mZ^2  -s12 - s14 - s24)=D4;
id 1/(-s12 - s14 - s24)=D5;

id 1/(p2*p2 - 2*p2*q3 + q3*q3)=-1/s23;
id 1/(p1*p1 + 2*p1*p2 + p2*p2)=-1/s12;




**Operaciones necesarias antes de realizar la traza de la linea férmionica (2)**
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
id p1.q3 = (s13-Mlepton^2)/2;
id p2.p2 = 0;
id p2.q3 = (s23-Mlepton^2)/2;
id q3.q3 = -Mlepton^2;
.sort


**Identidades para simplificar el resultado**
id d=4;
id N=3;
id s12 = - s13 - s23 + 2*Mlepton^2;


**Printeamos el resultado y se termina el programa**
print[BB*];
.end


