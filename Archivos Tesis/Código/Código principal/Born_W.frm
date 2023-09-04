**Este programa servirá para calcular el elemento de matriz asociada al Born del W**

**Se definen los elementos que se van a utilizar**


Symbols d,[2pi], m, g, e,[sqrt2],N, mW,Mlepton, D6, s12,s13,s14,s23,s24, sw,cw, cl, cq, Vud, Vqq;
Indices a1,a2,mu,nu,b1,b2,i1,i2,j1,j2,rho,theta,alpha;
Vectors p,r,k,q,p1,p2,p3,q1,q3;
CFunctions Deltafun,DeltaAdj;
Functions vrtx, Ub,U,Vb,V,W,Q,quark,electron,antiup,up,antilepton,lepton,neutrino,Gamma, Lepton, ZEULT1,ZEULT2;




**Se define el diagrama**
Global [d(1)]=
    Vb(1,p1)*
    vrtx(1,W(nu),quark(b1),quark(b2))*
    U(1,p2);

Global [d(2)]=
    Ub(2,q3)*
    vrtx(2,W(mu),lepton,neutrino)*
    V(2,p1+p2-q3)*
    W(mu,nu,p1+p2);

Global [d+(1)]=
    Ub(1,p2)*
    vrtx(1,W(theta),quark(b1),quark(b2))*
    V(1,p1);

Global [d+(2)]=
    -Vb(2,p1+p2-q3)*
    vrtx(2,W(rho),lepton,neutrino)*
    U(2,q3)*
    W(theta,rho,p1+p2);
 


**Definición de las reglas de Feynman utilizadas en la definición de los diagramas**
**Vertices**
id vrtx(1,W(mu?),quark(a1?),quark(a2?))= (i_*[2pi]^4)*i_*g*1/(2*[sqrt2])*Deltafun(a1,a2)*g_(1,mu)*g6_(1)*Vud;
id vrtx(2,W(mu?),lepton,neutrino)=(i_*[2pi]^4)*i_*g*1/(2*[sqrt2])*g_(2,mu)*g6_(2);

**Propagadores**
id W(mu?,nu?,p?)=(-i_/[2pi]^4)*d_(mu,nu)/(p^2+mW^2);
.sort



**Se hace explicito el cálculo que se quiere que Form haga con el diagrama definido**
Global [BB*] = (V(1,p1)*[d(1)]*[d+(1)]*ZEULT1*U(2,q3)*[d(2)]*[d+(2)]*ZEULT2);
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
id 1/(mW^2 + p1*p1+ 2*p1*p2 + p2*p2)=1/(mW^2-s12);
id 1/(mW^2-s12)=D6;
.sort



**Operaciones necesarias antes de realizar la traza de la linea férmionica (2)**
id d_(mu?,nu?)*d_(nu?,theta?)=d_(mu,theta);
id U(2,q3)*ZEULT2 = 1;
id U(2,q3)*Ub(2,q3) = -i_*g_(2,q3)+Mlepton;
id V(2,p?)*Vb(2,p?) = -i_*g_(2,p);


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
id [sqrt2]^-4=1/4;
id s12  =  Mlepton^2 - s13 - s23;


**Printeamos el resultado y se termina el programa**
Print [BB*];
.end

