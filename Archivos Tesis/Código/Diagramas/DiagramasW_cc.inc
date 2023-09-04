   
**Diagramas del W**

**Emisión desde el lepton**

Global [d1+(1)]=
    Ub(1,p2)*
    vrtx(1,W(theta),quark(b1),quark(b2))*
    V(1,p1);

Global [d1+(2)]=
    -Vb(2,p1+p2-q1-q3)*
    vrtx(2,W(rho),lepton,neutrino)*
    Lepton(2,(q3+q1),Mlepton)*
    vrtx(2,Gamma(phi),lepton,lepton)*
    U(2,q3)*
    W(rho,theta,p1+p2)*
    epsG(q1,phi);


**Emisión desde el W**

Global [d2+(1)]=
    Ub(1,p2)*
    vrtx(1,W(theta),quark(b1),quark(b2))*
    V(1,p1);

Global [d2+(2)]=
    -Vb(2,p1+p2-q1-q3)*
    vrtx(2,W(rho),lepton,neutrino)*
    U(2,q3)*
    W(rho,epsilon2,p1+p2-q1)*
    vrtx(W(p1+p2,epsilon1),W(-(p1+p2-q1),epsilon2),Gamma(-q1,phi))*
    epsG(q1,phi)*
    W(epsilon1,theta,p1+p2);

**Emisión desde anti-quark**

Global [d3+(1)]=
    Ub(1,p2)*
    vrtx(1,W(theta),quark(l1),quark(b2))*
    Q(1,l1,l2,(p1-q1))*
    vrtx(1,Gamma(phi),up(b1),up(l2))*
    V(1,p1);    
    
Global [d3+(2)]=
    -Vb(2,p1+p2-q1-q3)*
    vrtx(2,W(rho),lepton,neutrino)*
    U(2,q3)*
    W(rho,theta,p1+p2-q1)*
    epsG(q1,phi);


**Emisión desde quark** 
  
 Global [d4+(1)]=
    Ub(1,p2)*
    vrtx(1,Gamma(phi),down(l2), down(b2))*
    Q(1,l1,l2,(p2-q1))*
    vrtx(1,W(theta),quark(l1),quark(b1))*
    V(1,p1);

    Global [d4+(2)]=
    -Vb(2,p1+p2-q1-q3)*
    vrtx(2,W(rho),lepton,neutrino)*
    U(2,q3)*
    W(rho,theta,p1+p2-q1)*
    epsG(q1,phi);

**Ghost**

Global [d5+(1)]=
    Ub(1,p2)*
    vrtx(1,W(phi),quark(b1),quark(b2))*
    V(1,p1);
    
Global [d5+(2)]=
    -Vb(2,p1+p2-q1-q3)*
    vrtx(2,ghost,lepton(-1),neutrino)*
    U(2,q3)*
    W(phi,theta,p1+p2)*
    vrtx(W(theta),ghost(-1),Gamma(rho))*
    Ghost(p1+p2-q1)*
    epsG(q1,rho);


