
**Diagramas conjugados del W**

**Emisión desde el lepton**

Global [d1(1)]=
    Vb(1,p1)*
    vrtx(1,W(nu),quark(b1),quark(b2))*
    U(1,p2);
    
    
Global [d1(2)]=
    Ub(2,q3)*
    vrtx(2,Gamma(alpha),lepton,lepton)*
    Lepton(2,(q3+q1),Mlepton)*
    vrtx(2,W(mu),lepton,neutrino)*
    V(2,p1+p2-q1-q3)*
    W(mu,nu,p1+p2)*
    epsG(q1,alpha);


**Emisión desde el W**

Global [d2(1)]=
    Vb(1,p1)*
    vrtx(1,W(nu),quark(b1),quark(b2))*
    U(1,p2);

Global [d2(2)]=
    Ub(2,q3)*
    vrtx(2,W(mu),lepton,neutrino)*
    V(2,p1+p2-q1-q3)*
    W(beta1,nu,p1+p2)*
    vrtx(W(p1+p2,beta1),W(-(p1+p2-q1),beta2),Gamma(-q1,alpha))*
    W(mu,beta2,p1+p2-q1)*
    epsG(q1,alpha);


**Emisión desde anti-quark**

Global [d3(1)]=
    Vb(1,p1)*
    vrtx(1,Gamma(alpha),up(b1), up(i2))*
    Q(1,i1,i2,(p1-q1))*
    vrtx(1,W(nu),quark(i1),quark(b2))*
    U(1,p2);

Global [d3(2)]=
    Ub(2,q3)*
    vrtx(2,W(mu),lepton,neutrino)*
    V(2,p1+p2-q1-q3)*
    W(mu,nu,p1+p2-q1)*
    epsG(q1,alpha);


**Emisión desde quark**

Global [d4(1)]=
    Vb(1,p1)*
    vrtx(1,W(nu),quark(b1),quark(i1))*
    Q(1,i1,i2,(p2-q1))*
    vrtx(1,Gamma(alpha),down(i2), down(b2))*
    U(1,p2);

Global [d4(2)]=
    Ub(2,q3)*
    vrtx(2,W(mu),lepton,neutrino)*
    V(2,p1+p2-q1-q3)*
    W(mu,nu,p1+p2-q1)*
    epsG(q1,alpha);


**Ghost**

Global [d5(1)]=
    Vb(1,p1)*
    vrtx(1,W(mu),quark(b1),quark(b2))*
    U(1,p2);

    
Global [d5(2)]=
    Ub(2,q3)*
    vrtx(2,ghost,lepton(1),neutrino)*
    V(2,p1+p2-q1-q3)*
    W(mu,nu,p1+p2)*
    vrtx(W(nu),ghost(1),Gamma(alpha))*
    Ghost(p1+p2-q1)*
    epsG(q1,alpha);


