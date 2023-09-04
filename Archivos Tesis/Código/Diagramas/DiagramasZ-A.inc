
**Diagramas del Z**

**Emisión desde el lepton**

Global [d1(1)]=	Vb(1,p1)*
	vrtx(1,antiup(b1),up(b2),Z(nu))*
	U(1,p2);

Global[d1(2)]=
	Ub(2,q3)*
	vrtx(2,Gamma(alpha), lepton, lepton)*
	Lepton(2,(q3+q1))*
	vrtx(2,e,e,Z(mu))*
	V(2,p1+p2-q1-q3)*
	Z(p1+p2,mu,nu)*
	epsG(q1,alpha);
	

**Emisión desde el anti-lepton**

Global[d2(1)]=
	Vb(1,p1)*
	vrtx(1,antiup(b1),up(b2),Z(nu))*
	U(1,p2);

Global[d2(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,Z(mu))*
	Lepton(2,-(p1+p2-q3))*
	vrtx(2,Gamma(alpha), antilepton, antilepton)*
	V(2,p1+p2-q1-q3)*
	Z(p1+p2,mu,nu)*
	epsG(q1,alpha);


**Emisión desde el anti-quark**

Global [d3(1)]=
	Vb(1,p1)*
	vrtx(1,Gamma(alpha),antiup(b1),antiup(i2))*
	Q(1,i2,i1,-(p1-q1))*
	vrtx(1, antiup(i1),up(b2),Z(nu))*
	U(1,p2);

Global [d3(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,Z(mu))*
	V(2,p1+p2-q1-q3)*
	Z(p1+p2-q1,mu,nu)*
	epsG(q1,alpha);


**Emisión desde el quark**

Global [d4(1)]=
	Vb(1,p1)*
	vrtx(1, antiup(b1),up(i1),Z(nu))*
	Q(1,i2,i1,p2-q1)*
	vrtx(1,Gamma(alpha),up(b2),up(i2))*
	U(1,p2);

Global [d4(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,Z(mu))*
	V(2,p1+p2-q1-q3)*
	Z(p1+p2-q1,mu,nu)*
	epsG(q1,alpha);



**Diagramas del Gamma**

**Emisión desde el lepton**

Global [d13(1)]=
	Vb(1,p1)*
	vrtx(1,antiup(b1),up(b2),A(nu))*
	U(1,p2);

Global[d13(2)]=
	Ub(2,q3)*
	vrtx(2,Gamma(alpha), lepton, lepton)*
	Lepton(2,q3+q1)*
	vrtx(2,e,e,A(mu))*
	V(2,p1+p2-q1-q3)*
	A(p1+p2,mu,nu)*
	epsG(q1,alpha);


**Emisión desde el anti-leptón**

Global[d14(1)]=
	Vb(1,p1)*
	vrtx(1,antiup(b1),up(b2),A(nu))*
	U(1,p2);

Global[d14(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,A(mu))*
	Lepton(2,-(p1+p2-q3))*
	vrtx(2,Gamma(alpha), antilepton, antilepton)*
	V(2,p1+p2-q1-q3)*
	A(p1+p2,mu,nu)*
	epsG(q1,alpha);


**Emisión desde el anti-quark**

Global [d15(1)]=
	Vb(1,p1)*
	vrtx(1,Gamma(alpha),antiup(b1),antiup(i2))*
	Q(1,i2,i1,-(p1-q1))*
	vrtx(1, antiup(i1),up(b2),A(nu))*
	U(1,p2);

Global [d15(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,A(mu))*
	V(2,p1+p2-q1-q3)*
	A(p1+p2-q1,mu,nu)*
	epsG(q1,alpha);


**Emisión desde el quark**

Global [d16(1)]=
	Vb(1,p1)*
	vrtx(1, antiup(b1),up(i1),A(nu))*
	Q(1,i2,i1,p2-q1)*
	vrtx(1,Gamma(alpha),up(b2),up(i2))*
	U(1,p2);

Global [d16(2)]=
	Ub(2,q3)*
	vrtx(2,e,e,A(mu))*
	V(2,p1+p2-q1-q3)*
	A(p1+p2-q1,mu,nu)*
	epsG(q1,alpha);


