

**Diagramas del Z**

**Emisión desde el lepton**

Global[d1+(1)]=
	-Ub(1,p2)*
	vrtx(1,antiup(b1),up(b2),Z(rho))*
	V(1,p1);

Global[d1+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,e,e,Z(beta))*
	Lepton(2,(q3+q1))*
	vrtx(2,Gamma(epsilon), lepton, lepton)*
	U(2,q3)*
	Z(p1+p2,beta,rho)*
	epsG(q1,epsilon);


**Emisión desde el anti-lepton**

Global[d2+(1)]=
	-Ub(1,p2)*
	vrtx(1,antiup(b1),up(b2),Z(rho))*
	V(1,p1);

Global[d2+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,Gamma(epsilon), antilepton, antilepton)*
	Lepton(2,-(p1+p2-q3))*
	vrtx(2,e,e,Z(beta))*
	U(2,q3)*
	Z(p1+p2,beta,rho)*
	epsG(q1,epsilon);


**Emisión desde el anti-quark**

Global [d3+(1)]=
	-Ub(1,p2)*
	vrtx(1, antiup(l1),up(b2),Z(rho))*
	Q(1,l2,l1,-(p1-q1))*
	vrtx(1,Gamma(epsilon),antiup(b1),antiup(l2))*
	V(1,p1);

Global [d3+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,e,e,Z(beta))*
	U(2,q3)*
	Z(p1+p2-q1,beta,rho)*
	epsG(q1,epsilon);


**Emisión desde el quark**

Global [d4+(1)]=
	-Ub(1,p2)*
	vrtx(1,Gamma(epsilon),up(b2),up(l2))*
	Q(1,l2,l1,p2-q1)*
	vrtx(1, antiup(b1),up(l1),Z(rho))*
	V(1,p1);

Global [d4+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,e,e,Z(beta))*
	U(2,q3)*
	Z(p1+p2-q1,beta,rho)*
	epsG(q1,epsilon);



**Diagramas del Gamma**

**Emisión desde el lepton**

Global[d13+(1)]=
	-Ub(1,p2)*
	vrtx(1,antiup(b1),up(b2),A(rho))*
	V(1,p1);

Global[d13+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,e,e,A(beta))*
	Lepton(2,q3+q1)*
	vrtx(2,Gamma(epsilon), lepton, lepton)*
	U(2,q3)*
	A(p1+p2,beta,rho)*
	epsG(q1,epsilon);


**Emisión desde el anti-leptón**

Global[d14+(1)]=
	-Ub(1,p2)*
	vrtx(1,antiup(b1),up(b2),A(rho))*
	V(1,p1);

Global[d14+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,Gamma(epsilon), antilepton, antilepton)*
	Lepton(2,-(p1+p2-q3))*
	vrtx(2,e,e,A(beta))*
	U(2,q3)*
	A(p1+p2,beta,rho)*
	epsG(q1,epsilon);


**Emisión desde el anti-quark**

Global [d15+(1)]=
	-Ub(1,p2)*
	vrtx(1, antiup(l1),up(b2),A(rho))*
	Q(1,l2,l1,-(p1-q1))*
	vrtx(1,Gamma(epsilon),antiup(b1),antiup(l2))*
	V(1,p1);

Global [d15+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,e,e,A(beta))*
	U(2,q3)*
	A(p1+p2-q1,beta,rho)*
	epsG(q1,epsilon);


**Emisión desde el quark**

Global [d16+(1)]=
	-Ub(1,p2)*
	vrtx(1,Gamma(epsilon),up(b2),up(l2))*
	Q(1,l2,l1,p2-q1)*
	vrtx(1, antiup(b1),up(l1),A(rho))*
	V(1,p1);

Global [d16+(2)]=
	Vb(2,p1+p2-q1-q3)*
	vrtx(2,e,e,A(beta))*
	U(2,q3)*
	A(p1+p2-q1,beta,rho)*
	epsG(q1,epsilon);























