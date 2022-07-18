
I:
parameters
	real rho = 998.0;
	real l = 1.0;
	real d = 0.01;
	real pi = 3.14159;
equations
    state = int(p.e);
    p.f = state / (rho*l/(pi*(d^2)/4));

R:
parameters
	real mu = 1e-3;
	real d = 0.01;
	real pi = 3.14159;
	real l = 1.0;
equations
	p.e = (128*mu*l/(pi*d^4)) * p.f;

C:
parameters
	real l = 1.0;
	real d = 0.01;
	real pi = 3.14159;
	real B = 2.2e9;
equations
	state = int(p.f);
	p.e = state*(((pi*(d^2)/4)*l)/B);


Se:
parameters
	real effort = 1.0;
	real rho = 998.0;
	real g = 9.81;
variables
	real flow;
equations
	p.e = effort*rho*g;
	flow = p.f;
