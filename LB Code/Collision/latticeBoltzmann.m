function [f,omega,u,rho] = latticeBoltzmann(Q,D,w,e,dt,tau,c,rho)

%Declare your discrete density variables
f = sym("f",[Q 1]);

%Declare your weights
%w = sym("w",[Q 1]);

%Declare numeric parameters
%syms tau dt c;

%Declare fluid variables
%syms u rho;

%Assume incompressible flow with unity constant density
%rho = sum(f);

%Calculate the velocity
u = sym(zeros(D,1));
for d = 1:D
    for q = 1:Q
        u(d) = u(d)+ (c/rho)*(e(q,d).*f(q));
    end
end

u = expand(u);
u = collect(u,f);

%Calculate the equilibrium function
s = sym(zeros(Q,1));
feq = s;
eiu = s;
eiu2 = s;   
u2 = sym(0);

for q = 1:Q
    for d = 1:D
    eiu(q) = eiu(q)+sum(e(q,d).*u(d),1);
    end
end

eiu2 = eiu.^2;
u2 = sum(u.*u);

for q = 1:Q
s(q) = w(q)*((3/c)*eiu(q)+(9/(2*c^2))*(eiu2(q))-(3/(2*c^2))*u2);
feq(q) = w(q)*rho+rho*s(q);
end

%Calculate the collision term driving the differential equation
omega = -(dt/tau).*(f-feq);

rho = sum(f);

end