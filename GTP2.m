function xdot=GTP2(~,x,p)
% initialise the system 
xdot=zeros(2,1);

%Parameters
beta     = p.beta;
b       = p.b;
gamma     = p.gamma;
GT      = p.GT;
ell0    = p.ell0;
phi     = p.phi;
phi2     =p.phi2;
Gh      = p.Gh;
epsilon = p.epsilon;
hilln   = p.hilln;
hillp   = p.hillp;
m=   p.m;

G=x(1);
L=x(2);

restlength=(ell0-phi*G.^hillp./(Gh.^hillp+G.^hillp));
restlength2=(ell0-phi2*G.^hillp./(Gh.^hillp+G.^hillp));
squashing=beta.*(L.^m)./(restlength.^m-L.^m);
% ODE system
xdot(1) =((b+squashing+gamma.*(G.^hilln)/(1+G.^hilln)).*(GT-G)-G);
xdot(2) = -epsilon.*(L-restlength2);

end
