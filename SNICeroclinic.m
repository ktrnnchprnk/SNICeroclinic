% Non-central SNICeroclinic in plynomial model
% 10/2024
function xdot=SNICeroclinic(~,var,p)
% Initialise the system
xdot=zeros(2,1);
% Parameters
a=p.a;
b=p.b;
c=p.c;

eps=p.eps;

x=var(1);
y=var(2);

func_y = y.^3-3*y;
func_x= a.*x.^2+b.*x+c;
% ODE system
xdot(1) = eps.*(func_x-y);
xdot(2) = x-func_y;