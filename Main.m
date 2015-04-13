%{
Elizabeth Staker
AE 516

Main Script for Shock-tube problem.
Objective: Verify Euler solver using transient solutions for a 1D shock
tube problem and a transient 2D disturbance propagation using both a
central-difference and upwind flux scheme


%}
%Flux calculation Type
Types = {'Jameson','Van Leer'};
selection = 1;
type = Types{selection};
%Pipe Conditions
length = 1; %m
height = .2; %m

%LHS initial conditions
rhs_rho_inf = 1.2; %kg/m3
rhs_P_inf = 100000; %Pa (1 bar)

%RHS initial conditions
lhs_rho_inf = 2.4; %kg/m^3
lhs_P_inf = 200000; %Pa (2 bar)

%Boundary Conditions 
%invicid walls (TODO)

%gas constants
gamma = 1.4;

%Grid Initialization
n = 2;
m = 100; 
grid = zeros(n,m);
dx = length/(m-1);
dy = dx;

%Number of time steps 
ts = 100;

%Initialize arrays
%Conservative State (old and new)
u = zeros(m,n,4);
u_old = zeros(m,n,4);
%Primitive State [rho u v P]
v = zeros(m,n,4);
%X Fluxes
A = zeros(m,n,4);
%Y Fluxes
B = zeros(m,n,4);
%Residuals
R = zeros(m,n,4);
%artificial dissapation
D = zeros(m,n,4);

%Establish Initial Conditions
v(1:m/2,:,1) = lhs_rho_inf;
v(m/2+1:end,:,1) = rhs_rho_inf;
v(1:m/2,:,4) = lhs_P_inf;
v(m/2+1:end,:,4) = rhs_P_inf;
v(:,:,2) = 0;
v(:,:,3) = 0;

u = CFP(v);

%get time step
CFL=2*sqrt(2) - eps;
CFL_safe = 1.8;
c = sqrt(gamma*v(:,:,4)./v(:,:,1));
dt = CFL ./ (abs(v(:,:,2))/dx + abs(v(:,:,3))/dy + c.*sqrt(1./dx.^2 + 1./dy.^2));

global_dt = min(min(dt));

