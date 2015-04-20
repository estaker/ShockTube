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

%Free stream conditions (both sides)
u_inf = 0;
v_inf = 0;

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
gridhistory = zeros(n,m,ts);


%Initialize arrays
%Conservative State (old and new)
u = zeros(m,n,4);
u_old = zeros(m,n,4);
%Primitive State [rho u v P]
v = zeros(m,n,4);
%X Fluxes
A = zeros(m,n,4);
Ae = zeros(m,n,4);
Aw = zeros(m,n,4);
%Y Fluxes
B = zeros(m,n,4);
Bn = zeros(m,n,4);
Bs = zeros(m,n,4);
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
e = rand*.025;
gridhistory(:,:,1) = grid;
CFL = 1.8;
iu = zeros(size(u,1),size(u,2),size(u,3));

for i = 1:ts
    %Find Time Step
    c = sqrt(gamma*v(:,:,4)./v(:,:,1));
    dt = CFL ./ (abs(v(:,:,2))/dx + abs(v(:,:,3))/dy + c.*sqrt(1./dx.^2 + 1./dy.^2));
    global_dt = min(min(dt));

    %boundary condition treatment
    %LHS
    v(1,1,:) = [rhs_rho_inf 0 0 rhs_P_inf];
    %RHS
    v(end,1,:) = [lhs_rho_inf 0 0 lhs_P_inf];
    %boundary flux vectors
    u(1,1,:) = CFP(v(1,1,:));
    u(1,end,:) = CFP(v(1,end,:));
    A(1,1,:) = AFlux(v(1,1,:));
    B(1,1,:) = BFlux(v(1,1,:));
    A(end,1,:) = AFlux(v(end,1,:));
    B(end,1,:) = BFlux(v(end,1,:));
    
    %Start Solving for internal points
    Ae = .5*( A(1:end-2,:,:) + A(2:end-1,:,:));
    Aw = .5*( A(3:end,:,:) + A(2:end-1,:,:));
    
    %Settle for first order in B for now since we only have 2 points...
    
    %Calculate Artificial Dissapation
    ux = v(:,:,2);
    De = e*( .5*(ux(1:end-2,:) + ux(2:end-1,:)) + .5*(sos(v(1:end-2,:,:)) + sos(v(2:end-1,:,:)))) .* (u(2:end-1,:,:) - u(1:end-2,:,:));
    Dw = e*( .5*(ux(2:end-1,:) + ux(3:end,:)) + .5*(sos(v(2:end-1,:,:)) + sos(v(3:end,:,:)))) .* (u(3:end,:,:) - u(2:end-1,:,:));
    
    vy = v(:,:,3);
    %first order dispersion for n/s direction
    %Call dispersion in the y direction 0 for now, since that's what it'll
    %be
    
  
    
    
    
   
    
    
    
end
