%% BDF2 experiment
clear all; clc; rng(42);
addpath('functions');
addpath('data');
h = 0.001; %time step
phi = @(z) (1/h)*((z.^2)/2 - 2*z + 3/2); % conformal map of BDF2 method
dphi = @(z) (1/h)*(z-2); % derivative 
phiinv_m = @(z) 2-sqrt((2*h)*z+1); % inverse of phi

% The clamped beam example can be found in
% https://modelreduction.org/morwiki/Clamped_Beam
system = load("beam.mat");
A = (system.A);
n = size(system.A,1);
B = sparse(system.B); 
C = sparse(system.C);

figure()
unitcircle = exp(1i*linspace(0,2*pi,1000));
mapped_circ = phi(unitcircle);
plot(real(mapped_circ),imag(mapped_circ),'r-'); hold on
plot(real(eigs(A,n)), imag(eigs(A,n)), 'b.');
legend('$\varphi(\partial D)$','$\Lambda(A)$','interpreter','latex');

I = speye(n);
G = @(s) C*((s*I-A)\B);
dG = @(s) -(C/(s*I-A))*((s*I-A)\B);


% algorithm 1
r = 10;
init = 0.1*randn(r,1)+0.1i*randn(r,1); init = (init./abs(init)).*rand(r,1);
H = @(s) G(phi(s));
dH = @(s) dG(phi(s)).*dphi(s);
[Er,Ar,Br,Cr,Dr,~,~] = algorithm1(H,dH,r,init,1000,1e-6);

% data for AAA
dataAAA = 500;
Z = [exp(1i*pi*logspace(-7,0,dataAAA))]; Z = [Z,conj(Z)];
Z2 = [phi(exp(1i*pi*logspace(-7,0,dataAAA)))]; Z2 = [Z2,conj(Z2)];
Rdata = zeros(size(Z));
R2data = zeros(size(Z2));
for i = 1:size(Z,2)
    R2data(i) = G(Z2(i));
    Rdata(i) = H(Z(i));
end

% Apply AAA to unit disk
[R,pol,res,zer] = aaa(Rdata,Z,'degree',r,'lawson',0);
[R3,pol3,res3,zer3] = aaa(Rdata,Z);

% Apply AAA to boomerang
[R2,pol2,~,~] = aaa(R2data,Z2,'degree',r,'lawson',0);


%%
% evaluation on unit circle
w = exp(2i*pi*logspace(-7,0,1000));
Gr = @(s) Cr*((s*Er-Ar)\Br)+Dr;
R2_ = @(s) R2(phi(s));
G_eval = zeros(1,size(w,2));
Gr_eval = zeros(1,size(w,2));
RD = zeros(1,size(w,2));
RD3 = zeros(1,size(w,2));
RA = zeros(1,size(w,2));
for i = 1:1:size(w,2)
    Gr_eval(i) = Gr(w(i));
    G_eval(i) = H(w(i));
    RD(i) = R(w(i));
    RD3(i) = R3(w(i));
    RA(i) = R2_(w(i));
end

figure()
loglog(2*pi*logspace(-7,0,1000),abs(RD),'-.'); hold on; 
loglog(2*pi*logspace(-7,0,1000),abs(RA),':');
loglog(2*pi*logspace(-7,0,1000),abs(RD3),':');
loglog(2*pi*logspace(-7,0,1000),abs(Gr_eval),'--');
loglog(2*pi*logspace(-7,0,1000),abs(G_eval),'-');
legend('$R_D$','$R_A$','$R_D^L$','$G_r$','$G$','interpreter','latex');
title('Evaluation on the unit circle');
hold off


% error plot
figure()
loglog(2*pi*logspace(-7,0,1000),abs(G_eval-RD),'-.'); hold on
loglog(2*pi*logspace(-7,0,1000),abs(G_eval-RA),':');
loglog(2*pi*logspace(-7,0,1000),abs(G_eval-RD3),':');
loglog(2*pi*logspace(-7,0,1000),abs(G_eval-Gr_eval),'--');
legend('$R_D$','$R_A$','$R_D^L$','$G_r$','interpreter','latex');
hold off

%%
% Discrete dynamical system
time = 3e4;
x = zeros(n,time);
xr = zeros(r,time);
u = zeros(1,time); 
u(10) = 1; % discrete-time impulse

% computing the discrete-time dynamics
for k = 3:1:time
    x(:,k) = (I-(2*h/3)*A)\(((4/3)*x(:,k-1) - (1/3)*x(:,k-2)) + B*(2*h/3)*u(k));
    xr(:,k) = (Ar\Er)*xr(:,k-1) - (Ar\Br)*u(k);
end
y = C*x;
yr = Cr*xr + Dr*u;

% computing the l-infinity error bound
Gr = @(s) Cr*((s*Er-Ar)\(Br))+Dr;
funerror = @(z) abs((H(exp(1i*z))-Gr(exp(1i*z)))).^2;
H2D_bound = sqrt((1/(2*pi))*integral(funerror,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true));

figure()
subplot(1,2,1)
skip = 5;
plot(1:skip:time,y(1:skip:end),'k'); hold on
plot(1:skip:time,yr(1:skip:end),'r--');
legend('$y_k$','$\widehat{y}_k$','interpreter','latex');
title('Discrete-time dynamics');
subplot(1,2,2)
semilogy(1:1:time,abs(y-yr),'r--'); hold on
semilogy(1:1:time,ones(1,time)*H2D_bound,'b:');
legend('$|y_k-\widehat{y}_k|$','boundary','interpreter','latex');
title('Error and bound');

%% H2A norm
fomfom = @(z) abs(H(exp(1i*z))).^2;
H2fomfom = sqrt((1/(2*pi))*integral(fomfom,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true));

%AAA on disk
[R3,~,~,~] = aaa(Rdata,Z);
funerror0 = @(z) abs((H(exp(1i*z))-R3(exp(1i*z)))).^2;
H2D_AAA3 = sqrt((1/(2*pi))*integral(funerror0,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true))./H2fomfom;


r_range = 4:4:40;
H2D_IRKA = zeros(size(r_range));
H2D_AAA = zeros(size(r_range));
H2D_AAA2 = zeros(size(r_range));
sigma = [];
for i = 1:1:size(r_range,2)
    r = r_range(i)
    opts = size(sigma,1);
    init = 0.1*randn(r-opts,1)+0.1i*randn(r-opts,1); init = (init./abs(init)).*rand(r-opts,1); init = [sigma; init];
    [Er,Ar,Br,Cr,Dr,sigma,conv] = algorithm1(H,dH,r,init,1000,1e-6);
    while conv == 0
        init = 0.1*randn(r,1)+0.1i*randn(r,1); init = (init./abs(init)).*rand(r,1);
        [Er,Ar,Br,Cr,Dr,sigma,conv] = algorithm1(H,dH,r,init,1000,1e-6);
    end
    Gr = @(s) Cr*((s*Er-Ar)\Br)+Dr;
    funerror = @(z) abs((H(exp(1i*z))-Gr(exp(1i*z)))).^2;
    H2D_IRKA(i) = sqrt((1/(2*pi))*integral(funerror,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true))./H2fomfom;

    %AAA on disk
    [R,~,~,~] = aaa(Rdata,Z,'degree',r,'lawson',0);
    funerror2 = @(z) abs((H(exp(1i*z))-R(exp(1i*z)))).^2;
    H2D_AAA(i) = sqrt((1/(2*pi))*integral(funerror2,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true))./H2fomfom;

    %AAA on boomerang
    [R2,~,~,~] = aaa(R2data,Z2,'degree',r,'lawson',0);
    funerror3 = @(z) abs((H(exp(1i*z))-R2(phi(exp(1i*z))))).^2;
    H2D_AAA2(i) = sqrt((1/(2*pi))*integral(funerror3,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true))./H2fomfom;
end

figure()
semilogy(r_range,H2D_IRKA,'r-x'); hold on
semilogy(r_range,H2D_AAA,'b-o');
semilogy(r_range,H2D_AAA2,'k-+');
semilogy(r_range,ones(size(r_range))*H2D_AAA3,'g-x');
legend('Algorithm 1', 'AAA on $D$', 'AAA on $A$', 'AAA on $D$ Lawson','Interpreter','latex','FontSize',14);
xlabel('r');
ylabel('$H_2(A)$ error','interpreter','latex');
hold off

