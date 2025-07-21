%% Weird testing
clear all; clc; rng(42);
h = 1e-3; %time step
phi = @(z) (1/h)*((z.^4)/4 - (4/3)*(z.^3) + 3*(z.^2) - 4*z + 25/12); % conformal map BDF4 method (not A-stable)
dphi = @(z) (1/h)*((z.^3) - 4*(z.^2) + 6*z - 4);

% convection diffusion equation
% https://www.netlib.org/lyapack/guide.pdf
% https://morwiki.mpi-magdeburg.mpg.de/morwiki/Convection-Diffusion#cite_note-lyapack-1

A = fdm_2d_matrix(100, '100*x', '100*y', '0');
B = fdm_2d_vector(100, '.1<x<=.3');
C = fdm_2d_vector(100, '.7<x<=.9')';
n = size(A,1);

unitcircle = exp(1i*linspace(0,2*pi,1000));

I = speye(n);
G = @(s) C*((s*I-A)\B);
dG = @(s) -(C/(s*I-A))*(((s*I-A)\B));

% running algorithm 1
r = 10;
init = 0.1*randn(r,1)+0.1i*randn(r,1); init = (init./abs(init)).*rand(r,1);
H = @(s) G(phi(s));
dH = @(s) dG(phi(s)).*dphi(s);
[Er,Ar,Br,Cr,Dr,sigma,~] = algorithm1(H,dH,r,init,1000,1e-6);

% plotting reduced order poles
figure()
eigsAr = eig(Ar,Er);
plot(real(eigsAr),imag(eigsAr),'ko'); hold on
plot(real(unitcircle),imag(unitcircle),'k');
plot(real(sigma),imag(sigma),'k.');
legend('$\Lambda(\mathbf{A}_r)$','$D$','$\sigma$', 'interpreter', 'latex', 'FontSize',15, 'location', 'northwest')
axis equal


% plot poles FOM and A
figure()
partialA = phi(unitcircle);
eigenvalues = load('eigenvaluesA'); eigsA = eigenvalues.eigsA;
plot(real(eigsA),imag(eigsA),'k.'); hold on
plot(real(partialA),imag(partialA),'k','linewidth',2);
legend('$\Lambda(\mathbf{A})$', 'interpreter', 'latex', 'FontSize',15, 'location', 'northwest')
xlim([-8.5e4,1.5e4]);
ylim([-7e3,7e3]);
xlabel('Re($z$)','interpreter','latex');
ylabel('Im($z$)','interpreter','latex');


%%
w = exp(1i*linspace(0,2*pi,1000));
Gr = @(s) Cr*((s*Er-Ar)\Br)+Dr;
for i = 1:1:size(w,2)
    Gr_eval(i) = Gr(w(i));
    G_eval(i) = H(w(i));
end

figure()
subplot(2,1,1)
semilogy(linspace(0,2*pi,1000),abs(G_eval),'k-'); hold on
semilogy(linspace(0,2*pi,1000),abs(Gr_eval),'r--'); 
legend('$f\circ \varphi$','$G_r\circ\varphi$', 'interpreter', 'latex')
hold off
subplot(2,1,2)
semilogy(linspace(0,2*pi,1000),abs(G_eval-Gr_eval)./abs(G_eval),'r-');
legend('rel. error');
hold off



% Discrete setting
time = 1e2;
x = zeros(n,time);
xr = zeros(r,time);
u = zeros(1,time); 
% u(7:time) = 1; % step input
u(10) = 1; % impulse input
% u(1:time) = sin(2*linspace(0,2*pi,time)); % sine input

for k = 5:1:time
    x(:,k) = (I-(12*h/25)*A)\((48/25)*x(:,k-1) - (36/25)*x(:,k-2) + (16/25)*x(:,k-3) - (3/25)*x(:,k-4) + B*(12*h/25)*u(k));
    xr(:,k) = (Ar\Er)*xr(:,k-1) - (Ar\Br)*u(k);
end
y = C*x;
yr = Cr*xr + Dr*u;

Gr = @(s) Cr*((s*Er-Ar)\(Br))+Dr;
funerror = @(z) abs((H(exp(1i*z))-Gr(exp(1i*z)))).^2;
H2D_bound = sqrt((1/(2*pi))*integral(funerror,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true));
    

figure()
subplot(1,2,1)
plot(1:1:time,y,'k'); hold on
plot(1:1:time,yr,'r--');
subplot(1,2,2)
semilogy(1:1:time,abs(y-yr),'r--'); hold on
semilogy(1:1:time,ones(1,time)*H2D_bound,'b:');

%% H2 norm
fomfom = @(z) abs(H(exp(1i*z))).^2;
sigma = [];
r_range = 4:2:14;
for i = 1:1:size(r_range,2)
    r = r_range(i)
    init = 0.1*randn(r,1)+0.1i*randn(r,1); init = (init./abs(init)).*rand(r,1);
    [Er,Ar,Br,Cr,Dr,sigma,~] = algorithm1(H,dH,r,init,10000,1e-6);
    Gr = @(s) Cr*((s*Er-Ar)\(Br))+Dr;
    funerror = @(z) abs((H(exp(1i*z))-Gr(exp(1i*z)))).^2;
    H2D_IRKA(i) = sqrt((1/(2*pi))*integral(funerror,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true));%./H2fomfom;
    
end

figure()
semilogy(r_range,H2D_IRKA,'r-x');
legend('$H^2(D)$ error TF-IRKA','Interpreter','latex','FontSize',14);
hold off


