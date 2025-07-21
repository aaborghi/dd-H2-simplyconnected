%% construct heat transfer function
clear all; close all; clc;
rng(42);

G = @(s) exp(-sqrt(s));
dG = @(s) -exp(-sqrt(s))./(2*sqrt(s));

% construct conformal map
phi = @(z) (1+z)./(1-z);
dphi = @(z) 2./(1-z).^2;
dphisqrt= @(z) (sqrt(2)./(1-z));
ddphi = @(z) 4./(1-z).^3;
phiinv = @(z) (z-1)./(z+1);


% algorithm 1
r = 10;
init = 0.1*randn(r,1)+0.1i*randn(r,1); init = (init./abs(init)).*rand(r,1); 
H = @(s) G(phi(s)).*dphisqrt(s);
dH = @(s) dG(phi(s)).*(dphi(s).*dphisqrt(s))...
         + G(phi(s)).*0.5.*(1./dphisqrt(s)).*ddphi(s);
[Er,Ar,Br,Cr,Dr,sigma,~] = algorithm1(H,dH,r,init,1000,1e-6);

% classic TFIRKA
init = 10+rand(r/2,1) + 10i*rand(r/2,1); init = [init;conj(init)];
[Er2,Ar2,Br2,Cr2,sigma2,~] = TFIRKA(G,dG,r,@(x) -conj(x),init,1000,1e-6);



figure()
plot_phase(H, [4.3,-4.3,4.3,-4.3]); hold on
p(1) = plot(real(eig(Ar,Er)),imag(eig(Ar,Er)),'ko','DisplayName','$\lambda_j$','markersize',10,'linewidth',2);
p(2) = plot(real(exp(1i*linspace(0,2*pi,100))),imag(exp(1i*linspace(0,2*pi,100))),'k-','linewidth',2);
p(3) = plot(real(sigma),imag(sigma),'k.','DisplayName','$\sigma$','markersize',15,'linewidth',2);
legend(p([1,3]),'Interpreter','latex','FontSize',20);
hold off
c = colorbar;


%% bode plots

% on imaginary axis
w = logspace(-2,3,1000);

Gr = @(s) (Cr*((phiinv(s)*Er-Ar)\(Br))+Dr).*((1-phiinv(s))./sqrt(2));
Gr2 = @(s) (Cr2*((s*Er2-Ar2)\(Br2)));
for i = 1:1:size(w,2)
    Gr_eval(i) = Gr(w(i)*1i);
    Gr2_eval(i) = Gr2(w(i)*1i);
    G_eval(i) = G(w(i)*1i);
end

figure()
subplot(2,1,1);
loglog(w,abs(G_eval),'k-'); hold on
loglog(w,abs(Gr_eval),'r--'); 
loglog(w,abs(Gr2_eval),'b:'); 
legend('G','algorithm 1', 'TF-IRKA','Interpreter','latex','FontSize',16);
subplot(2,1,2);
loglog(w,abs(G_eval-Gr_eval),'r-'); hold on
loglog(w,abs(G_eval-Gr2_eval),'b-'); 
legend('algorithm 1', 'TF-IRKA','Interpreter','latex','FontSize',16);
hold off

% on unit circle
w = exp(1i*linspace(0,2*pi,1000));
Gr2 = @(s) (Cr2*((phi(s)*Er2-Ar2)\(Br2))).*dphisqrt(s);
Gr = @(s) Cr*((s*Er-Ar)\(Br)); %+Dr
for i = 1:1:size(w,2)
    Gr_eval(i) = Gr(w(i));
    Gr2_eval(i) = Gr2(w(i));
    G_eval(i) = H(w(i));
end

figure()
subplot(2,1,1);
semilogy(linspace(0,2*pi,1000),abs(G_eval),'k-'); hold on
semilogy(linspace(0,2*pi,1000),abs(Gr_eval),'r--'); 
semilogy(linspace(0,2*pi,1000),abs(Gr2_eval),'b:'); hold off
legend('G','algorithm 1', 'TF-IRKA','Interpreter','latex','FontSize',16);
subplot(2,1,2);
semilogy(linspace(0,2*pi,1000),abs(G_eval-Gr_eval),'r-'); hold on
semilogy(linspace(0,2*pi,1000),abs(G_eval-Gr2_eval),'b-'); hold off
legend('algorithm 1', 'TF-IRKA','Interpreter','latex','FontSize',16);
hold off


%% E2D computation 
fomfom = @(z) abs(H(exp(1i*z)))';
H2fomfom = sqrt((1/(2*pi))*integral(fomfom,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true));

sigma = [];
sigma2 = [];

r_range = 4:2:18;
for i = 1:1:size(r_range,2)
    r = r_range(i)

    % algorithm1 E2D
    opts = size(sigma,1);
    init = 0.1*randn(r-opts,1)+0.1i*randn(r-opts,1); init = (init./abs(init)).*rand(r-opts,1); init = [sigma; init];
    [Er,Ar,Br,Cr,Dr,sigma,conv] = algorithm1(H,dH,r,init,1000,1e-6);
    %repeat until convergence
    while conv == 0
        init = 0.1*randn(r,1)+0.1i*randn(r,1); init = (init./abs(init)).*rand(r,1);
        [Er,Ar,Br,Cr,Dr,sigma,conv] = algorithm1(H,dH,r,init,1000,1e-6);
    end
    Gr = @(s) Cr*((s*Er-Ar)\(Br)) + Dr;
    funerror = @(z) abs((H(exp(1i*z))-Gr(exp(1i*z)))).^2;
    E2D_tfirka(i) = sqrt((1/(2*pi))*integral(funerror,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true))./H2fomfom;
    
    % classic TFIRKA
    opts2 = size(sigma2,1);
    init2 = 1 + rand(r-opts2,1) + 10i*rand(r-opts2,1); init2 = [sigma2; init2];
    [Er2,Ar2,Br2,Cr2,sigma2,conv2] = TFIRKA(G,dG,r,@(x) -conj(x),init2,1000,1e-6);
    %repeat until convergence
    if conv2 == 0
        init2 = 1 + rand(r/2,1) + 10i*rand(r/2,1); init2 = [init2;conj(init2)];
        [Er2,Ar2,Br2,Cr2,sigma2,conv2] = TFIRKA(G,dG,r,@(x) -conj(x),init2,1000,1e-6);
    end
    Gr = @(s) (Cr2*((phi(s)*Er2-Ar2)\(Br2))).*dphisqrt(s);
    funerror = @(z) abs((H(exp(1i*z))-Gr(exp(1i*z)))).^2;
    H2Cp_tfirka(i) = sqrt((1/(2*pi))*integral(funerror,0,2*pi,'RelTol',1e-8,'AbsTol',1e-12,'ArrayValued',true))./H2fomfom;
    
end

figure()
semilogy(r_range,E2D_tfirka,'r-x'); hold on
semilogy(r_range,H2Cp_tfirka,'b-o');
legend('$E^2(D)$ error Algorithm 1', '$E^2(D)$ error TFIRKA','Interpreter','latex','FontSize',14);
hold off


