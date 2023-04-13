clear 
clc
tic

T_st = 15.6;
p_st = 101.325;  % [kPa]
NG = WorkingFluid('Haynesville1.mix',T_st,p_st);
rho_st = NG.density();

d = 50; 
% dd = d*1e-4;
% d = d+dd;
Dx = 1; % useful in function
H = 30;
L = 106.5;
N = 1000;
x = linspace(0,1,N+1); % coordinates
Nf = 10;

% BD
h1 = 0;
h2 = 1;
h3 = 0; 
k1 = 1;
k2 = 0;
k3 = 0;

Ti = (300-32)/1.8; 
pi = 10200*6.895;  % [kPa]
pf = 1000*6.895;

%   initial alpha
NG = WorkingFluid('Haynesville1.mix',Ti,pi);
miu0 = NG.dynamicViscosity();
cg0 = 1e-3*NG.isothrmCompress(); % 1/Pa
rho0 = NG.density();
Zg0 = NG.compFactor();

% 
NG.pressure = pf;
miuf = NG.dynamicViscosity();
rhof = NG.density();
shale.k = 0.5e-6*(9.869233e-13);  % [D] to [m^2]
shale.e = 0.06;                                                                  
shale.Sg = 0.75;
shale.Sw = 1-shale.Sg;

% Langmuir model
pl = 710*6.895;
vl = 119*3.121e-2/9.072e2;
shale.ea0 = 0;
% shale.ea0 = ea_L(pi,pl,vl,cg0,rho0,rho_st);
ai = shale.k/(shale.e*shale.Sg+shale.ea0)/miu0/cg0;
tau = d^2/ai/3600/24;  % unit is day
M = Nf*4*L*H*d*(shale.e*shale.Sg+shale.ea0)*rho0;

% 
dt = 0.005;
% dt = floor(dt*1e5)/1e5; 
t = 0:dt:(2000*dt);
Nt = length(t);
i = 1:N+1;
P_old(i) = 1;
P(:,1) = P_old';
P(:,Nt) = zeros(N+1,1);

% 
e = ones(N+1,1);
At = diag(e)/dt; 
At(1,1) = 0;
At(N+1,N+1) = 0;

% 
for step = 2:Nt
    if ~iscolumn(P_old)
        P_old = P_old';   
    end
    P_rhs = P_old/dt;
    err = 1;
    k = 0;
    % pb = 6.895*4025.1*((step-1)*1.94)^(-0.114);
    P_rhs(1,1) = h3;
    P_rhs(N+1,1) = k3;

    while err > 1.0e-8
        k = k + 1;
        % 
        p = P_old.*(pi-pf)+pf;
        NG.pressure = p;
        niu = NG.kineticViscosity();
        miu = NG.dynamicViscosity();
        cg = 1e-3*NG.isothrmCompress(); % [1/Pa]
        rho = NG.density();
%         shale.ea = ea_L(p,pl,vl,cg,rho,rho_st);
        shale.ea = 0;
        a_den = (shale.e*shale.Sg+shale.ea);
        a_k = shale.k./a_den;
        a0 = a_k/miu0/cg0;
        a = niu.*a_k./(miu.*cg.*a0.*a0);
        D = a0./niu;

        % 
        [Ad,~] = diff_matrix_2(a,D,h1,h2,k1,k2,N,Dx);
        P_new = (Ad+At)\P_rhs;
        errs = abs(P_old - P_new);
        err = max(errs);
        P_old = P_new;
    end
    P(:,step) = P_old;
    if mod(step,10)==0
    disp(['The ' num2str(step) ' th time step uses iteration k = ' num2str(k)])
    figure(1)
    plot(x, P(:,step),'blue')
    hold on
    end
end
  delta_P = (P(2,:)-P(1,:))*(pi-pf);
%    delta_P = (P(2,:)-P(1,:));
%    p2 = P_old(2,1)*(pi-pf)+pf;  
    m_rate = Nf*4*H*L*N*shale.k*1000*delta_P*rhof/(d*miuf);  
    m_rate = m_rate';
 %  L2 = length(m_rate)-1;
  %  m_rate = (m_rate(2:end))+0.1*randn(L2,1);
  m_cum = cumsum(m_rate)*dt*tau*3600*24;
  RF = m_cum/M;
  figure(2)
  plot(t,RF);
  
% end
%     title('dimensionless P')
%     xlabel('dimensionless time')
%     ylabel('production m (kg/s)')
%     title('RF')
%     xlabel('dimensionless time')
%    m_rate2 = m_rate(2:end)/rho_st/2.8317e4*3600*24;  % [kg/s] to [MMscf/D]
    toc
    tE = toc/60;
    disp(tE)
