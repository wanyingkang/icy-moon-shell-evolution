%% based on thinice_spherical_dominant_coupled_linear_analysis5.m
% + correct wrong Lf
% linear analysis is only correct for n=1
addpath(genpath('~/Dropbox/matlab/code/'))
close all

%% control

T=100e6; % yr 
dt=1000;    % yr
Nt=T/dt;
showinstant=1e6/dt; %showinstant=1;
spheric=1;
linearana=1;
explore_moon=1;
explore_regime=0; 
skipping=1;
loadpre=0;
savefile=0;
if explore_regime==1
    linearana=1;
    loadpre=1;
    savefile=1;
end
if explore_moon==1
    explore_regime=0;
    skipping=1;
    loadpre=0;
    savefile=0;
end
warning('off')
mkani=0; casename='H22Ts100gamma12'; iani=1;
fileappendix='_wide_withpred';
moons={'enceladus','europa','ganymede','titan','callisto'};
%moons={'enceladus'};

%% loop over moons
for imoon=1:length(moons)
moon=char(moons(imoon));
fprintf([moon,'\n\n'])
%% physics constants
switch moon
    case 'ganymede'
    a=5268e3;
    g=1.428; g0=g;
    e=0.0013;
    Om=2*pi/(7.15*86400);
    Ts=110;
    rhob=1936;
    Hbar_rng=[150e3 475e3 800e3];
    % >150 https://www.nature.com/articles/417419a
    % ~300, ~108, https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JE004570

    case 'enceladus'
    a=2.52e5;
    g=0.113; g0=g;
    e=0.0047;
    Om=2*pi/(86400*1.37);           % orbital frequency
    Ts=80;   % surface temperature
    rhob=1610; % ambient density
    Hbar_rng=[14e3 21e3 28e3];  % ice depth (m)
    % Beuthe's setup for enceladus
    %rhoice=1000; rhowater=1000; g=0.1135; Hbar=23000; Ts=59;

    case 'europa'
    a=1.5608e6;                 % planetary radius, m
    g=1.315; g0=g;
    e=0.0094; 
    Om=2.048e-5;                % orbital frequency
    Ts=102;
    rhob=3013;
    Hbar_rng=[1e3 13e3 25e3];
    % https://www.sciencedirect.com/science/article/abs/pii/S0019103515000688
    % 19-25 https://www.nature.com/articles/417419a
    % 20-25 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2003JE002099
    % <15 https://www.sciencedirect.com/science/article/pii/S0019103507000656?casa_token=cVGC4vAO1UYAAAAA:NO_ecbqZ4C1q4vA36hb6oj9srCK67osXob65sp6r7hIMNRZNYz18YjJieZwOjdkMOpM21GD-dw
    
    case 'callisto'
    a=2410e3;                 % planetary radius, m
    g=1.235; g0=g;
    e=0.0074; 
    Om=2*pi/(86400*16.69);                % orbital frequency
    Ts=134;
    rhob=1834;
    Hbar_rng=[135e3 143e3 150e3];
    % https://link.springer.com/article/10.1007/s11208-005-0043-0
    
    case 'titan'
    a=2574e3;                 % planetary radius, m
    g=1.352; g0=g;
    e=0.0288; 
    Om=2*pi/(86400*15.95);                % orbital frequency
    Ts=94;
    rhob=1880;
    Hbar_rng=[50e3 125e3 200e3];
    % http://eps.berkeley.edu/~djheming/publications/nature12400.pdf
    
end

%% topography amplification factor
p_alpha=-1.5;
betaH2H0=1;
betaH0F=1;
upperbound=-p_alpha*betaH0F*(betaH2H0/2+1)-1;
lowerbound=betaH2H0*betaH0F/4;%*(Hbar-Hbar_conv)/Hbar;
betaH0Fmin=1/(-p_alpha*(betaH2H0/2+1)-betaH2H0/4);
upperbound_H0Fmin=-p_alpha*betaH0Fmin*(betaH2H0/2+1)-1;
lowerbound_H0Fmin=betaH2H0*betaH0Fmin/4;

%% physics parameter & grid
Nx=200;
x1=[-1:2/(Nx):1]'; 
x1=x1(1:end-1); 
dx1=2/(Nx);

lf=2; % tidal forcing degree
L=a*pi;                     % size of the domain, pi*Re
x=L.*x1; dx=L.*dx1;
theta=x./a;
k1=2*pi/L/2;
f0=2*2*pi/(86400*1.37);           % coriolis coef.
F20=(a*Om)^2*e*(3/2);
F22=(a*Om)^2*e*(sqrt(1^2+7^2)/8);
rhowater=1050;
rhoice=950;
yrtosec=365*86400;

Lf=330000;  % J/kg fusion energy
mue=3.5e9; % elastic shear modulus
nu0=0.33; % poisson ratio
n=1; % power in glen's law 
if n==1
    eta_melt=1e14; % viscosity at melting
    eta_max_Q=Inf; % almost does not make any difference
    eta_max_heat=eta_max_Q;
end
if n==3
    eta_melt=1/3.2e-24/1000000; % creep prefactor
end

Tm=273;  % ice melting temperature
Ea=59.4e3; % activation energy J/mol
Rg=8.31;   % gas constant J/K/mol
kappa0=651; % W/m, kappa=kappa0/T

Hbar=Hbar_rng(2);
Hbar_conv=8000;

%% spherical harmonics
r2d=180/pi;
N20=sqrt(4*pi/5);
N22=sqrt(96*pi/5);
N40=sqrt(4*pi/9);
N00=sqrt(4*pi/1);
P2=3/2*cos(k1.*x).^2-1/2;
P3=(2.5.*cos((pi/L).*x).^3-1.5.*cos((pi/L).*x));
P4=35/8*cos(k1.*x).^4-30/8*cos(k1.*x).^2+3/8;
P0=1;
if spheric
    wgt=abs(sin(k1.*x))/mean(abs(sin(k1.*x)));
else
    wgt=1;
end
dprimef=@(l) -(l-1).*(l+2); % eigen value of delta'=laplacian+2
ef=@(l1,m1,l2,m2,l) sqrt((2.*l1+1).*(2.*l2+1).*(2.*l+1)./4./pi).*Wigner3j([l1,l2,l],[m1,m2,-m1-m2]).*Wigner3j([l1,l2,l],[0,0,0]).*(-1).^(m1+m2);

%% central difference
d_dx=diag(ones(Nx-1,1)./dx./2,1)-diag(ones(Nx-1,1)./dx./2,-1)-diag(1/dx/2,Nx-1)+diag(1/dx/2,-(Nx-1));
d_dx2=(-2*diag(ones(1,Nx),0)+1*diag(ones(1,Nx-1),1)+1*diag(ones(1,Nx-1),-1)+1*diag(ones(1,1),Nx-1)+1*diag(ones(1,1),-(Nx-1)))/dx^2;
if spheric
cot1=cot(k1.*x);
cot1(1)=0; cot1(Nx/2+1)=0;
else
cot1=0;
end

%% regime to explore
if explore_regime || explore_moon
    if loadpre
        load(sprintf('symbrk_a%dpalpha%.1f%s.mat',a/1000,-p_alpha,fileappendix))
        N_Hbars=length(Hbars);
        N_Tss=length(Tss);
        tauQ=zeros(N_Hbars, N_Tss);
        tauH=tauQ;
        tauF=tauQ;
    else
        if explore_regime
        Hbars=Hbar_rng(1):(Hbar_rng(3)-Hbar_rng(1))/10:Hbar_rng(3); 
        Tss=(Ts-20:10:Ts+60);
        end
        if explore_moon
        Hbars=Hbar_rng; 
        Tss=Ts;
        end
        N_Hbars=length(Hbars);
        N_Tss=length(Tss);
        symbrk=zeros(N_Hbars, N_Tss);
        symbrk_range=zeros(N_Hbars, N_Tss,2);
        tauQ=zeros(N_Hbars, N_Tss);
        tauH=tauQ;
        tauF=tauQ;
    end
    
else
    Hbars=Hbar; N_Hbars=1;
    Tss=Ts; N_Tss=1;
end


%% loop over all Ts Hbar combinations
for iHbar=1:N_Hbars
for iTs=1:N_Tss
Ts=Tss(iTs);
Hbar=Hbars(iHbar);
disp(sprintf('Hbar=%d, Ts=%d...',Hbar,Ts))
%% steady state tidal heating profile (membrane only and mix/bend)
% visco-elastic constants
mu0=mue/(log(Tm/Ts))^(0+1)*integral(@(T) 1./(1 - (1i*mue/Om)./min(eta_melt.*exp((Ea/Rg).*(-1/Tm + 1./T)),eta_max_heat)).*(log(Tm./T) - log(Tm/Ts)).^0./T,Ts,Tm);
mu1=mue/(log(Tm/Ts))^(1+1)*integral(@(T) 1./(1 - (1i*mue/Om)./min(eta_melt.*exp((Ea/Rg).*(-1/Tm + 1./T)),eta_max_heat)).*(log(Tm./T) - log(Tm/Ts)).^1./T,Ts,Tm);
mu2=mue/(log(Tm/Ts))^(2+1)*integral(@(T) 1./(1 - (1i*mue/Om)./min(eta_melt.*exp((Ea/Rg).*(-1/Tm + 1./T)),eta_max_heat)).*(log(Tm./T) - log(Tm/Ts)).^2./T,Ts,Tm);
muinv=mu2-mu1^2/mu0;
alphainv=1/(2*(1+nu0)*mu0*Hbar);
Dinv=2*muinv/(1-nu0)*Hbar^3;
alpha0=alphainv;
D0=Dinv;
chii0=-imag(mu1/mu0)*Hbar/a;
% solve thin ice model, (F,w)
L0=[a^3*dprimef(lf),D0*dprimef(lf)*(dprimef(lf)-1+nu0)+1*a^4*rhowater*g*(1-3/(2*lf+1)*rhowater/rhob);...
    alpha0*dprimef(lf)*(dprimef(lf)-1-nu0),-1/a*dprimef(lf)];
R0=[a^4*rhowater;0];
sol=L0\R0;
F0=sol(1); w0=sol(2);
% calculate heat generation
UAmps=[F20*N20;F22*N22];
ms=[0,2];
efs=[ef(lf,ms(1),lf,-ms(1),0),ef(lf,ms(1),lf,-ms(1),2),ef(lf,ms(1),lf,-ms(1),4);...
    ef(lf,ms(2),lf,-ms(2),0),ef(lf,ms(2),lf,-ms(2),2),ef(lf,ms(2),lf,-ms(2),4)];
dprimes=[dprimef(0),dprimef(2),dprimef(4)];

Hmem0=sum(real((UAmps.^2).*(-Om/2*imag(alpha0)).*(dprimef(lf)*dprimef(lf)*F0*conj(F0).*efs...
    +(1+nu0)/4*(dprimes.^2+2*(1-dprimef(lf)-dprimef(lf))*dprimes+((dprimef(lf)-dprimef(lf))^2+2*(dprimef(lf)+dprimef(lf))-8)).*efs*F0*conj(F0))));
Hbend0=sum(real((UAmps.^2).*(Om/2*imag(D0)/a^4).*(dprimef(lf)*dprimef(lf)*w0*conj(w0).*efs...
    +(1-nu0)/4*(dprimes.^2+2*(1-dprimef(lf)-dprimef(lf))*dprimes+((dprimef(lf)-dprimef(lf))^2+2*(dprimef(lf)+dprimef(lf))-8)).*efs*w0*conj(w0))));
Hmix0=sum(real((UAmps.^2).*(-Om/2*chii0/a/4).*((dprimes.^2+2*(1-dprimef(lf)-dprimef(lf))*dprimes+((dprimef(lf)-dprimef(lf))^2+2*(dprimef(lf)+dprimef(lf))-8))).*2.*real(efs*w0*conj(F0))));

tidal_mem0=Hmem0(1)*(P0/N00)+Hmem0(2)*(P2/N20)+Hmem0(3)*(P4/N40);
tidal_mixbend0=(Hmix0(1)+Hbend0(1))*(P0/N00)+(Hmix0(2)+Hbend0(2))*(P2/N20)+(Hmix0(3)+Hbend0(3))*(P4/N40);
mixbendfactor0=1;
tauH(iHbar,iTs)=(Lf*rhoice*Hbar)/(Hmem0(1)*(P0/N00)+(Hmix0(1)+Hbend0(1))*(P0/N00));

%% Beuthe analytical
Phi_0=3/5*(3+4)*e^2;
Phi_2=-3/7*(3+8)*e^2*P2;
Phi_4=9/140*(27+16)*e^2*P4;
Phi_A=Phi_0+Phi_2+Phi_4;
Phi_B=Phi_0+Phi_2/2-(2/3)*Phi_4;
Phi_C=Phi_0-Phi_2+(1/6)*Phi_4;
h2=0.0448-0.0015i;
f2=(a*abs(h2)/g*Om^2)^2;
Hmem_beuthe=f2*Om/(5+nu0)^2*imag(1/alpha0)*(4*(1-nu0)*Phi_A+6*(1+nu0)*Phi_C);
Hmix_beuthe=-f2*2*Om/(5+nu0)*chii0*real(1/alpha0)*(4*Phi_A-6*Phi_C);
Hbend_beuthe=f2*Om*imag(D0)/a^2*(4*(1+nu0)*Phi_A+6*(1-nu0)*Phi_C);
showcomp=0;

%% plot heating
if explore_regime==0
figure
plot(-theta*r2d,tidal_mem0,'k','LineWidth',3)
hold on
plot(-theta*r2d,tidal_mixbend0,'k--','LineWidth',3)
plot(-theta*r2d,theta.*0,'k:')
xlim([0,180])
xticks([0:30:180])
%xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('\theta')
title('$\mathbf{\mathcal{H}_0}  \ (W/m^2)$','Interpreter','latex')
set(gca,'FontSize',20)
legend({'Membrane','Mix+Bend'})
if showcomp
    plot(-theta, (Hmem_beuthe),'ro','DisplayName','Membrane from Beuthe19');
    plot(-theta, (Hmix_beuthe+Hbend_beuthe),'rx','DisplayName','Mix+bend from Beuthe19');
    plot(-theta, (Hmem_beuthe)*1.25,'bo','DisplayName','1.25 x Membrane from Beuthe19');
    plot(-theta, (Hmix_beuthe+Hbend_beuthe)*1.25,'bx','DisplayName','1.25 x Mix+bend from Beuthe19');
    ylim([-0.5e-3,4.5e-3])
    % we drop the self-gravity term and thus slightly overestimate H.
pause
end
end

%% ice flow parameter
%eta_mean=integral(@(T) eta_melt.*exp((Ea./Rg).*(1./T-1./Tm))./T,Ts,Tm)./log(Tm/Ts);
ice_flow0=2*((rhowater-rhoice)*rhoice/rhowater*g)^n*Hbar^(2*n+1)/((log(Tm/Ts))^(n+2))*integral2(@(Tp,T) 1./min(eta_melt.*exp(-Ea/Rg/Tm+(Ea/Rg)./Tp),eta_max_Q).*((log(Tm./Tp)).^n)./Tp./T.*(Tp<=T),Ts,Tm,Ts,Tm);
tauQ(iHbar,iTs)=a^2/ice_flow0;

%% ice conductivity 
F1=kappa0*log(Tm/Ts); % F=F1/Hbar
tauF(iHbar,iTs)=(Lf*rhoice*Hbar)/(F1/Hbar);

%% amplify factor for tidal heating
if explore_regime==1
symbrk_range(iHbar,iTs,2)=mean(wgt.*F1./Hbar)/mean(wgt.*(tidal_mem0+mixbendfactor0.*tidal_mixbend0));
gamma0=0.7*symbrk_range(iHbar,iTs,2);
maxiter=100; % max interation
else
gamma0=38; 
% 36.9-38.1 for Hbar=24000; Ts=80; palpha=-1.5; eta_melt=1e14;
maxiter=1;
end

if ~skipping
for igamma0=1:maxiter 
disp(sprintf('gamma=%3.2f......',gamma0));
%% start loop
hbarp=zeros(Nx,Nt);
if ~linearana
hbarp(1:Nx/2+1,1)=0.05.*(rand(Nx/2+1,1)-0.5);
hbarp(Nx/2+2:end,1)=hbarp(Nx/2:-1:2,1);
%hbarp(:,1)=-0.3*cos(2.*theta)+0.1.*cos(3.*theta);
%hbarp(:,1)=hbarp(:,1)-0.3*cos(2.*theta);
end

for it=2:Nt
    %% tendency
    heat_loss=F1./Hbar./(1+hbarp(:,it-1))./(rhoice*Lf*Hbar);
    heat_loss_mean=mean(wgt.*heat_loss);
    
    tidal_melting0=gamma0*(tidal_mem0.*(1+hbarp(:,it-1)).^p_alpha+mixbendfactor0.*tidal_mixbend0)./(rhoice*Lf*Hbar);
    tidal_melting0_mean=mean(wgt.*tidal_melting0);
    
    heat_factor=heat_loss_mean/tidal_melting0_mean;
    if it==Nt
    fprintf('F/H=%.3f\n',heat_factor)
    end
    %tidal_melting=tidal_melting0.*(heat_factor);
    tidal_melting=tidal_melting0.*(1);
    
    ice_flow_div=ice_flow0.*((n+2).*((1+hbarp(:,it-1)).^(n+1)).*(d_dx*hbarp(:,it-1)).^(n+1) ...
        + n.*((1+hbarp(:,it-1)).^(n+2)).*(d_dx*hbarp(:,it-1)).^(n-1).*(d_dx2*hbarp(:,it-1))...
        + cot1./a.*(1+hbarp(:,it-1)).^(n+2).*(d_dx*hbarp(:,it-1)).^n);
   
    %hbarp(:,it)=hbarp(:,it-1)+(-tidal_melting+heat_loss+ice_flow_div).*(dt*yrtosec);
    hbarp(:,it)=hbarp(:,it-1)+(-tidal_melting+tidal_melting0_mean+heat_loss-heat_loss_mean+ice_flow_div).*(dt*yrtosec);
    
    %% convection
     hbarp_conv=-(Hbar-Hbar_conv)/Hbar;
     if any(hbarp(:,it)<hbarp_conv)
     hbarp(:,it)=max(hbarp(:,it),hbarp_conv);
     if linearana
         break;
     end
     end
     
    %% remove mean
     hbarp(:,it)=hbarp(:,it)-mean(hbarp(:,it).*wgt);
    
    %% plot
    if explore_regime==0 && showinstant~=0 && mod(it,showinstant)==0
    hf=figure(1);
    clf
    yyaxis left
    %plot(-theta*r2d,hbarp(:,it),'k-','LineWidth',2)
    % ylabel('h''')
    plot(-theta*r2d,(1+hbarp(:,it))*Hbar/1e3,'k-','LineWidth',2)
    ylabel('H (km)')
    set(gca,'ycolor','k')
    yyaxis right
    %plot(-theta*r2d,(tidal_melting-tidal_melting0_mean).*(yrtosec),'r--','LineWidth',2)
    %hold on
    %plot(-theta*r2d,ice_flow_div.*(yrtosec),'b--','LineWidth',2)
    %plot(-theta*r2d,(heat_loss-heat_loss_mean).*(yrtosec),'g--','LineWidth',2)
    %ylabel('dh''/dt (yr^{-1})')
    plot(-theta*r2d,(tidal_melting-tidal_melting0_mean).*(yrtosec*1e6)*(Hbar/1e3),'r--','LineWidth',2)
    hold on
    plot(-theta*r2d,ice_flow_div.*(yrtosec*1e6)*(Hbar/1e3),'b--','LineWidth',2)
    plot(-theta*r2d,(heat_loss-heat_loss_mean).*(yrtosec*1e6)*(Hbar/1e3),'g--','LineWidth',2)
    ylabel('dH/dt (km/Myr)')
    set(gca,'ycolor','k')
    set(gca,'FontSize',20)
    xlabel('\theta')
    xlim([0,180])
    xticks([0:30:180])
    title([sprintf('%3.1f',dt*(it-1)/1e6),' Myr'])
    if mkani
        if iani==1
            mkdir(casename)
        end
        saveas(hf,sprintf([casename,'/tmp%03d.tiff'],iani),'tiff')
        iani=iani+1;
    end
    %pause
    end
end

if linearana
%% linear analysis
    hp=hbarp(:,it);
    hp(Nx/4+1:Nx/2+1)=hp(Nx/4+1:-1:1); hp(Nx/2+1:Nx)=hp(1:Nx/2);
    if ~all(isfinite(hp)) || max(diff(hp))>0.1 || max(abs(hp(3:end)-hp(1:end-2)))>0.1 || any(hp<hbarp_conv+0.02)
        symbrk_range(iHbar,iTs,2)=gamma0;
        if symbrk_range(iHbar,iTs,2)-symbrk_range(iHbar,iTs,1)<0.05
            disp(sprintf('Hbar=%d, Ts=%d: fail!!!',Hbar,Ts))
            break;
        else
            gamma0=mean(symbrk_range(iHbar,iTs,:));
            continue;
        end
    end
    M=zeros(Nx,Nx);
    digits(10)
    hpcub=(1+vpa(hp)).^3;
    % hpn2=(1+hp).^(n+2);
    % hpn1=(1+hp).^(n+1);
    dhp=d_dx*vpa(hp);
    d2hp=d_dx2*vpa(hp);


    % M= diag(-gamma0.*tidal_mem0./(Lf*rhoice*Hbar).*p_alpha.*(1+hp).^(p_alpha-1))...
    %     +diag((-F1/Hbar/(Lf*rhoice*Hbar)).*(1+hp).^(-2))...
    %     +(ice_flow0*(n+2)).*diag((d_dx*(hpn1)).*dhp.^n+hpn1.*n.*dhp.^(n-1).*d2hp)...
    %     +(ice_flow0*(n+2)).*diag(hpn1.*dhp.^n)*d_dx...
    %     +(ice_flow0*n).*diag((d_dx*(hpn2)).*dhp.^(n-1)+hpn2.*(n-1).*dhp.^max(n-2,0).*d2hp).*d_dx...
    %     +(ice_flow0*n).*diag(hpn2.*dhp.^(n-1)).*d_dx2...
    %     +(ice_flow0*(n+2)).*diag(hpn1.*dhp.^n.*(cot1./a))...
    %     +(ice_flow0*n).*diag(hpn2.*dhp.^(n-1).*(cot1./a)).*d_dx;

    % M= diag(-gamma0.*tidal_mem0./(Lf*rhoice*Hbar).*p_alpha.*(1+hp).^(p_alpha-1))...
    %     +diag((-F1/Hbar/(Lf*rhoice*Hbar)).*(1+hp).^(-2))...
    %     +(2*ice_flow0).*diag((d_dx*hpcub))*d_dx...
    %     +ice_flow0*diag(hpcub)*d_dx2...
    %     +ice_flow0*diag(d_dx2*hpcub)...
    %     +ice_flow0*diag(cot1./a.*(d_dx*hpcub))...
    %     +ice_flow0*diag(cot1./a.*hpcub)*d_dx;

    M= diag(-gamma0.*tidal_mem0./(Lf*rhoice*Hbar).*p_alpha.*(1+vpa(hp)).^(p_alpha-1))...
        +diag((-F1/Hbar/(Lf*rhoice*Hbar)).*(1+vpa(hp)).^(-2))...
        +(2*ice_flow0).*diag((d_dx*vpa(hpcub)))*d_dx...
        +ice_flow0*diag(vpa(hpcub))*d_dx2...
        +ice_flow0*diag(d_dx2*vpa(hpcub))...
        +ice_flow0*diag(cot1./a)*d_dx*diag(vpa(hpcub));

    M=double(M-((wgt./Nx)'*M));
    
    
    [eigvec,grow]=eigs(M,5,'largestreal');
    grow=diag(grow);
    % exclude the nonphysical mode (east and west hemispheres should be symmetric by assumption)
    for ieig=1:5
    if eigvec(Nx/2-2,ieig)*eigvec(Nx/2+4,ieig)>0 
        break;
    end
    end
end

%% next guess 
if explore_regime==1
    if grow(ieig)>1e-16 && grow(ieig)<1e-12 && eigvec(1,ieig)*eigvec(Nx/2+1,ieig)<0
        symbrk(iHbar,iTs)=1;
        disp(sprintf('Hbar=%d, Ts=%d: successful!!! gamma=%3.5f',Hbar,Ts,gamma0));
        break;
    elseif grow(ieig)<1e-16
        symbrk_range(iHbar,iTs,1)=gamma0;
        if igamma0==1 % if first iteration, refine the range so that lateral variation of thickness is considered in tidal heating calculation
            symbrk_range(iHbar,iTs,2)=gamma0*heat_factor;
        end
    elseif grow(ieig)>1e-12
        symbrk_range(iHbar,iTs,2)=gamma0;
    end
    if symbrk_range(iHbar,iTs,2)-symbrk_range(iHbar,iTs,1)<0.05
        disp(sprintf('Hbar=%d, Ts=%d: fail!!!',Hbar,Ts))
        break;
    else
        gamma0=mean(symbrk_range(iHbar,iTs,:));
    end
end
%pause
end
if linearana    
    %% plot unstable mode
    % construct an even function from eigvec(:,2) and eigvec(:,3)
    figure
    subplot(1,2,1)
    yyaxis left
    plot(-theta*r2d,(hp+1)*Hbar/1e3,'k-','LineWidth',2)
    ylabel('H_{eq} (km)')
    set(gca,'ycolor','k')
    %ylim([-0.6 0.2])
    %yticks([-0.6:0.2:0.4])
    ylim([0, Hbar/1e3*1.2])
    yyaxis right
    plot(-theta*r2d,tidal_melting.*(yrtosec*1e6).*(Hbar/1e3),'r--','LineWidth',2)
    hold on
    plot(-theta*r2d,ice_flow_div.*(yrtosec*1e6).*(Hbar/1e3),'b--','LineWidth',2)
    plot(-theta*r2d,heat_loss.*(yrtosec*1e6).*(Hbar/1e3),'g--','LineWidth',2)
    set(gca,'ycolor','k')
    set(gca,'FontSize',20)
    xlabel('\theta')
    xlim([0 180])
    xticks([0:30:180])
    %xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
    %ylabel('dh''/dt (yr^{-1})')
    ylabel('dH/dt (km/Myr)')
    hold off
    title('equilibrium state')
    
    subplot(1,2,2)
    plot(-theta*r2d,eigvec(:,ieig)./max(abs(eigvec(:,ieig))),'k-','LineWidth',2)
    hold on
    text(35.,0.9,['growth rate=',sprintf('%3.3g',grow(ieig)*yrtosec*1e6),'/Myr'],'FontSize',18)

    xlabel('\theta')
    ylabel('growing structure')
    set(gca,'Fontsize',20)
    xlim([0 180]) 
    xticks([0:30:180])
    %xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
    set(gca,'YAxisLocation','right')
    title('most unstable mode')
    hold off
    set(gcf,'Position',[1470         735         883         333])
    drawnow
else
    if mkani
        path1 = getenv('PATH');
        path1 = [path1 ':/usr/local/Cellar'];
        path1 = [path1 ':/usr/local/bin'];
        setenv('PATH', path1);
        cd(casename);
        system(['ffmpeg -y -r 1 -i "tmp%03d.tiff" -c:v libx264 -r 20 -vf scale=1200:900 -pix_fmt yuv420p ' casename '.mp4']);
        system(['rm -f *.tiff'])
        cd ..;
    end
    %%
    figure(2)
    clf
    %yyaxis left
    plot(-theta*r2d,(1+hbarp(:,end))*Hbar/1e3,'k-','LineWidth',3)
    ylabel('H (km)')
    set(gca,'ycolor','k')
    %ylim((1+[-1 0.4])*Hbar/1e3)
    ylim([0 max((1+hbarp(:,end))*Hbar/1e3)*1.1])
    %yticks([-1:0.2:0.4])
    %yyaxis right
    hold on
    plot(-theta*r2d,(1+hbarp(:,1))*Hbar/1e3,'k--','LineWidth',1)
    %ylabel('h(initial) (km)')
    %set(gca,'ycolor','k')
    set(gca,'FontSize',20)
    legend({'final state','initial state'},'Location','southwest','Fontsize',20)
    xlabel('\theta')
    %ylim((1+[-1 0.4])*Hbar/1e3)
    %yticks([-1:0.2:0.4]/10)
    xlim([0, 180])
    xticks([0:30:180])
    %xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
    title([sprintf('%3.1f',T/1e6),' Myr'])
    set(gcf,'Position',[1000         500         600         400])
end % linearana
end % skipping

end
end
%% save file
if savefile
save(sprintf('symbrk_a%dpalpha%.1f%s.mat',a/1000,-p_alpha,fileappendix),'symbrk','symbrk_range','Tss','Hbars','tauF','tauH','tauQ')
end

if explore_regime   
    %% time scale ratios
    tauF_tauQ=tauF./tauQ;
    tauF_tauQ_symbrk=tauF_tauQ;
    tauF_tauQ_symbrk(tauF_tauQ>upperbound)=NaN;
    tauF_tauQ_symbrk(tauF_tauQ<lowerbound)=NaN;
    
    %%
    figure
    imagesc(Tss,Hbars/1000, symbrk)
    box on
    set(gca,'Ydir','normal')
    colormap(gray)
    caxis([0 1])
    hold on
    %contour(Tss,Hbars/1000, tauF_tauQ,[0.5,0.5],'LineColor','r','LineWidth',6,'LineStyle','--')
    contour(Tss,Hbars/1000, tauF_tauQ,[lowerbound upperbound],'LineColor','r','LineWidth',6,'LineStyle','--')
    xlabel('T_s (K)')
    ylabel('H_0 (km)')
    set(gca,'FontSize',20)
    title(['p_\alpha=',sprintf('%.1f',p_alpha)])
    ylim([15 36])
    yticks([15:3:36])
    text(1.3*min(Tss),1.25*min(Hbars)/1000,'Symmetry Breaking Allowed', 'FontSize',20)
%    %%
%     figure
%     contourf(Tss,Hbars/1000,tauF_tauQ_symbrk,30,'LineColor','none')
%     colormap(jet)
%     caxis([lowerbound upperbound])
%     hold on
%     contour(Tss,Hbars/1000, symbrk,1,'Color','k','LineWidth',3)
%     xlabel('T_s (K)')
%     ylabel('H_0 (km)')
%     set(gca,'FontSize',20)
%     title(['p_\alpha=',sprintf('%.1f',p_alpha)])
%     text(1.3*min(Tss),1.25*min(Hbars)/1000,'Symmetry Breaking Allowed', 'FontSize',20)

end % explore regime

if explore_moon
    figure(15)
    hold on
    maxgamma=50;
    tauFH_mean=min(tauF(2)./tauH(2),1);
    tauFH_max=min(maxgamma*max(tauF./tauH),1);
    tauFH_min=min(min(tauF./tauH),1);
    tauFQ_mean=tauF(2)./tauQ(2);
    tauFQ_max=max(tauF./tauQ);
    tauFQ_min=min(tauF./tauQ);
    errorbar(tauFH_mean, tauFQ_mean, tauFQ_mean-tauFQ_min, tauFQ_max-tauFQ_mean...
        , tauFH_mean-tauFH_min, tauFH_max-tauFH_mean,'Marker', 'o','Color','k','LineWidth',2);
    text(tauFH_mean*1.2, tauFQ_mean*2, moon,'FontSize',20)
    xlabel('$\tau_{\mathcal{F}}/\tau_{\mathcal{T}}$','Interpreter','latex')
    ylabel('$\tau_{\mathcal{F}}/\tau_{\mathcal{Q}}$','Interpreter','latex')
    set(gca,'FontSize',20)
    set(gca,'yscale','log','xscale','log')
    box on
    if imoon==length(moons)
        xs=[betaH0Fmin,1,1,betaH0Fmin];
        ys=[lowerbound_H0Fmin,lowerbound,upperbound,upperbound_H0Fmin];
        patch(xs,ys,'red','FaceAlpha',.3,'EdgeColor','none')
    end
    %title(['p_\alpha=',sprintf('%.1f',p_alpha)])
end

end % moon