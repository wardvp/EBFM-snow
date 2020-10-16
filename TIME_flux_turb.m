%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TURBULENT SENSIBLE HEAT FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SHF,LHF] = TIME_flux_turb(C,Tsurf,IN,cond)
% Monin-Obukhov theory
a = 1;
b = 0.667;
c = 5;
d = 0.35;
kappa = 0.40;
z_u = 2.0;
z_tq = 2.0;
z0_u = 1d-3;
visc = 1.46d-5;

u_star_init = kappa.*IN.WS(cond)./log(z_u/z0_u);
Rey = u_star_init.*z0_u./visc;
b0(Rey<=0.135,1) = 1.25;
b0(Rey>0.135 & Rey<2.5,1) = 0.149;
b0(Rey>=2.5,1) = 0.317;
b1(Rey<=0.135,1) = 0.0;
b1(Rey>0.135 & Rey<2.5,1) = -0.550;
b1(Rey>=2.5,1) = -0.565;
b2(Rey<=0.135,1) = 0.0;
b2(Rey>0.135 & Rey<2.5,1) = 0.0;
b2(Rey>=2.5,1) = -0.183;
z0_tq_init = z0_u.*exp(b0+b1.*log(Rey)+b2.*log(Rey).^2);
t_star_init = kappa.*(IN.T(cond)-Tsurf)./(log(z_tq./z0_tq_init));
SHF_init = IN.Dair(cond).*C.Cp.*u_star_init.*t_star_init;

maxit=5;
for it=1:maxit
    if it == 1
        u_star = u_star_init;
        SHF = SHF_init;
    end
    SHF0 = SHF;
    
    L = IN.Dair(cond).*C.Cp.*u_star.^3.*IN.T(cond)./(kappa.*C.g.*SHF);
    cond_s = L>0;
    Psi_m = zeros(length(IN.Dair(cond)),1);
    Psi_h = zeros(length(IN.Dair(cond)),1);
    b0 = zeros(length(IN.Dair(cond)),1);
    b1 = zeros(length(IN.Dair(cond)),1);
    b2 = zeros(length(IN.Dair(cond)),1);
    
    Psi_m(cond_s) = -a.*z_u./L(cond_s) - b.*(z_u./L(cond_s)-c./d).*exp(-d.*z_u./L(cond_s)) - b.*c./d;
    Psi_h(cond_s) = -(1+2./3.*a.*z_tq./L(cond_s)).^1.5 - b.*(z_tq./L(cond_s)-c./d).*exp(-d.*z_tq./L(cond_s)) - b.*c./d + 1;
    
    y_m = (1-16.*z_u./L(~cond_s)).^0.25;
    y_h = (1-16.*z_tq./L(~cond_s)).^0.25;
    Psi_m(~cond_s) = log((0.5.*(1+y_m.^2))) + 2*log(0.5.*(1+y_m))-2.*atan(y_m) + pi./2;
    Psi_h(~cond_s) = 2.*log((0.5.*(1+y_h.^2)));
    
    Rey = u_star.*z0_u./visc;
    b0(Rey<=0.135,1) = 1.25;
    b0(Rey>0.135 & Rey<2.5,1) = 0.149;
    b0(Rey>=2.5,1) = 0.317;
    b1(Rey<=0.135,1) = 0.0;
    b1(Rey>0.135 & Rey<2.5,1) = -0.550;
    b1(Rey>=2.5,1) = -0.565;
    b2(Rey<=0.135,1) = 0.0;
    b2(Rey>0.135 & Rey<2.5,1) = 0.0;
    b2(Rey>=2.5,1) = -0.183;
    z0_tq = z0_u.*exp(b0+b1.*log(Rey)+b2.*log(Rey).^2);
        
    t_star = kappa.*(IN.T(cond)-Tsurf)./(log(z_tq./z0_tq)-Psi_h);
    u_star = kappa.*IN.WS(cond)./(log(z_u./z0_u)-Psi_m);
    
    SHF = IN.Dair(cond).*C.Cp.*u_star.*t_star;
    
    dSHF = abs(SHF-SHF0);
    SHF((dSHF>1 | u_star<0) & it==maxit) = SHF_init((dSHF>1 | u_star<0) & it==maxit);
    Psi_h((dSHF>1 | u_star<0) & it==maxit) = 0;
    Psi_m((dSHF>1 | u_star<0) & it==maxit) = 0;
    z0_tq((dSHF>1 | u_star<0) & it==maxit) = z0_tq_init((dSHF>1 | u_star<0) & it==maxit);
end

qsurf =  C.VP0.*C.eps./IN.Pres(cond).*exp(C.Ls/C.Rv/273.15*(1.0-273.15./Tsurf)) .* (Tsurf<273.15) ...
    + C.VP0.*C.eps./IN.Pres(cond).*exp(C.Lv/C.Rv/273.15*(1.0-273.15./Tsurf)) .* (Tsurf>=273.15);

q_star = kappa.*(IN.q(cond)-qsurf)./(log(z_tq./z0_tq)-Psi_h);
u_star = kappa.*IN.WS(cond)./(log(z_u./z0_u)-Psi_m);

LHF = IN.Dair(cond).*C.Ls.*u_star.*q_star .* (Tsurf<273.15) ...
    + IN.Dair(cond).*C.Lv.*u_star.*q_star .* (Tsurf>=273.15);

end