function [EffPout,hdim,Sdim,Qdim] = subglacialhydrology(xgrid,q,EffPend,Qup)
% input from temperate zone: xgrid, q, effective pressure at end

% parameters
ell = max(xgrid)*1000; % m
F = 650; f = 1/sqrt(F);
rhoi = 917; % kg/m^3
g = 9.806; % m/s^2
slope = 10^(-3);
Psi0 = rhoi*g*slope; % N/m^3
etai = 10^13; % Pa*s
Latent = 334000; % J/kg
G = 0.06;
etaw = 10^(-3);
k0 = 10^(-11)/3;
W = 10*1000;
R = 2*10^(-3);
ub = 10^(-7);
alpha = 4/3;
beta = 3/2;

Qtotal = Qup + q*W*ell;
% Nondimensional scales
Qscale = max(Qtotal);
hscale = ((Qscale*etaw)./(k0*W*Psi0)).^(1/3);
Sscale = (Qscale./(f*(Psi0.^(beta-1)))).^(1/alpha);
Nscale = (f*(Sscale^(alpha-1))*(Psi0^beta)*etai)./(rhoi*Latent);
delta = Nscale./(Psi0*ell);
curlyG = (etai./(hscale*Nscale))*((G./(rhoi*Latent))+R*ub);

startx = (1.01*min(xgrid)./max(xgrid));
fluxfun = @(x) interp1(xgrid./max(xgrid),Qtotal./Qscale,max(x,startx));

[x,y] = ode45(@channelfun,[1 0],EffPend./Nscale);
EffPout = interp1(x,y*Nscale,xgrid./max(xgrid)); % Pa

[h_out,~,S_out] = outputfuns(fluxfun(x),y);
hdim = interp1(x,h_out*hscale,xgrid./max(xgrid));
Sdim = interp1(x,S_out*Sscale,xgrid./max(xgrid));
Qdim = interp1(x,fluxfun(x)*Qscale,xgrid./max(xgrid));

    function dydx = channelfun(x,y)
        Q = fluxfun(x);
        N = y(1);
        
        h = curlyG./N;
        Qc = Q - (h.^3);
        if Qc>0
            S = ((Qc.^(beta./(beta-1)))./N).^((beta-1)./(alpha+beta-1));
            dNdx = (((Qc./(S.^(alpha))).^(1/(beta-1)))-1)./delta;
        else
            S = 0;
            dNdx = ((Q/(h.^3)) - 1)./delta;
        end
        dydx = dNdx;
    end

    function [h_out,Qc_out,S_out] = outputfuns(Q,N)
        for i = 1:length(Q)
            h(i) = curlyG./N(i);
            Qc(i) = Q(i) - (h(i).^3);
            if Qc(i)>0
                S(i) = ((Qc(i).^(beta./(beta-1)))./N(i)).^((beta-1)./(alpha+beta-1));
            else
                S(i) = 0;
            end
        end
        h_out = h;
        Qc_out = Qc;
        S_out = S;
    end
end