% Parameters
params.rhoi = 917; % kg/m^3
params.rhow = 1000; % kg/m^3
params.H = 1000; % m
params.grav = 9.806; % acceleration due to gravity m/(s^2)
params.sinalpha = 10^(-3); % 1
params.psi0 = params.rhow*params.grav*params.sinalpha; % Pa
params.Tstar = 263.15; % Junction temperature (K)
params.Tmelt = 273.15; % melting temperature (K)
params.Tsurface = 273.15-26; % surface temperature (K)
params.DeltaT = params.Tmelt-params.Tsurface; % temperature difference (K)
params.A = 24*10^(-25); % 1/(s*Pa^n)
params.A0 = 3.5*10^(-25); % softness prefactor (1/(s*Pa^n))
params.Ainv = 4*10^(-24); % inversion softness (1/(s*Pa^n))
params.n = 3; % glen's law exponent
params.K = 2.1; % W/(m*K)
params.K0 = 9.83; % thermal conductivity prefactor W/(m*K)
params.cp = 2050; % J/(kg*K)
params.R = 8.314; % ideal gas constant
params.muw = 10^(-3); % Pa*s
params.Latent = 334000; % J/kg
params.accumulation = 0.1/(3600*24*365); % m/s
params.edotxy = [0.00136,0.0202]./(3600*365*24); % linear fit shear strain rate from Bindschadler shear margin : (1/2)*(du/dy)
params.edotxyscale = params.edotxy(2); % shear strain rate scale from Bindschadler shear margin : (1/2)*(du/dy)
params.etascale = (params.edotxyscale.^((1-params.n)/params.n))/(2*(params.A0.^(1/params.n))); % Pa*s
params.uvel = 10^(-7); % downstream velocity (m/s) and (x in km) % % (150 m/yr average from inversion) % %
params.k0 = 10^(-12); % permeability (m^2)
params.alpha = 7/3; % Permeability exponent (Hewitt uses 2, Schoof uses 7/3)
params.timescale = params.rhoi*params.cp*(params.H^2)/(params.K); % seconds
params.W = 2*(params.A0.^(-1/params.n))*(params.edotxyscale.^((params.n+1)/params.n)); % shear heating scale from Bindschadler shear margin

% nondimensional parameters
params.Pe = -(params.rhoi*params.cp*params.accumulation*params.H)/params.K; % Peclet number
params.epsilon = (params.rhoi*params.cp*params.DeltaT)/(params.rhow*params.Latent); % Shear heating
params.N0 = (params.etascale*params.K)./(params.rhoi*params.cp*params.H^2); % DIMENSIONAL effective pressure scale
params.qscale = params.K*params.epsilon/(params.rhoi*params.cp*params.H); % DIMENSIONAL m/s
params.Br = params.W*(params.H^2)/(params.K*params.DeltaT); % Shear heating (Brinkman)
params.kappa = params.k0*params.rhoi*params.cp*params.H*(params.epsilon.^(params.alpha-1))*(params.rhow-params.rhoi)*params.grav./(params.K*params.muw); % Scaled darcy flux
params.delta = params.N0./((params.rhow-params.rhoi)*params.grav*params.H); % Compaction pressure scale
params.g = -1; % gravity direction

% ice softness functions
fun.Qc = @(T) ((115000-60000)*(tanh(T-263.15)+1)./2)+60000;
fun.Gamma = @(T) fun.Qc(T)./(params.R*T);
fun.A = @(T,P,E) E*params.A0*exp(-fun.Gamma(T).*(1-(T./params.Tstar))).*(1. + 2.4*P);
fun.eta = @(T,P,E) ((fun.A(T,P,E)/params.A0).^(-1/params.n)).*...
    ((polyval(params.edotxy./params.edotxyscale,xgrid')).^((params.n+1)/params.n)); % nondimensional viscosity function
fun.shearheat = @(T,P,E) 2*(fun.A(T,P,E).^(-1/params.n)).*...
    ((polyval(params.edotxy,xgrid')).^((params.n+1)/params.n))*(params.H^2)./(params.DeltaT*params.K);