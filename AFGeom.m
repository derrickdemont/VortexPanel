% Airfoil curve fitting 
% from Kulfan 'Universal Parametric Geometry Representation' 2008
% Inputs: af: airfoil coordinates (from airfoil database)
%         Np: number of nodes generated
%         NB: order of Bernstein poly (historically have used 5)
% Outputs: coords: generated interpolated airfoil coordinates (x,y)
%                  of size(Np*2 -1, 2) marching CW from bottom TE

function coords = AFGeom(af,Np,NB)
    N = length(af);
    
    psi = flip(af(1:ceil(N/2),1));
    ZETAu = flip(af(1:ceil(N/2),2));
    ZETAl = af(ceil(N/2):N,2);
    
    % Shape coeff - for round nose and pointed aft end airfoil
    N1 = 0.5;
    N2 = 1.0;
    
    % Class function
    C = psi.^(N1) .* (1-psi).^(N2);
    
    % Trailing edge thickness
    Delta_zeta_u = af(1,2);
    Delta_zeta_l = af(N,2);
    
    % Bernstein Poly
    i = linspace(0,NB,NB+1);
    K = factorial(NB)./(factorial(i).*factorial(NB-i));
    S = K.*psi.^i .*(1-psi).^(NB-i);
    BigS = C.*S;
    
    Bu = ZETAu - psi*Delta_zeta_u;
    Bl = ZETAl - psi*Delta_zeta_l;
    
    Au = lsqr(BigS,Bu, 1e-8);
    Al = lsqr(BigS,Bl, 1e-8);
    
    % Use Au, Al for larger number of points - cosine spacing for x
    psif = 0.5*(1- cos(linspace(0,pi,Np)))';
    Cf = psif.^(N1) .* (1-psif).^(N2);
    
    Sf = K.*psif.^i .*(1-psif).^(NB-i);
    BigSf = Cf.*Sf;
    
    zetau = BigSf*Au + psif*Delta_zeta_u;
    zetal = BigSf*Al + psif*Delta_zeta_l;
    
    coords = [[flip(psif); psif(2:Np)] [flip(zetal); zetau(2:Np)]];
