% Get n number of coordinates of NACA airfoil using NACA 4-digit funciton
% Inputs : fourdig : string of naca 4-digit airfoil (ex: '2412')
%                n : number of desired coordinates
% Outputs : X,Y : vectors of x,y coordinates

function [X,Y] = naca4(fourdig, n)

m = str2double(fourdig(1))*.01;
p = str2double(fourdig(2))*.1;
t = str2double(fourdig(3:4))*.01;

theta = linspace(0,pi,n);
x = (0.5*(1-cos(theta)));   % cosine clustering

if p ==0
    y =0.*x;   
else
    y = m/p^2 * (2*p.*x - x.^2);
end
yc = (m/(1-p)^2)*((1-2*p) + 2*p.*x - x.^2);
px = round(p*(n-1)+1);

yc(1:px) = y(1:px);

yt = (t/0.2)*(0.2969.*x.^0.5 - 0.126.*x - 0.3516.*x.^2 + 0.2843.*x.^3 - 0.1015.*x.^4);

theta = atan2(diff(yc), diff(x));
theta = [theta theta(end)];
xu = x - yt.*sin(theta);
xl = x + yt.*sin(theta);
yu = yc + yt.*cos(theta);
yl = yc - yt.*cos(theta);


X = [flip(x) x]';
Y = [flip(yu) yl]';

X = [X(1:n); X(n+2:end)];
Y = [Y(1:n); Y(n+2:end)];


end