% Vortex Panel Potential Flow, Airfoil
% From Phillips and Katz & Plotkin 
% Inputs:       alfa : angle of attack
%               af   : airfoil geometry
% Outputs       Cp   : Pressure coefficient along x/c
%               Cl   : section lift coefficient
%               Cd   : section drag coefficient (pressure drag)
%               Cmle : section pitching moment coefficient about LE
%               Cmc4 : section pitching moment coefficient about c/4
clear; clc; close all;


% Freestream
v = 1;
alfa =0;
Vinf = [v*cosd(alfa) v*sind(alfa)];


% Pt coords
af = readmatrix("C:\Users\derri\Box\Derrick Wiberg\AeroTools\AirfoilGeom\naca2412.txt");
% af = flip(af);          % March CW about af surface starting at bottom TE
af = AFGeom(af,20,5);

dx = diff(af(:,1));
dy = diff(af(:,2));

% Linear interpolate between airfoil points
ctrlp = [0.5*dx+af(1:end-1,1) 0.5*dy+af(1:end-1,2)];
N = length(ctrlp);

a = atan2(dy,dx);                       % Local angle of panel   
s = vecnorm([dx dy],2,2);            % Panel lengths
   
n = [-sin(a) cos(a)];               % panel normal unit vector 

t = [cos(a) sin(a)];                % panel tangent unit vector


figure
plot(af(:,1),af(:,2),'o')
hold on
plot(ctrlp(:,1),ctrlp(:,2),'o')
axis equal

for i=1:N
      quiver(ctrlp(:,1),ctrlp(:,2),n(:,1),n(:,2),'Color','b')
      quiver(ctrlp(:,1),ctrlp(:,2),t(:,1),t(:,2),'Color','r')
end
hold off


A = zeros(N+1);
B = zeros(N+1);

for i = 1:N
    for j = 1:N

            % In panel coords
            x = ctrlp(i,1);
            y = ctrlp(i,2);

            
                x0 = af(j,1);
                y0 = af(j,2);

                r = [x-x0 y-y0];
                xp = dot(r,t(j,:));
                yp = dot(r,n(j,:));


            Phi = atan2(yp*s(j),yp^2+xp^2 - xp*s(j));
            Psi = 0.5*log((xp^2 + yp^2)/((xp-s(j))^2 + yp^2));

            if i == j
                upa = 0.25;
                upb = 0.25;
                wpa = -0.15916;
                wpb = 0.15916;
            else
                upa = (1/(2*pi*s(j)))*((s(j)-xp)*Phi+yp*Psi);
                upb = (1/(2*pi*s(j)))*(xp*Phi-yp*Psi);
                wpa = (1/(2*pi*s(j)))*(yp*Phi - (s(j) - xp)*Psi - s(j));
                wpb = (1/(2*pi*s(j)))*(-yp*Phi-xp*Psi + s(j));
            end

            trns2 = [cos(-a(j)) sin(-a(j)); -sin(-a(j)) cos(-a(j))]*[upa; wpa];
            ua = trns2(1);
            wa = trns2(2);

            trns3 = [cos(-a(j)) sin(-a(j)); -sin(-a(j)) cos(-a(j))]*[upb; wpb];
            ub(i,j) = trns3(1);
            wb(i,j) = trns3(2);



                if j ==1
                    aij = dot([ua wa],n(i,:));
                    bij = dot([ua wa],t(i,:));
                else
                    aij = dot([ua+ub(i,j-1) wa+wb(i,j-1)],n(i,:));
                    bij = dot([ua+ub(i,j-1) wa+wb(i,j-1)],t(i,:));
                end


            A(i,j) = aij;
            B(i,j) = bij;

    end
end

for i=1:N
    A(i,N+1) = dot([ub(i,N) wb(i,N)],n(i,:));
    B(i,N+1) = dot([ub(i,N) wb(i,N)],t(i,:));
end
A(N+1,1) = 1;
A(N+1,N+1) = 1;


Vinf = repmat(Vinf,N,1);
RHS = dot(-Vinf,n,2);

RHS(N+1) = 0;

gamma = A\RHS;

tmp = B*gamma;
V = tmp(1:N) + cosd(alfa)*cos(a) + sind(alfa)*sin(a);
% Cl = sum(V.*s);

Cp = 1 - V.^2;

figure
plot(ctrlp(:,1),Cp)
xlabel('x/c'); ylabel('Cp');
set(gca, 'YDir','reverse');

Cn = sum(-Cp.*dx);
Ca = sum(Cp.*dy);
Cl = Cn*cosd(alfa) - Ca*sind(alfa);
Cd = Cn*sind(alfa) + Ca*cosd(alfa);
Cmle = sum(Cp.*dx.*ctrlp(:,1) + Cp.*dy.*ctrlp(:,2));
Cmc4 = sum(Cp.*dx.*(ctrlp(:,1)-0.25) + Cp.*dy.*ctrlp(:,2));








        
