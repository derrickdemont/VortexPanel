
% Vortex Panel Potential Flow, Airfoil
% From Katz & Plotkin 
clear; clc; close all;


% Freestream
v = 1;
alfa =5;
Vinf = [v*cosd(alfa) v*sind(alfa)];


% Pt coords
af = readmatrix("naca0010.txt");
af = flip(af);          % FOR SOME REASON, THIS CODES REQUIRES THAT YOU START FROM LOWER TRAILING EDGE AND MARCH CW

dx = diff(af(:,1));
dy = diff(af(:,2));

% Linear interpolate between airfoil points
ctrlp = [0.5*diff(af(:,1))+af(1:end-1,1) 0.5*diff(af(:,2))+af(1:end-1,2)];
M = length(ctrlp);
N =M+1;

theta = atan2(dy,dx);                      % Local angle of panel


figure
plot(af(:,1),af(:,2),'o')
hold on
plot(ctrlp(:,1),ctrlp(:,2),'o')
axis equal


A = zeros(N+1);

for i = 1:M
    for j = 1:M
           xt = ctrlp(i,1) - af(j,1);
           zt = ctrlp(i,2) - af(j,2);
           x2t = af(j+1,1) - af(j,1);
           z2t = af(j+1,2) - af(j,2);

           x = xt*cos(theta(j)) +zt*sin(theta(j));
           z = -xt*sin(theta(j)) + zt*cos(theta(j));
           x2 = x2t*cos(theta(j)) +z2t*sin(theta(j));
           z2 = 0;

           if i == 1
               DL(j) = x2;
           end

           % Find R, thetas
           R1 = sqrt(x^2 + z^2);
           R2 = sqrt((x-x2)^2 + z^2);
           Th1 = atan2(z,x);
           Th2 = atan2(z,x-x2);

           % Compute vels in j ref frame
        if i == j
            U1L = 0.25;
            U2L = 0.25;
            W1L = -0.15916;
            W2L = 0.15916;
        else
            U1L = -(z*log(R2/R1)+x*(Th2-Th1)-x2*(Th2-Th1))/(6.28319*x2);
            U2L = (z*log(R2/R1)+x*(Th2-Th1))/(6.28319*x2);
            W1L = -((x2-z*(Th2-Th1))-x*log(R1/R2)+x2*log(R1/R2))/(6.28319*x2);
            W2L = ((x2-z*(Th2-Th1))-x*log(R1/R2))/(6.28319*x2);
        end

        % Transform local vels into global ref frame
        U1 = U1L*cos(-theta(j)) + W1L*sin(-theta(j));
        U2 = U2L*cos(-theta(j)) + W2L*sin(-theta(j));
        W1 = -U1L*sin(-theta(j)) + W1L*cos(-theta(j));
        W2 = -U2L*sin(-theta(j)) + W2L*cos(-theta(j));

        % Compute Influence matrix
        if j ==1
            A(i,1) = -U1*sin(theta(i)) + W1*cos(theta(i));
            holda = -U2*sin(theta(i)) + W2*cos(theta(i));
            B(i,1) = U1*cos(theta(i)) + W1*sin(theta(i));
            holdb = U2*cos(theta(i)) + W2*sin(theta(i));
        elseif j == M
            A(i,M) = -U1*sin(theta(i)) + W1*cos(theta(i)) + holda;
            A(i,N) = -U2*sin(theta(i)) + W2*cos(theta(i));
            B(i,M) = U1*cos(theta(i)) + W1*sin(theta(i)) + holdb;
            B(i,N) = U2*cos(theta(i)) + W2*sin(theta(i));
        else
            A(i,j) = -U1*sin(theta(i)) + W1*cos(theta(i)) + holda;
            holda = -U2*sin(theta(i)) + W2*cos(theta(i));
            B(i,j) = U1*cos(theta(i)) + W1*sin(theta(i)) + holdb;
            holdb = U2*cos(theta(i)) + W2*sin(theta(i));
        end
          
    end
    A(i,N+1) = cosd(alfa)*sin(theta(i)) - sind(alfa)*cos(theta(i));
end


A(N,1) = 1;
A(N,N) = 1;

% Solve for gamma
RHS = A(1:N,N+1);
newA = A(1:N,1:N);
G = newA\RHS;

% Convert vortex strengths to Vtangential
N = M+1;
Cl = 0;
for i=1:M
    vel = 0;
    for j=1:N
        vel = vel+ B(i,j)*G(j);
    end
    V(i) = vel + cosd(alfa)*cos(theta(i)) + sind(alfa)*sin(theta(i));
    Cl = Cl + V(i)*DL(i);
    Cp(i) = 1 - V(i)^2;
end


figure
plot(ctrlp(:,1),Cp)
xlabel('x/c'); ylabel('Cp');
set(gca, 'YDir','reverse');



