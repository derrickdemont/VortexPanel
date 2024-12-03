
% Source Potential Flow, Arbitrary shape
% From MCCORMICK p 50-54
clear; clc; close all;


% Freestream
Vinf = 1;
alfa =10;
V = [Vinf*cosd(alfa) Vinf*sind(alfa)];



% Pt coords

    % CYL
    I = 80;      % Number of panels - counts CCW
    T = 0:360/I:360;
    xpts = cosd(T);
    ypts = sind(T);



% Panel center coords
xp = diff(xpts)/2 + xpts(1:end-1);
yp = diff(ypts)/2 + ypts(1:end-1);
THETA = atan2(yp,xp);



V = repelem(V,I,1);

% S (length of panels)

S = (diff(xpts).^2 + diff(ypts).^2).^(0.5);

% Normal unit vectors

% or from book..
n = transpose([diff(ypts)./S; -diff(xpts)./S]);

% Radii
R = zeros(I^2,2);
k=0;
A = zeros(I,I);

% Create A matrix (based on all velocities summing to be tangential to
% surface)
for i=1:I
    for j=1:I
        k=k+1;
        R(k,:) = [xp(i)-xp(j) yp(i)-yp(j)];
        if i ~=j
            crk = R(k,:)/(2*pi*rssq(R(k,:))^2);
            A(i,j) = dot(crk,n(i,:));
        end
    end
end


A = A + diag(0.5*S.^(-1));      % Add Diagonals - influence of panel I on panel I 



% Set up Equation
B = dot(V,-n,2);

Q = A^(-1)*B;    % Solve for source strength

% Solve for velocity at each panel


V = [Vinf*cosd(alfa) Vinf*sind(alfa)];

for i=1:I
    VQII = (Q(i)/(2*S(i))).*n(i,:);

    VQIJ = zeros(I,2);
    for j=1:I

        R = [xp(i)-xp(j) yp(i)-yp(j)];
        if i~=j
            VQIJ(j,:) = Q(j).*R./(2*pi*rssq(R)^2);
        end
    end
    VR = V + sum(VQIJ,1) + VQII;
    VRmag(i) = rssq(VR);
    Cp(i) = 1- VRmag(i)^2;
end



%Panels for visualization

figure
plot(xpts,ypts)
xlabel('x'); ylabel('y')
hold on
scatter(xp,yp)

% Cp distribution
for i=1:I
    if THETA(i) < 0
        THETA(i) = THETA(i)+2*pi;
    end
end

figure
plot(THETA*180/pi, Cp)
xlabel('Angle [°]'); ylabel('Cp')

