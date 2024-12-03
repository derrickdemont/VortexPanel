
% Source point Potential Flow, Arbitrary shape - for nonlifting bodies
% (like symmetric airfoil)
% From MCCORMICK p 50-54
clear; clc; close all;


% Freestream
Vinf = 1;
alfa = 0;
V = [Vinf*cosd(alfa) Vinf*sind(alfa)];



% Pt coords
en = 60;
[EX,EY] = naca4('0012',en);

    X = [];
    Y = [];
    N=2;
    for i=1:length(EX)-1
        ex = linspace(EX(i),EX(i+1),N);
        ey = linspace(EY(i),EY(i+1),N);
        X = [X ex(1:end-1)];
        Y = [Y ey(1:end-1)];
    end
    
    xpts = [X'; EX(end)];
    ypts = [Y'; EY(end)];
%     xpts = [xpts(1:en); xpts(en+2:end)];
%     ypts = [ypts(1:en); ypts(en+2:end)];


%     CORDS = readmatrix('naca0010.xls');
%     xpts = [CORDS(:,1); CORDS(1,1)];
%     ypts = [CORDS(:,2); CORDS(1,2)];
%     xpts = CORDS(:,1);
%     ypts = CORDS(:,2);


% Panel center coords
xp = diff(xpts)/2 + xpts(1:end-1);
yp = diff(ypts)/2 + ypts(1:end-1);

% X = [];
% Y = [];
% N=4;
% for i=1:length(xpts)-1
%     ex = linspace(xpts(i),xpts(i+1),N);
%     ey = linspace(ypts(i),ypts(i+1),N);
%     X = [X ex(1:end-1)];
%     Y = [Y ey(1:end-1)];
% end

% xp = X';
% yp = Y';

I = length(xp);

%Panels for visualization

figure
plot(xpts,ypts)
xlabel('x'); ylabel('y')
hold on
scatter(xp,yp)

V = repelem(V,I,1);

% S (length of panels)

S = (diff(xpts).^2 + diff(ypts).^2).^(0.5);


% Normal unit vectors
n = [diff(ypts)./S -diff(xpts)./S];

% S = repelem(S,N-1,1);
% n = repelem(n,N-1,1);

% Radii
R = zeros(I^2,2);
k=0;
A = zeros(I,I);

% Create A matrix (based on all velocities summing to be tangential to
% surface)
for i=1:I
    if i == 20
        HEY = 1;
    end
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






% Determine upper, lower surface
xu = [];
xl = [];
Cpu = [];
Cpl = [];

for i=1:I
    if yp(i) < 0
        xl = [xl; xp(i)];
        Cpl = [Cpl; Cp(i)];
    else
        xu = [xu; xp(i)];
        Cpu = [Cpu; Cp(i)];
    end
end

% Cp distribution

figure
scatter(xu, Cpu)
hold on
scatter(xl,Cpl)
xlabel('x/c'); ylabel('Cp')
legend('upper', 'lower')
set(gca, 'YDir','reverse')



