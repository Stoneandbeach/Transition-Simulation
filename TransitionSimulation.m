%Transition between 1s and 2p in a laser field
%Please allow some time for program to run

clear all
close all

%Define constants
e = 1.602e-19;
h = 6.626e-34;
hbar = h/(2*pi);
a0 = 5.292e-11;
Ryd = 1.09678e7; %Using Rh for hydrogen
c = 299792458;

%Input variables
%This controls the transition process
E0 = 1e6; %Laser field amplitude
Om = a0 / (2*sqrt(2)) * e * E0 / hbar; %Rabi frequency
d = 0; %1e10; %Detuning
Ga = 0; %1e10; %Decay, for example spontaneous emission
tstart = 0;
tend = 2.20e-10 / 2;
numpoints = 100;

%Frequencies for the states
w1 = 2*pi*c*Ryd; %1s
w2 = 2*pi*c*Ryd / 4; %2p

%Bloch equations
udot = @(t, R) -Ga/2 * R(1) + d * R(2);
vdot = @(t, R) -d * R(1) - Ga/2 * R(2) + Om * R(3);
wdot = @(t, R) -Om * R(2) - Ga * (R(3) + 1);

ODEFUN = @(t, R) [udot(t, R); vdot(t, R); wdot(t, R)];
TSPAN = linspace(tstart, tend, numpoints);
R0 = [0; 0; 1];

%Solve the Bloch equations
[T, R] = ode45(ODEFUN, TSPAN, R0);

u = R(:, 1);
v = R(:, 2);
w = R(:, 3);

%Coefficients for wavefunction
r12 = (u + 1i*v) / 2;
r21 = (u - 1i*v) / 2;
r11 = (w + 1) / 2;
r22 = (1 - w) / 2;

%Eigenfunctions for 1s and 2p
phi1 = @(x, y) sqrt(1/(4*pi)) * (1/a0)^(3/2) * 2 * exp(-sqrt(x^2 + y^2) / a0);
phi2 = @(x, y) sqrt(3/(4*pi)) * y / sqrt(x^2 + y^2) * (1/(2*a0))^(3/2) * ...
    2 * (1 - sqrt(x^2 + y^2)/a0) * exp(-sqrt(x^2 + y^2) / a0);

%Probability distribution
P =  @(tind, t, x, y) r11(tind) * abs(phi1(x, y))^2 + r22(tind) * abs(phi2(x, y))^2 + ... 
    (r12(tind) * exp(1i * (w2 - w1) * t) + r21(tind) * exp(1i * (w1 - w2) * t)) * phi1(x, y) * phi2(x, y);
    
x = linspace(-1 * a0, 1 * a0, 500);
[X, Y] = meshgrid(x);

%Store the plots of the probability distribution of electron, and the Bloch vector
fig = figure;
set(fig, 'Visible', 'off', 'position', [1 1 600 1200])
m = 0;
[h1, h2, h3] = sphere;

for t = 1 : length(T)
    %Transition
    for j = 1 : length(x)
        for k = 1 : length(x)
            Z(j, k) = P(t, T(t), X(j, k), Y(j, k));
        end
    end
    subplot(2,1,1)
    s = surf(Z);
    s.EdgeColor = 'none';
    colormap hot;
    view(2)
    drawnow
    
    %Bloch sphere
    subplot(2,1,2)
    quiver3(0, 0, 0, u(t), v(t), w(t))
    hold on
    mesh(h1, h2, h3, zeros(length(h1)), 'FaceAlpha', 0.1)
    if t > 2
        for m = 2 : t
            plot3([u(m) u(m - 1)], [v(m) v(m - 1)], [w(m) w(m - 1)], 'blue', 'LineWidth', 2)
        end
    end
    xline(0, '-', 'v');
    yline(0, '-', 'u');
    
    hold off
    axis([-1 1 -1 1 -1 1])
    view(25, 25)
    storeMovie(t) = getframe(fig);
end

%Display the plots
[height, width, p] = size(storeMovie(1).cdata);
fig2 = figure;
axis off
set(fig2, 'position', [50 50 width height])
movie(fig2, storeMovie)

%Once the script has run, the movie may be shown again by simply executing the following line:
%movie(fig2, storeMovie)