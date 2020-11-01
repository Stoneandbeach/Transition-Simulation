%Simulation of the transition between 1s and 2p of a hydrogen atom in a
%laser field, displaying probability distribution of electron and evolution
%of the Bloch vector.
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
tend = 2.20e-10;
numpoints = 200;

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

linepoints = 500;
x = linspace(-1 * a0, 1 * a0, linepoints);
[X, Y] = meshgrid(x);
Z = zeros(linepoints, linepoints);

for j = 1 : linepoints
    for k = 1 : linepoints
        Phis1(j, k) = phi1(X(j, k), Y(j, k));
        Phis2(j, k) = phi2(X(j, k), Y(j, k));
    end
end
Phi1abs = abs(Phis1).^2;
Phi2abs = abs(Phis2).^2;

Zs = zeros(linepoints, linepoints, numpoints);
for t = 1 : length(T)
    Zs(:, :, t) = r11(t) * Phi1abs + r22(t) * Phi2abs + ...
        (r12(t) * exp(1i * (w2 - w1) * T(t)) + r21(t) * exp(1i * (w1 - w2) * T(t))) .* Phis1 .* Phis2;
end

fig = figure;
set(fig, 'Visible', 'on', 'position', [1 1 500 1000])

pl1 = subplot(2,1,1);
s = surf(pl1, Z);
s.EdgeColor = 'none';
colormap hot;
view(2)

pl2 = subplot(2,1,2);
axis([-1 1 -1 1 -1 1])
view(25, 25)
hold(pl2)
q = quiver3(pl2, 0, 0, 0, u(1), v(1), w(1));
pl3 = plot3(pl2, [u(1) u(1)], [v(1) v(1)], [w(1) w(1)], 'blue', 'LineWidth', 2);
[h1, h2, h3] = sphere;
mesh(pl2, h1, h2, h3, zeros(length(h1)), 'FaceAlpha', 0.1)
xline(0, '-', 'v');
yline(0, '-', 'u');
hold(pl2)

%Store the plots of the probability distribution of electron, and the Bloch vector
for t = 1 : length(T)
    s.ZData = Zs(:, :, t);
    
    %Bloch sphere
    q.UData = u(t);
    q.VData = v(t);
    q.WData = w(t);
    if t > 2
        hold on
        plot3(pl2, [u(t) u(t - 1)], [v(t) v(t - 1)], [w(t) w(t - 1)], 'blue', 'LineWidth', 2);
        hold off
    end

    drawnow
    storeMovie(t) = getframe(fig);
end

%Display the plots
[height, width, p] = size(storeMovie(1).cdata);
fig2 = figure;
axis off
set(fig2, 'position', [50 50 width height])
buttonReplay = uicontrol;
buttonReplay.String = 'Replay';
buttonReplay.Callback = @(~, ~) movie(fig2, storeMovie);
movie(fig2, storeMovie)