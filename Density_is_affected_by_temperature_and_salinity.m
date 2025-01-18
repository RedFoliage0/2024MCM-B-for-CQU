t = linspace(-2, 30, 100);
s = linspace(0, 40, 100);
[T, S] = meshgrid(t, s);
rho0 = 999.842594 + 6.793952e-2 .* T - 9.095290e-3 .* T.^2 + ...
       1.001685e-4 .* T.^3 - 1.120083e-6 .* T.^4 + 6.536336e-9 .* T.^5;
A = 8.24493e-1 - 4.0899e-3 .* T + 7.6438e-5 .* T.^2 - ...
    8.2467e-7 .* T.^3 + 5.3875e-9 .* T.^4;
B = -5.72466e-3 + 1.0227e-4 .* T - 1.6546e-6 .* T.^2;
C = 4.8314e-4;
rho = rho0 + A .* S + B .* S.^(3/2) + C .* S.^2;
figure;
surf(T, S, rho, 'EdgeColor', 'none');
colormap('parula');
alpha(0.6);
colorbar;
xlabel('Temperature (T) [°C]','FontSize', 28);
ylabel('Salinity (S) [PSU]','FontSize', 28);
zlabel('Density (\rho) [kg/m^3]','FontSize', 28);
grid on;
view(45, 30);
