s = 38;
t = linspace(-2, 30, 500);
rho0 = 999.842594 + 6.793952e-2 .* t - 9.095290e-3 .* t.^2 + ...
       1.001685e-4 .* t.^3 - 1.120083e-6 .* t.^4 + 6.536336e-9 .* t.^5;
A = 8.24493e-1 - 4.0899e-3 .* t + 7.6438e-5 .* t.^2 - ...
    8.2467e-7 .* t.^3 + 5.3875e-9 .* t.^4;
B = -5.72466e-3 + 1.0227e-4 .* t - 1.6546e-6 .* t.^2;
C = 4.8314e-4;
rho = rho0 + A .* s + B .* s.^(3/2) + C .* s.^2;
figure;
plot(t, rho, 'b-', 'LineWidth', 1.5); 
grid on;
xlabel('Temperature (T) [°C]','FontSize', 28);
ylabel('Density (\rho) [kg/m^3]','FontSize', 28);
xlim([min(t), max(t)]);
ylim([min(rho)-0.1, max(rho)+0.1]);
