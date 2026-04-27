close all
plot(tsim, x_nonl)
hold on
plot(tsim, x_ref)
xlabel("Time (s)")
ylabel("Location (m)")
legend("X Position","X Reference")
figure
plot(tsim, theta_nonlin)
xlabel("Time (s)");
ylabel("Angle (rad)");