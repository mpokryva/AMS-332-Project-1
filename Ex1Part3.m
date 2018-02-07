close
n = 100;
s = (0 : 1 : n);
h = 2;
vMax = [2.0, 5.0, 10.0]; % 5.0 mM/s
k = 20.0; % 20.0 mM
y = zeros(3, n); % 3 x 100 array
for i = 1 : length(vMax)
    for j = 1 : length(s)
        sh = s(j) .^ h;
        y(i, j) = (vMax(i) * sh) / (k .^ h + sh);
        fprintf("y(%d, %d): %d\n", i, j, y(i, j));
    end
end
for i = 1 : size(y, 1)
  hold on
  plot(s, y(i, :));
  hold off
end
title("Reaction Velocity vs Concentration of Autoregulatory Gene");
xlabel("Substrate concentration [S] (mM)");
ylabel("Reaction velocity (mM/s)");
legend("vMax=2.0", "vMax=5.0", "vMax=10.0");
    