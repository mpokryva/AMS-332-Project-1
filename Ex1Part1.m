close
n = 100;
s = (0 : 1 : n);
h = [1, 2, 10];
vMax = 5.0; % 5.0 mM/s
k = 20.0; % 20.0 mM
y = zeros(3, n); % 3 x 100 array
for i = 1 : length(h)
    for j = 1 : length(s)
        sh = s(j) .^ h(i);
        y(i, j) = (vMax * sh) / (k .^ h(i) + sh);
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
legend("h=1", "h=2", "h=10");
    