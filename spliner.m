	
X = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
Y = [0.01 -0.01 0.02 0.01 0 0.23 0.42 0.6 0.95];

N = length(X);
h = (X(N) - X(1)) / (N - 1);
M = diag(ones(1, N - 3), 1) ...
  + diag(4 * ones(1, N - 2))...
  + diag(ones(1, N - 3), - 1);
for i = 1:N - 2
    z(i) = (6 / h ^ 2) * (Y(i + 2) - 2 * Y(i + 1) + Y(i));
end
w = M \ z';

sigma = [0; w; 0];

for i = 1:N - 1
    a(i) = (sigma(i + 1) - sigma(i)) / (6 * h);
    b(i) = sigma(i) / 2;
    c(i) = (Y(i + 1) - Y(i)) / h ...
         - h / 6 * (2 * sigma(i) + sigma(i + 1));
    d(i) = Y(i);
end

r = 21;
hh = h / r; %size of subinterval

x = X(1):hh:X(N);

W = [];
for i = 1:N - 1
    for j = r * (i - 1) + 1:r * i
        W(j) = a(i) * (x(j) - X(i)) ^ 3 ...
             + b(i) * (x(j) - X(i)) ^ 2 ...
             + c(i) * (x(j) - X(i)) + d(i);
    end
end
W(r * (N - 1) + 1) = Y(N);
grid on
plot(X, Y, 'o')
hold on
plot(x, W)
hold off