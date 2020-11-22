un = [-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
in = [0.01 -0.01 0.02 0.01 0 0.23 0.42 0.6 0.95];

i = @(u) interpolate(u, un, in);
r = @(u) i(u) / u;

ud = [-1:0.01:1];
id = [];
ad3 = [];
ad5 = [];
spl = [];
for i = 1:length(ud)
    id(i) = ctv(interpolate(un, in), ud(i));
    ad3(i) = ctv(aprox(un, in, 3), ud(i));
    ad5(i) = ctv(aprox(un, in, 5), ud(i));
    hold on
end
plot(un, in, 'o', ud, id, ud, ad3, ud, ad5)
grid on

%interpolacja
function A = interpolate(X, Y)
    M = vandermonde(X);
    A = M \ Y';
end

%aproksymacja stopnia n
function A = aprox(X, Y, n)
    M = vandermonde(X);
    M = M(:, 1:n);
    MTM = M'*M;
    MTY = M'*Y';
    A = MTM \ MTY;
end

%macierz Vandermonde'a
function M = vandermonde(X)
    N = length(X);
    for w = 1:N
        for k = 1:N
            A(w, k) = X(w) ^ (k - 1);
        end
    end
    M = A;
end

%na podstawie macierzy współczynników i x wylicza wartość wielomianu W(x)
function y = ctv(A, x)
    s = 0;
    for i = 1:length(A)
        s = s + A(i) * x ^ (i - 1);
    end
    y = s;
end