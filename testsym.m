syms x y
y = cos(x);
D = diff(y,x);
E = subs(D,x,[1 1 1 1 1 1]);

syms f(x)
f(x) = x^2;
fnew = subs(f,x,2)