function y = Ein(x)

gammac  = 0.5772156649015328606;
y       = -expint(x) - gammac - log(x);
y(x==0) = 0;