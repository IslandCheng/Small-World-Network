function dif_i=gradient_f(x,i)
s=0;
for j=1:50
    s=s+x(j);
end
s=s-x(i);
a_i=50;
c=5;
d=8;
b_i=50-0.2*i;
dif_i=2*a_i*(x(i)-b_i)+2*c*x(i)+c*s+d;
end