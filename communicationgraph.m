%communication graph
A=load('A.mat').A;

figure
G1=graph(A)
plot(G1,'XData',X,'YData',Y)
savefig('figure6.fig')