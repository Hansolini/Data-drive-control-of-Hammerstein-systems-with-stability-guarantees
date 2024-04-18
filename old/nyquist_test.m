s = tf('s');

G = 1/(s^2 + s +1);
K1 = 1 + 1/s;
K2 = 3 + s;



[A1,A2] = meshgrid(linspace(1,10,30));
Y1 = zeros(size(A1));
Y2 = zeros(size(A1));

for i = 1:size(A1,1)
    for j = 1:size(A1,2)
        Y1(i,j) = max(real(zero((A1(i,j)*K1 + A2(i,j)*K2)*G + 1)));
        Y2(i,j) = max(real(zero(max(A1(i,j),A2(i,j))*(K1 + K2)*G + 1)));
    end
end
%%
figure(1);
clf; grid on; hold on;
surf(A1,A2,Y2 - Y1,'FaceColor','interp')
surf(A1,A2,0*Y2,'FaceColor','red')


%view([45,15])

%%
G1 = s^2 + 2*s + 1;
G2 = s^2 + s + 3;
a = 2; b = 3;
clc;
disp(max(real(zero(G1 + a*s + b)))')
disp(max(real(zero(G1 + b*s + b)))')
disp(max(real(zero(G1 + a*s + b))) < max(real(zero(G1 + b*s + b))))
disp('------------------------')
disp(max(real(zero(G2 + a*s + b)))')
disp(max(real(zero(G2 + b*s + b)))')
disp(max(real(zero(G2 + a*s + b))) < max(real(zero(G2 + b*s + b))))


