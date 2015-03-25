load('lightertraining.mat')
densi = reshape(densipoints,3,10,[]);
r = densi(1,:,:);
g = densi(2,:,:);
b = densi(3,:,:);
[n1, xout1] = hist(r(:));
[n2, xout2] = hist(g(:));
[n3, xout3] = hist(b(:));
figure(1);subplot(1,3,1);
bar(xout1, n1, 'r');
a = axis;
subplot(1,3,2);
bar(xout2, n2, 'g');
a = max(a,axis);
subplot(1,3,3);
bar(xout3, n3, 'b');
a = max(a,axis);

for i = 1:3
subplot(1,3,i);
axis(a);
end