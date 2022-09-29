
initial = readmatrix('phase_distribution_40_reverse.txt');
second = readmatrix('ph40_clean.txt');
figure;
hold on
plot(initial(:,1), initial(:,2))
plot(second(:,1), second(:,2))
