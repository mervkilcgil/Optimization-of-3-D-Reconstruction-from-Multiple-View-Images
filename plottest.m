clear all
close all
clc
myFolder = 'C:\Users\e1943\Documents\merve\iam566\566 project\matlab';
ef1 = fullfile(myFolder, 'energyfunctions_2.txt');
Tef1 = readtable(ef1);
tef1 = Tef1{1:61,1:end};


e1 = fullfile(myFolder, 'errors_2.txt');
Te1 = readtable(e1);
te1 = Te1{1:61,1:end};


iteration1 = linspace(1,length(tef1),length(tef1));


figure (10)
plot(iteration1, tef1)
title('Energy function change for 1 image');
xlabel('iteration number')
ylabel('Energy function value')
set(gca, 'YScale', 'log')
figure (11)
plot(iteration1,te1)
title('Error change for 1 image');
xlabel('iteration number')
ylabel('Error')
set(gca, 'YScale', 'log')
