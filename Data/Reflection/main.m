close all; clear; clc;
Data 	= 	load('Data.dat');

figure()
hold on
plot(Data(:,1),Data(:,2).^2)
plot(Data(:,1),Data(:,3).^2)
hold off
ylim([0 1])

figure()
hold on
plot(Data(:,1),Data(:,4).^2)
plot(Data(:,1),Data(:,5).^2)
hold off
ylim([0 1])