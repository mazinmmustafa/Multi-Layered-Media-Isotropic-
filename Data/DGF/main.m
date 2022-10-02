close all; clear; clc;
Data 	= 	load('Data.dat');

figure()
plot(Data(:,1),log10(Data(:,2)))

figure()
plot(Data(:,1),log10(Data(:,3)))

figure()
plot(Data(:,1),log10(Data(:,4)))

figure()
plot(Data(:,1),log10(Data(:,5)))

figure()
plot(Data(:,1),log10(Data(:,6)))

figure()
plot(Data(:,1),log10(Data(:,7)))

figure()
plot(Data(:,1),log10(Data(:,8)))

figure()
plot(Data(:,1),log10(Data(:,9)))

figure()
plot(Data(:,1),log10(Data(:,10)))