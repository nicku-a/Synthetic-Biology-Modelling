clc; clear all; close all;

IPTG = csvread('IPTG.csv');
ARAB = csvread('ARAB.csv');
sol = csvread('AND_gate_solution.csv');

surf(IPTG,ARAB,sol)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([10^(-7) 10^(-2)])
ylim([10^(-7) 10^(-2)])
xlabel('IPTG (M)')
ylabel('ARAB (M)')
zlabel('Norm. fluo/OD600')
colorbar