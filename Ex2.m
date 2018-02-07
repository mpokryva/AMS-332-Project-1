close
%%% Part 1 %%%
totalTime = 20;
deltaT = 0.01;
[mu, omega, chiProt, chiRna] = deal(1); % (s^-1)
k = 0.33; % (mM)
initProt = 0.5;
initRna = 0.5;
figureNum = 1;
[time, prot, rna] = protRnaConcent(totalTime, deltaT, initProt, initRna, ...
    chiProt, chiRna, mu, omega, k);
plotProtAndRna(time, prot, rna, figureNum, "Protein=" + initProt + " " + "RNA=" + initRna);
figureNum = figureNum + 1;

%%% Part 2 %%%


init1 = 0.0;
init2 = [0.2, 0.5];
[time, prot, rna] = protRnaConcent(totalTime, deltaT, initProt, initRna, ...
    chiProt, chiRna, mu, omega, k);

% Plot prot = 0.
for i = 1 : length(init2)
   [time, prot, rna] = protRnaConcent(totalTime, deltaT, init1, init2(i), ...
       chiProt, chiRna, mu, omega, k);
plotProtAndRna(time, prot, rna, figureNum, "Protein=" + init1 + " RNA=" + init2(i));
figureNum = figureNum + 1;
end
   
% Plot rna = 0.
for i = 1 : length(init2)
   [time, prot, rna] = protRnaConcent(totalTime, deltaT, init2(i), init1, ...
       chiProt, chiRna, mu, omega, k);
    plotProtAndRna(time, prot, rna, figureNum, "Protein=" + init2(i) + " RNA=" + init1);
    figureNum = figureNum + 1;
end


%%% Part 3 %%%

init = (0 : 0.2 : 1.4); % Initial concentrations.
prots = zeros(length(init), totalTime/deltaT + 1);
rnas = zeros(length(init), totalTime/deltaT + 1);
for i = 1 : length(init) % i = protein.
    for j = 1 : length(init) % j = rna.
        hold on
        [time, prots(i,:), rnas(j,:)] = protRnaConcent(totalTime, deltaT, init(i), init(j), ...
            chiProt, chiRna, mu, omega, k);
        figure(figureNum);
        plot(prots(i,:), rnas(j,:));
        xlabel("Protein Concentration (mM)");
        ylabel("RNA Concentration (mM)");
        hold off
    end
end


function y = dtProt(prot, rna, chiProt, omega)
    y = (omega * rna) - (chiProt * prot);
end

function y = dtRna(prot, rna, chiRna, mu, k)
    y = (mu * prot .^ 2) / (k .^ 2 + prot .^ 2) - (chiRna * rna);
end

function [time, prot, rna] = protRnaConcent(totalTime, deltaT, initProt, ... 
initRna,chiProt, chiRna, mu, omega, k)
    prot = zeros(1, totalTime/deltaT);
    rna = zeros(1, totalTime/deltaT);
    prot(1) = initProt;
    rna(1) = initRna;
    time = zeros(1, totalTime/deltaT);
    % Forward Euler Algorithm
    for i = 1 : totalTime/deltaT
        deltaProt = dtProt(prot(i), rna(i), chiProt, omega) * deltaT;
        prot(i+1) = prot(i) + deltaProt;
        deltaRna = dtRna(prot(i), rna(i), chiRna, mu, k) * deltaT;
        rna(i+1) = rna(i) + deltaRna;
        time(i+1) = time(i) + deltaT;
    end
end

function plotProtAndRna(time, prot, rna, figureNum, plotTitle)
    figure(figureNum);
    plot(time, prot, time, rna);
    X_LABEL = "Time (s)";
    Y_LABEL = "Concentration (mM)";
    xlabel(X_LABEL);
    ylabel(Y_LABEL);
    labels = ["Protein", "RNA"];
    legend(labels);
    title(plotTitle);
end