locationerr1 = zeros(T,1);
locationerr2 = zeros(T,1);
% locationerr3 = zeros(T,1);
locationerr4 = zeros(T,1);
T=100;
for j = 1:T
    locationerr1 = locationerr1+[Decomposed_cost1{j}.localisation]';
    locationerr2 = locationerr2+[Decomposed_cost2{j}.localisation]';
    % locationerr3 = locationerr3+[Decomposed_cost3{j}.localisation]';
    locationerr4 = locationerr4+[Decomposed_cost4{j}.localisation]';
end
locationerr1 = locationerr1/T;
locationerr2 = locationerr2/T;
% locationerr3 = locationerr3/T;
locationerr4 = locationerr4/T;
figure(14);
clf(14);

hold on;
plot(1:T,locationerr1,"r","LineWidth",1.5);
plot(1:T,locationerr2,"b","LineWidth",1.5);
% plot(1:T,locationerr3,"g","LineWidth",1.5);
plot(1:T,locationerr4,"k--","LineWidth",1.5);


