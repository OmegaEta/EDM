function [x,y] = irregularcircle(r,m,theta,color)

if ~isempty(r{1})
figure(1)
theta = [theta []];
r_ = evaluateFourierSeries(r{1}.a(1:10),r{1}.b(1:10),theta);
[x,y] = pol2cart(theta,r_);

c = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'w'];
plot(x+m(1),y+m(2),strcat('-',c(mod(color-1,7)+1)),'linewidth',1,'HandleVisibility','off');
hold on;
else
    x=[];
    y=[];
end

end