function drawGGIW(GGIW,color,fig)
% clarify
% [x,y] = Sigmacircle(cx,cy,P,N,color)
d=2;
H = [1 0 0 0;0 1 0 0];
cx = GGIW.m(1);
cy = GGIW.m(2);
P = GGIW.V/(GGIW.v - 2*d -2) ;%+ H*GGIW.P*H';
N = 2;
N2 = 20;

num = size(cx,1);
for i = 1 : num
    if size(P(:,:,i),1) == 2
        sqrtP = N*sqrtm_2by2(P(:,:,i));
    else
        sqrtP = N*sqrtm(P(:,:,i));
    end

    phi = 0:pi/N2:2*pi;
    xy = sqrtP*[cos(phi); sin(phi)];
    xy = xy + [cx(i)*ones(1,length(phi)) ; cy(i)*ones(1,length(phi))];

    x = xy(1,:)';
    y = xy(2,:)';
    
    c = ['r' 'g' 'b' 'c' 'm' 'y' 'k'];
    %    ºì   ÂÌ  À¶ ÇàÂÌ Ñóºì »Æ  ºÚ 
    figure(fig);
    plot(x,y,strcat('-',c(mod(color-1,7)+1)),'linewidth',1);
    hold on;

end
end

