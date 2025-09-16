function  drawgroundTruth(groundTruth)
% °´Ê±¿Ì»­³ögroundTruth

N = size(groundTruth);
figure(1);
clf(1);
box on;
axis([-200,200,-200,200]);
for  i = 1:N
    [i]
    axis([-200,200,-200,200]);
    hold on;
    cx = groundTruth{i}.x(1,:)';
    cy = groundTruth{i}.x(2,:)';
    P = groundTruth{i}.X;
    num = size(cx,1);
    for j = 1 : num
        n=2;
        if size(P(:,:,j),1) == 2
            sqrtP = n*sqrtm_2by2(P(:,:,j));
        else
            sqrtP = n*sqrtm(P(:,:,j));
        end
        N2 = 20;
        phi = 0:pi/N2:2*pi;
        xy = sqrtP*[cos(phi); sin(phi)];
        xy = xy + [cx(j)*ones(1,length(phi)) ; cy(j)*ones(1,length(phi))];

        x = xy(1,:)';
        y = xy(2,:)';
        %    ºì   ÂÌ  À¶ ÇàÂÌ Ñóºì »Æ  ºÚ  °×
        h{j} = plot(x,y,'k','linewidth',1);
        plot(cx,cy,'r.')
%         h = fill(x,y,'g');
%         set(h,'facealpha',0.1);
    end
    pause(0.1);
    for j=1:num
        delete(h{j});
    end
    
end

end

