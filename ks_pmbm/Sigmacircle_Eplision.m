function r = Sigmacircle_Eplision(cx,cy,P,N,model)

num = size(cx,1);
for i = 1 : num
    if size(P(:,:,i),1) == 2
        sqrtP = N*sqrtm_2by2(P(:,:,i));
    else
        sqrtP = N*sqrtm(P(:,:,i));
    end

    phi = linspace(0,2*pi,model.direction);
    xy = sqrtP*[cos(phi); sin(phi)];
    xy = xy + [cx(i)*ones(1,length(phi)) ; cy(i)*ones(1,length(phi))];

    x = xy(1,:)';
    y = xy(2,:)';

    A = [x';y'] - mean([x';y'],2);
    r_ = vecnorm(A, 2, 1);
    r.r = FFT2FS_new(r_);

    Cov_irregular(:,:,:) = (reshape(A, 2, 1, []) .* permute(reshape(A, 2, 1, []), [2, 1, 3]));
    r.shape = Cov_irregular;

%     %plot(x,y,'k-','linewidth',1);
%     plot(x,y,'-','color',[0.7 0.7 0.7],'linewidth',1);
%     fill(xy(1,:),xy(2,:),[0.9 0.9 0.9]);

%     c = ['r' 'g' 'b' 'c' 'm' 'y' 'k' 'w'];
%     %    ºì   ÂÌ  À¶ ÇàÂÌ Ñóºì »Æ  ºÚ  °×
%     plot(x,y,strcat('-',c(mod(color-1,7)+1)),'linewidth',1);
%     hold on;
%     pbaspect([1 1 1]);
%     axis([-200 200 -200 200])
end
end
