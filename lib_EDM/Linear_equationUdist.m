% Initialize, where each row of (a, b, c) represents a line.
% The function order is ax + by + c > 0, < 0, = 0.
% There are a total of three lines: two parallel lines and one line that passes through the 'center' point and is perpendicular to them.
% 'dist' represents the distance between the two parallel lines.
% 'center' represents the Object's position center.
% 'angle_' represents the angle at which the two parallel lines are located.
% 'block' represents the number of regions.

function [range_linear] = Linear_equationUdist(blokdist,center,angle)
    block=length(angle);
    range_linear.a=zeros(3,block);range_linear.b=zeros(3,block);range_linear.c=zeros(3,block);

   % Initialize, where each row of (a, b, c) represents a line.
    range_linear_.a=zeros(3,1);range_linear_.b=zeros(3,1);range_linear_.c=zeros(3,1);
    plusorminus=[];

    for num_block=1:block
        angle_ = mod(angle(num_block),2*pi);
       
        %Baseline on the axis.,center(x0,y0)。
        if angle_==pi/2
            %(1) x-(x0-d) > 0
            %(2) x-(x0+d) < 0
            %(3)        x > y0
            range_linear_.a(1)=1;range_linear_.b(1)=0;
            range_linear_.c(1)=-center(1)+blokdist/2;
            range_linear_.a(2)=1;range_linear_.b(2)=0;
            range_linear_.c(2)=-center(1)-blokdist/2;
            range_linear_.a(3)=0;range_linear_.b(3)=1;
            range_linear_.c(3)=-center(2);
            %plusorminus表示三个式子的正负情况，0为正，1为负
            plusorminus=[0,1,0];
%                 equal2
        elseif angle_==pi
            %(1) y-(y0-d) > 0
            %(2) y-(y0+d) < 0
            %(3)        y < x0
            range_linear_.a(1)=0;range_linear_.b(1)=1;
            range_linear_.c(1)=-center(2)+blokdist/2;
            range_linear_.a(2)=0;range_linear_.b(2)=1;
            range_linear_.c(2)=-center(2)-blokdist/2;
            range_linear_.a(3)=1;range_linear_.b(3)=0;
            range_linear_.c(3)=-center(1);
            plusorminus=[0,1,1];
%                 equal2
        elseif angle_==(3/2)*pi
            %(1) x-(x0-d) > 0
            %(2) x-(x0+d) < 0
            %(3)        x < y0
            range_linear_.a(1)=1;range_linear_.b(1)=0;
            range_linear_.c(1)=-center(1)+blokdist/2;
            range_linear_.a(2)=1;range_linear_.b(2)=0;
            range_linear_.c(2)=-center(1)-blokdist/2;
            range_linear_.a(3)=0;range_linear_.b(3)=1;
            range_linear_.c(3)=-center(2);
            plusorminus=[0,1,1];
%                 equal2
        elseif angle_==2*pi || angle_==0
            %(1) y-(y0-d) > 0
            %(2) y-(y0+d) < 0
            %(3)        y > x0
            range_linear_.a(1)=0;range_linear_.b(1)=1;
            range_linear_.c(1)=-center(2)+blokdist/2;
            range_linear_.a(2)=0;range_linear_.b(2)=1;
            range_linear_.c(2)=-center(2)-blokdist/2;
            range_linear_.a(3)=1;range_linear_.b(3)=0;
            range_linear_.c(3)=-center(1);
            plusorminus=[0,1,0];
%                 equal2
        end
    
        %Reference line off the axis.
        if angle_>0 && angle_<pi/2
            %1>0,2<0,3>0
            %(1) ax1+y1+c1>0
            %(2) ax2+y2+c2<0
            %(3) ax3+y3+c3>0
            range_linear_.a(1)=-tan(angle_);range_linear_.b(1)=1;
            range_linear_.c(1)=tan(angle_)*(center(1)+sin(angle_)*(blokdist/2))-(center(2)-cos(angle_)*(blokdist/2));
            range_linear_.a(2)=-tan(angle_);range_linear_.b(2)=1;
            range_linear_.c(2)=tan(angle_)*(center(1)-sin(angle_)*(blokdist/2))-(center(2)+cos(angle_)*(blokdist/2));
            range_linear_.a(3)=-1/(-tan(angle_));range_linear_.b(3)=1;
            range_linear_.c(3)=-(-1/-tan(angle_))*center(1)-center(2);
            plusorminus=[0,1,0];
%                 equal1
        elseif angle_>pi/2 && angle_<pi
            %(1) ax1+y1+c1>0
            %(2) ax2+y2+c2<0
            %(3) ax3+y3+c3>0
            range_linear_.a(1)=-tan(angle_);range_linear_.b(1)=1;
            range_linear_.c(1)=tan(angle_)*(center(1)-cos(angle_-pi/2)*(blokdist/2))-(center(2)-sin(angle_-pi/2)*(blokdist/2));
            range_linear_.a(2)=-tan(angle_);range_linear_.b(2)=1;
            range_linear_.c(2)=tan(angle_)*(center(1)+cos(angle_-pi/2)*(blokdist/2))-(center(2)+sin(angle_-pi/2)*(blokdist/2));
            range_linear_.a(3)=-1/-tan(angle_);range_linear_.b(3)=1;
            range_linear_.c(3)=-(-1/-tan(angle_)*center(1)+center(2));
            plusorminus=[0,1,0];
%                 equal1
        elseif angle_>pi && angle_<(3/2)*pi
            %(1) ax1+y1+c1>0
            %(2) ax2+y2+c2<0
            %(3) ax3+y3+c3<0
            range_linear_.a(1)=-tan(angle_);range_linear_.b(1)=1;
            range_linear_.c(1)=tan(angle_)*(center(1)+sin(angle_-pi)*(blokdist/2))-(center(2)-cos(angle_-pi)*(blokdist/2));
            range_linear_.a(2)=-tan(angle_);range_linear_.b(2)=1;
            range_linear_.c(2)=tan(angle_)*(center(1)-sin(angle_-pi)*(blokdist/2))-(center(2)+cos(angle_-pi)*(blokdist/2));
            range_linear_.a(3)=-1/-tan(angle_);range_linear_.b(3)=1;
            range_linear_.c(3)=-(-1/-tan(angle_)*center(1)+center(2));
            plusorminus=[0,1,1];
%                 equal1
        elseif angle_>(3/2)*pi && angle_<2*pi
            %(1) ax1+y1+c1>0
            %(2) ax2+y2+c2<0
            %(3) ax3+y3+c3<0
            range_linear_.a(1)=-tan(angle_);range_linear_.b(1)=1;
            range_linear_.c(1)=tan(angle_)*(center(1)-cos(angle_-(3/2)*pi)*(blokdist/2))-(center(2)-sin(angle_-(3/2)*pi)*(blokdist/2));
            range_linear_.a(2)=-tan(angle_);range_linear_.b(2)=1;
            range_linear_.c(2)=tan(angle_)*(center(1)+cos(angle_-(3/2)*pi)*(blokdist/2))-(center(2)+sin(angle_-(3/2)*pi)*(blokdist/2));
            range_linear_.a(3)=-1/-tan(angle_);range_linear_.b(3)=1;
            range_linear_.c(3)=-(-1/-tan(angle_)*center(1)+center(2));
            plusorminus=[0,1,1];
%                 equal1
        end
        range_linear.a(1:3,num_block)=range_linear_.a;
        range_linear.b(1:3,num_block)=range_linear_.b;
        range_linear.c(1:3,num_block)=range_linear_.c;
        range_linear.plusorminus(:,:,num_block)=plusorminus;
    end
end

