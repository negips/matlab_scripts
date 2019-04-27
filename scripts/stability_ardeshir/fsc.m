function dy = fsc(t,y,mfac)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dy=zeros(5,1);
betaH=2*mfac/(mfac+1);
dy(1)=y(2);
dy(2)=y(3);
dy(2)=y(3);
dy(3)=-(y(1)*y(3)+betaH*(1-y(2)*y(2)));
dy(4)=y(5);
dy(5)=-y(1)*y(5);
end

