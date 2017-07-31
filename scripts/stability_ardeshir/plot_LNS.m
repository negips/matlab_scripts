function []=plot_LNS(eigval,eigvec,y,header)
% []=plot_LNS(eigval,eigvec,y,header) 
% This function plots the spectrum and selected eigenfunctions
%
% INPUT
%   eigval: array containing the eigenvalues
%   eigvec: matrix corresponding the eigenvectors
%   y: normal coordinates
%   header: string to appear as the header of plot
%
% (c) Ardeshir Hanifi 2014
%
eglist=find(abs(eigval)<10);
eigval=eigval(eglist);
eigvec=eigvec(:,eglist);

N=length(y);

% Choose eigenvalue to plot

disp('Right click to pick the mode. Left click to stop');
button=1;
while button==1
    
    figure(1)
    plot(real(eigval),imag(eigval),'o',[0 1.2],[0 0],'r-');
    %ylim([-1,0.1]);
%    xlim([0 0.2]);
%    ylim([-0.2 0])   
    title(header);
    ylabel('imag(\omega)')
    xlabel('real(\omega)')
    
    ylim([-1,0.1]);
    [xp,yp,button] = ginput(1);
    a=xp+sqrt(-1)*yp;
    [c,locate]=min(abs(eigval-a));
    
    u=eigvec(1:N,locate);
    v=eigvec(N+1:2*N,locate);
    w=eigvec(2*N+1:3*N,locate);
    p=eigvec(3*N+1:4*N,locate);
    
    figure(2)
%    plot(abs(u),y,abs(v),y,abs(w),y,abs(p),y);
%    legend('abs(u)','abs(v)','abs(w)','abs(p)')
    plot(real(u),y,real(v),y,real(w),y,real(p),y);
    legend('u','v','w','p')
    title(strcat(header,{', omega= '},{num2str(eigval(locate))}));
    ylabel('y')

end

end



