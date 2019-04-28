function []=plot_OS(eigval,eigvec,y,header)
% []=plot_OS(eigval,eigvec,y,header)
% This function plots the spectrum and the corresponding eigenvectors 
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

% Choose eigenvalue to plot
disp('Right click to pick the mode. Left click to stop');
button=1;
while button==1
    
    figure(1)
    plot(real(eigval),imag(eigval),'o',[0 1],[0 0],'r-');
    ylim([-1,0.1]);
    title(header);
    ylabel('imag(\omega)')
    xlabel('real(\omega)')
         
    [xp,yp,button] = ginput(1);
    a=xp+sqrt(-1)*yp;
    [c,locate]=min(abs(eigval-a));
    
    v=eigvec(:,locate);
    figure(2)
    plot(real(v),y,imag(v),y);
    legend('real(v)','imag(v)')
    title(strcat(header,{', omega= '},{num2str(eigval(locate))}));
    ylabel('y')

end
end



