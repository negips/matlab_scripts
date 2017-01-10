% scr 
function varargout = reshapenek3D(data,nelx,nely,nelz)
% function varargout = reshapenek3D(data,nelx,nely)

%
% Reshape data from nekread to a meshgrid
%
%   [meshgrid1,meshgrid2,...] = reshapenek(data,nelx,nely)
%

% get dimension and check number of elements
[nel,N2,nfld] = size(data); N = int16(N2^(1/3));

if nel ~= nelx*nely*nelz
    disp('Error: nel ~= nelx*nely*nelz.');
%     disp('Error: nel ~= nelx*nely.');

    return
end

% check output
if nfld < nargout
    disp('Error: nfld < outputs.');
    return
end

% reshape data
for ifld = 1:min([nfld,nargout])
    mesh = zeros((N-1)*nelz+1,(N-1)*nely+1,(N-1)*nelx+1);
%     mesh = zeros((N-1)*nely+1,(N-1)*nelx+1);
   nel2D =nelx*nely; 
    for iel = 1:nel
        
        zslice = floor((iel-1)/nel2D);
        
        ielx = floor((iel-1-zslice*nel2D)/nely) + 1;
        iely = mod((iel-zslice*nel2D)-1,nely) + 1;
        ielz = zslice + 1;

        ii = (0:N-1) + (N-1)*(ielx-1) + 1;
        jj = (0:N-1) + (N-1)*(iely-1) + 1;
        kk = (0:N-1) + (N-1)*(ielz-1) + 1;
        
        mesh(kk,jj,ii) = reshape(data(iel,:,ifld),N,N,N);
        
    end
    
    varargout{ifld} = mesh;
end