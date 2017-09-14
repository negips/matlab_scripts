function [gl_soln] = gs_gather(lx,periodic)

     [nx nz] = size(lx);
      dof = nz*(nx-1)+1;
      xg=zeros(dof,1);
      xl=0*lx;

      for nel=1:nz
        if nel==1
              xg(1:nx)=lx(:,nel);
              count=nx;
              xl=lx(:,nel);
        else
             xg(count)=(xg(count)+lx(1,nel))/2;
             xg(count+1:count+nx-1)=lx(2:end,nel);
             count=count+nx-1;
             xl(:,nel) = lx(:,nel);
             xl(1,nel) = (xl(1,nel)+lx(end,nel-1))/2;
             xl(end,nel-1)=xl(1,nel);
        end
      end
      if (periodic)
        xg(1)=(xg(1)+xg(end))/2;
        xg(end) = xg(1);
      end
      gl_soln=xg;

      if (periodic)
        xl(end,nz)=(xl(end,nz)+xl(1,1))/2;
        xl(1,1)=xl(end,nz);
      end

      return

