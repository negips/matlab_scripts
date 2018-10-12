function rea = StupidChange(rea,shift)

      if shift<=0
        disp(['Invalid Shift: ', num2str(shift)])
        return
      end

      nelg=rea.mesh.nelg;
      ncurve=rea.mesh.Ncurve;
%     Shift curved sides      
      for i=1:ncurve
        cs=rea.mesh.curveedge(i);
        cs=cs+shift;
        if cs>4
          cs=cs-4;
        elseif cs<0
          cs=cs+4;
        end
        rea.mesh.curveedge(i)=cs;
      end

%     Shift connecting element information
      nfaces=4;
      for i=1:nelg
        on=zeros(nfaces,1);
        ind=zeros(nfaces,1);
        for j=1:nfaces
          of(j)=rea.mesh.cbc(j,i).onface+shift;
          if of(j)>4
            of(j)=of(j)-4;
          elseif of(j)<0
            of(j)=of(j)+4;
          end

          ind(j)=ind(j)+shift;
          if ind(j)>4
            ind(j)=ind(j)-4;
          elseif ind(j)<0
            ind(j)=ind(j)+4;
          end
        end

        for j=1:nfaces
          jj=ind(j);
          rea.mesh.cbc(jj,i).onface=of(j);
        end  
      end

end   %
