function [conn] = ElemConnectivity(elno,nelx,nely,ii,jj,xperiodic,yperiodic)

%         Element numbers of neighbors
%         [left bottom right top]

          conn = zeros(1,4);

%         Connectivity
%         Along x
          if ii>1
               conn(1) = elno-1;
          else
               if xperiodic
                    conn(1) = elno+nelx-1;
               else
                    conn(1) = 0;
               end
          end
          if ii<nelx
               conn(3) = elno+1;
          else
               if xperiodic
                    conn(3) = elno-nelx+1;
               else
                    conn(3) = 0;
               end
          end

%         Along y
          if jj>1
               conn(2) = elno-nelx;
          else
               if yperiodic
                    conn(2) = (nely-1)*nelx+ii;
               else
                    conn(2) = 0;
               end
          end
          if jj<nely
               conn(4) = elno+nelx;
          else
               if yperiodic
                    conn(4) = ii;
               else
                    conn(4) = 0;
               end
          end
          return

%%        End of connectivity       
