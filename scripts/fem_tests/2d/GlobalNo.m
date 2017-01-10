function [lglmap glnew] = GlobalNo(elno,Nx,Ny,nelx,nely,xperiodic,yperiodic,gltot,tmpEl)

     gllel = zeros(Nx+1,Ny+1);
     glno=gltot;
     for jj=0:Ny
     for ii=0:Nx
          if mod(elno,nelx)==1 && ii==0                % Left Boundary
               if elno==1                              % first row.
                    glno=glno+1;
                    lglmap(ii+1,jj+1)=glno;
               elseif (elno-(nelx*(nely-1))==1) && jj==Ny        % Last row. top left vertex.
                    connecti=tmpEl(elno).conn(4);
                    if connecti>0
                         lglmap(ii+1,jj+1)=tmpEl(connecti).lglmap(ii+1,1);
                    else
                         glno=glno+1;
                         lglmap(ii+1,jj+1)=glno;
                    end
               elseif jj==0                            % lower vertex
                    connecti=tmpEl(elno).conn(2);
                    if connecti>0
                         lglmap(ii+1,jj+1)=tmpEl(connecti).lglmap(ii+1,Ny+1);
                    end
               elseif jj>0
                    glno=glno+1;
                    lglmap(ii+1,jj+1)=glno;
               end
               continue;
          end

          if jj==0 && elno<=nelx                       % First row. Lower Boundary
               if ii==0                                % Lower vertex
                    connecti=tmpEl(elno).conn(1);
                    lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(Nx+1,jj+1);
               elseif ii==Nx && elno==nelx             % Lower right vertex of the domain.
                    connecti=tmpEl(elno).conn(3);
                    if connecti>0
                         lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(1,jj+1);
                    else
                         glno=glno+1;
                         lglmap(ii+1,jj+1)=glno;
                    end
               else    
                    glno=glno+1;
                    lglmap(ii+1,jj+1)=glno;
               end
               continue;
          end

          if jj==Ny && (elno-(nelx*(nely-1))>0)         % Last row. Upper Boundary
               connecti = tmpEl(elno).conn(4);
               if connecti>0
                    lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(ii+1,1);
               else
                    glno=glno+1;
                    lglmap(ii+1,jj+1) = glno;
               end
               continue;
          end

          if (ii==Nx) && (mod(elno,nelx)==0)           % Last Column. Right boundary
               if elno==nelx                           % First row
                    connecti = tmpEl(elno).conn(3);
                    if connecti>0
                         lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(1,jj+1);
                    else
                         glno=glno+1;
                         lglmap(ii+1,jj+1) = glno;
                    end
               elseif jj>0
                    connecti = tmpEl(elno).conn(3);
                    if connecti>0
                         lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(1,jj+1);
                    else
                         glno=glno+1;
                         lglmap(ii+1,jj+1) = glno;
                    end
               else
                    connecti = tmpEl(elno).conn(2);
                    if connecti>0
                         lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(ii+1,Ny+1);
                    else
                         glno=glno+1;
                         lglmap(ii+1,jj+1) = glno;
                    end
               end
               continue;
          end

          if ii==0
               connecti = tmpEl(elno).conn(1);
               lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(Nx+1,jj+1);
               continue;
          end

          if jj==0
               connecti = tmpEl(elno).conn(2);
               lglmap(ii+1,jj+1) = tmpEl(connecti).lglmap(ii+1,Ny+1);
               continue;
          end

          glno = glno +1;
          lglmap(ii+1,jj+1) = glno;
          
     end
     end
     glnew=glno;

     return
     

               
          
