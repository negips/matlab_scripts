function [bigmass bigconv bigconvd bigconvxd bigforc bigconvxd_new bigconvyd_new bigconvd_new biglapl velvec gno nreps El] = AssembleBig6(El,Nx,Ny,nelx,nely,nelv,ifplots)

%%   Build combined matricies

     [lx ly] = size(El(1).lglmap);

     nbasis = lx*ly;

     vecn = zeros(nbasis,nelv);
     gno = zeros(nbasis,1);
     velvec = zeros(nbasis,1);
     nreps = zeros(nbasis,1);

     nn = 0;
     for elno=1:nelv
          if elno==1
               n2=0;
               for jj=0:Ny
                    for ii=0:Nx
                         n2=n2+1;
                         nn=nn+1;
                         vecn(n2,elno) = n2;
                         gno(nn) = El(elno).lglmap(ii+1,jj+1);
                         nreps(nn) = 1;
          
                         velvec(nn)=El(elno).un(ii+1,jj+1);
                    end
               end
               continue;
          end

          n3=0;
          for jj=0:Ny
               for ii=0:Nx
                    n3=n3+1;
                    glno = El(elno).lglmap(ii+1,jj+1);
                    ix = find(abs(gno) == glno);
                    if isempty(ix)
                         nn=nn+1;
                         vecn(n3,elno) = nn;
                         gno(nn) = glno;
                         nreps(nn) = 1;

                         velvec(nn) = El(elno).un(ii+1,jj+1);
                    else
                         vecn(n3,elno) = -ix;
                         nreps(ix) = nreps(ix) + 1;
                    end

               end
          end
     
     end

     for elno=1:nelv
       El(elno).TotDegFreedom = nn;
     end  

     bigmass = zeros(nn,nn);
     bigconv = zeros(nn,nn);
     bigconvd = zeros(nn,nn);
     bigconvxd = zeros(nn,nn);
     bigforc = zeros(nn,nn);
     bigconvxd_new = zeros(nn,nn);
     bigconvyd_new = zeros(nn,nn);
     bigconvd_new = zeros(nn,nn);
     biglapl = zeros(nn,nn);

     for elno=1:nelv
          for ii=1:nbasis
               ix = abs(vecn(ii,elno));
               for jj=1:nbasis
                    iy = abs(vecn(jj,elno));
                    entry = El(elno).convall(ii,jj);
                    bigconv(ix,iy) = bigconv(ix,iy)+entry;

                    entry = El(elno).convxd(ii,jj);
                    bigconvxd(ix,iy) = bigconvxd(ix,iy)+entry;

                    entry = El(elno).convalld(ii,jj);
                    bigconvd(ix,iy) = bigconvd(ix,iy)+entry;

                    entry = El(elno).mass(ii,jj);
                    bigmass(ix,iy) = bigmass(ix,iy)+entry;

                    entry = El(elno).forc(ii,jj);
                    bigforc(ix,iy) = bigforc(ix,iy)+entry;

                    entry = El(elno).convxd_new(ii,jj);
                    bigconvxd_new(ix,iy) = bigconvxd_new(ix,iy)+entry;

                    entry = El(elno).convyd_new(ii,jj);
                    bigconvyd_new(ix,iy) = bigconvyd_new(ix,iy)+entry;

                    entry = El(elno).convalld_new(ii,jj);
                    bigconvd_new(ix,iy) = bigconvd_new(ix,iy)+entry;

                    entry = El(elno).lapl(ii,jj);
                    biglapl(ix,iy) = biglapl(ix,iy)+entry;
               end
          end
     end


     return
