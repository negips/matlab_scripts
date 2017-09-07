function [bigmass bigconv bigconvd bigconvxd velvec gno nreps nn] = AssembleBig(El,Nx,Ny,nelx,nely,nelv,ifplots)

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

     bigmass = zeros(nn,nn);
     bigconv = zeros(nn,nn);
     bigconvd = zeros(nn,nn);
     bigconvxd = zeros(nn,nn);

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
               end
          end
     end

     if ifplots     
          figure
          spytol(El(1).convalld,6,'r')

          h=figure;
          spytol(bigconvd,6)

          filename = 'spy_dealiased_N4_nelv4_uniform';
          destn = './plots';
%         SaveFig(h, filename, destn, 1)
     end
%%   Eigenvalues:
     deltat = 0.001;
     e1 = eig(inv(bigmass)*bigconvd);
     lambda_deltat = e1;
%     clines = load('bdfk-neutral-curve.mat');

     lambdar = real(e1)*deltat;
     lambdai = imag(e1)*deltat;

     if ifplots
          hstab= figure;
          plot(lambdar,lambdai,'*', 'MarkerSize', 12)
          hold on
%          plot(clines.cline3(1,2:end),clines.cline3(2,2:end), 'r')

          %xlim([min(lambdar) max(lambdar)]);
          %ylim([min(lambdai) max(lambdai)]);
     end
%%

     return
