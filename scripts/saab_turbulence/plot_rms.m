%     Plotting saved rms values

      clear
      clc
      close all

      load rmsuv
      load rmsyn
%     x0 = [0.4 0.5 0.6 0.7 0.75 0.8 0.85 0.9];
      kred=0.4;
      U0=1.0;
      chord=1;
      semichord=chord/2;
      omega=kred*U0/semichord;
      Tosc=2*pi/omega;
      ptch_start=6.0;
      ptch_amp=1.0;
      ini_aoa=3.4;
      ini_phase=-pi/2;

      nt = length(T);

      ind = 6;
      gray = [0.8 0.8 0.8];
      val  = [];
      yind = [1:40];
      tstamps = [];
      figure(1)

      ic=0;
      for i=1:nt
        clc
        i

        if i>20
%          delete(pl(i-20))
        end  
        t=T(i);
        theta(i) = ini_aoa + ptch_amp*sin(omega*(t-ptch_start)+ini_phase);
        omg   = omega*ptch_amp*cos(omega*(t-ptch_start)+ini_phase);

        if  theta(i)>=3.9 && theta(i)<=4.1 && t>32 && omg>=0
          ic=ic+1;
          tstamps = [tstamps t];
          y    = yn{i}(yind,ind);
          vt   = uv{i}(yind,ind);

          if mod(i,1)==0
            pl(i) = plot(y,vt,'Color', gray); hold on
            pause(0.01);
          end  

          if ic==1
            val = vt;
          else
            val = val + vt;
          end
        end  

      end  

      val = val/ic;

      plot(y,val, 'Color', 'k', 'LineWidth', 2)


      figure
      plot(T,theta, '.')





