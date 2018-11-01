%     Testing the reading of the re2 file

      clear
      clc
      close all

      disp('Test reading re2')

      endian= 'le';
      fname2='new.re2';
      [fre2,message] = fopen(fname2,'r',['ieee-' endian]);
      if fre2 == -1
        disp(message) 
        return 
      end

      cr = fread(fre2,80,'char')';
      cr = char(cr)

      etag=fread(fre2,1,'float32')      

      
      nelg=3344;
      ndim=2;
      xy=zeros(2^ndim,2,nelg);
%     skip x,y,z      
      for i=1:nelg
        grp=fread(fre2,1,'float64');
        for j=1:ndim
          data=fread(fre2,2^ndim,'float64')';
          xy(:,j,i)=data;
        end
      end

%      for i=1:1
%        x=xy(:,1,i);
%        x=x(:);
%        y=xy(:,2,i);
%        y=y(:);
%        plot(x,y, '.'); hold on
%        pause(1.01)
%      end  

%     ncurves
      ncurve=fread(fre2,1,'float64')'
      for i=1:ncurve
        par{1}=fread(fre2,1,'float64');%el
        par{2}=fread(fre2,1,'float64');%edge
        par{3}=fread(fre2,1,'float64');%p1
        par{4}=fread(fre2,1,'float64');%p2
        par{5}=fread(fre2,1,'float64');%p3
        par{6}=fread(fre2,1,'float64');%p4
        par{7}=fread(fre2,1,'float64');%p5
        ctype=fread(fre2,1,'float64');%ctype
        par{8}=char(ctype);
      end  

%     nbc
      nbc=fread(fre2,1,'float64')'

      for i=1:nbc
        par{1}=fread(fre2,1,'float64');%el
        par{2}=fread(fre2,1,'float64');%face
        par{3}=fread(fre2,1,'float64');%connectsto
        par{4}=fread(fre2,1,'float64');%onface
        par{5}=fread(fre2,1,'float64');%p1
        par{6}=fread(fre2,1,'float64');%p2
        par{7}=fread(fre2,1,'float64');%p3
        bc=fread(fre2,8,'char*1');%bc
        par{8}=char(bc);
      end  

      fclose(fre2);
