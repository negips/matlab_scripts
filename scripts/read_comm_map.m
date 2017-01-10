% Building Communication Map

clear
clc

filename = ['comm_map_00256.out'];
fid = fopen(filename);
s = str2double(filename(10:14));

sndmap = ones(s) -3;
recmap = ones(s) -3;
recnos = ones(s,1) -3;

cnt = 0;
snd_max = 0;
rec_max = 0;

while ~feof(fid)

    l = fgets(fid);
    cnt = cnt+1;   
    
    lin = mod(cnt,3);
    
    if lin == 1
        A = sscanf(l,'%*s %d %d %d %d');
        NID = A(1);
        SND_CNT = A(2);
        RECNO = A(3);
        MAXREC = A(4);
        if SND_CNT>snd_max
          snd_max = SND_CNT;
        end
        if RECNO>rec_max
          rec_max=RECNO;
        end        
                  
        recnos(NID+1) = RECNO;  
    elseif lin == 2
        readformat = '%*s %d';
        for i=2:SND_CNT
            readformat = [readformat ' %d'];
        end      
        PROCID = sscanf(l,readformat);

        l1 = length(PROCID);
        
        while l1<SND_CNT
            readformat = '%d';
            for i=2:(SND_CNT-l1)
                readformat = [readformat ' %d'];
            end
            l2 = fgets(fid);
            A2 = sscanf(l2,readformat);
            PROCID = [PROCID; A2];
            l1 = length(PROCID);
        end           
          

    elseif lin==0
        readformat = '%*s %d';
        for i=2:SND_CNT
            readformat = [readformat ' %d'];
        end         
        PROCPOS = sscanf(l,readformat);
        l1 = length(PROCPOS);
        
        while l1<SND_CNT
            readformat = '%d';
            for i=2:(SND_CNT-l1)
                readformat = [readformat ' %d'];
            end
            l2 = fgets(fid);
            A2 = sscanf(l2,readformat);
            PROCPOS = [PROCPOS; A2];
            l1 = length(PROCPOS);
        end

         for j = 1:SND_CNT
               if recmap(PROCID(j)+1,PROCPOS(j))<0
                    recmap(PROCID(j)+1,PROCPOS(j)) = NID;
                    sndmap(NID+1,j) = PROCID(j);
               else
                    recmap(PROCID(j)+1,PROCPOS(j)) = -NID;
                    sndmap(NID+1,j) = -PROCID(j);
               end 
         end

    end
                
end

fclose(fid);

recmap_trim = recmap(:,1:rec_max);
recmap_trim = [transpose(0:(s-1)) recmap_trim];


sndmap_trim = sndmap(:,1:snd_max);
sndmap_trim = [transpose(0:(s-1)) sndmap_trim];
recnos = [transpose(0:(s-1)) recnos];







