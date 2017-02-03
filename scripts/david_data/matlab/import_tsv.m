function [pos, t, events] = import_tsv(fname)
% [pos,t,events] = import_tsv(fname)
%
% Load a TSV text file written by the QualiSys Motion Tracking system.
%
% fname:  Path to TSV file
% pos:    MxN matrix of position data, one row per timestep. Each marker
%         position takes up three columns (x,y,z)
% t:      Mx1 array of equidistant time values starting at 0.0
% events: Array of events, where each event is a struct with name string,
%         frame index and time of occurence. Empty if there are no events
%         recorded in the file header.

    fid = fopen(fname, 'r');    
    nmark = 0;
    f = 0;
    pos = [];
    events = [];
    while ~feof(fid)
        s = fgetl(fid);
        if nmark == 0
            [x, n] = sscanf(s, 'NO_OF_MARKERS\t%d', 1);
            if n == 1
                nmark = x
            end
        elseif f == 0
            [x, n] = sscanf(s, 'FREQUENCY\t%d', 1);
            if n == 1
                f = x
            end
        elseif ~isempty(strfind(s, 'EVENT'))
            n = 1;
            idx = 6;
            words = {};
            while n == 1
                [w,n,~,offset] = sscanf(s(idx:end),'%s',1);
                if n >= 1 
                    words = horzcat(words, w);
                    idx = idx + offset;
                end
            end
            ev.name = strjoin(words, ' ');
            ev.frame = sscanf(words{end-1}, '%d', 1);
            ev.time = sscanf(words{end}, '%f', 1);
            events = vertcat(events, ev);
        else
            [a,n] = sscanf(s, '%f', 3*nmark+2);
            if n == 3*nmark+2
                [p,nread] = fscanf(fid, '%f');
                nr = fix(nread/(3*nmark+2))
                p = reshape(p, [3*nmark+2, nr]);
                pos = [a'; p'];
            end
        end
    end
        
    nrow = size(pos, 1);
    t = (0:(nrow-1))' / f;
    pos = pos(:, 3:end);

end
