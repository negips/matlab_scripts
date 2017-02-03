function umtaps_merge(espfile, tsvfile, h5file)
% umtaps_merge(espfile, tsvfile, h5file)
% 
% Merge experimental pressure and position data and save to single HDF5
% file.


    % delay from trigger signal to first data point
    espdelay = 0.0; % 87.5e-3;

    if ~isempty(espfile)
        
        fid = fopen(espfile, 'r');
        [n, ~] = fread(fid, 1, 'int32', 0, 'ieee-be');
        [a, ~] = fread(fid, [n, inf], 'double', 0, 'ieee-be');
        fclose(fid);
        
        a = a';
        [nt, ncol] = size(a);
        h5create(h5file, '/esp/time', [nt,1], ...
                 'ChunkSize', [nt,1], 'Deflate', 1);
             
        % time of trigger event is t = 0
        h5write(h5file, '/esp/time', a(:,2)*2e-3 + espdelay );
        h5create(h5file, '/esp/pressure', [nt,ncol-2],  ...
                 'ChunkSize', [nt,1], 'Deflate', 1);
        h5write(h5file, '/esp/pressure', a(:, 3:end));
        
        taps = [   2.7651230e-03   8.8610830e-03
                   5.9005880e-03   1.3450120e-02
                   1.1148280e-02   1.9036800e-02
                   2.0215680e-02   2.6200210e-02
                   3.6381630e-02   3.5717800e-02
                   6.1963950e-02   4.7028130e-02
                   9.3964250e-02   5.7858090e-02
                   1.2901250e-01   6.7299890e-02
                   1.6551640e-01   7.5313560e-02
                   2.0284090e-01   8.2004990e-02
                   2.4073510e-01   8.7499130e-02
                   2.7903930e-01   9.1901840e-02
                   3.1763440e-01   9.5296990e-02
                   3.5643010e-01   9.7742050e-02
                   3.9536040e-01   9.9281890e-02
                   4.3434210e-01   9.9939280e-02
                   4.7330680e-01   9.9713170e-02
                   5.1216120e-01   9.8581870e-02
                   5.5081010e-01   9.6487950e-02
                   5.8917820e-01   9.3352220e-02
                   6.2718390e-01   8.9083890e-02
                   6.6470510e-01   8.3582460e-02
                   7.0148630e-01   7.6737500e-02
                   7.3714000e-01   6.8334350e-02
                   7.7221100e-01   5.7949820e-02
                   8.0790660e-01   4.5474600e-02
                   8.4430000e-01   3.1988270e-02
                   8.7875770e-01   2.0599030e-02
                   9.1301630e-01   1.1827610e-02
                   9.4763080e-01   5.3954460e-03
                   9.7950000e-01   1.6198060e-03
                   9.4471790e-01  -9.6501950e-03
                   9.0777710e-01  -1.3673110e-02
                   8.6932840e-01  -1.6508920e-02
                   8.2999530e-01  -1.8550180e-02
                   7.9003460e-01  -2.0094400e-02
                   7.4967370e-01  -2.1368990e-02
                   7.0910890e-01  -2.2517780e-02
                   6.6852570e-01  -2.3577410e-02
                   6.2795070e-01  -2.4537990e-02
                   5.8736170e-01  -2.5398090e-02
                   5.4672300e-01  -2.6174260e-02
                   5.0602170e-01  -2.6897440e-02
                   4.6527660e-01  -2.7599900e-02
                   4.2456280e-01  -2.8316080e-02
                   3.8398790e-01  -2.9007740e-02
                   3.4351900e-01  -2.9535160e-02
                   3.0299100e-01  -2.9939340e-02
                   2.6266500e-01  -3.0171730e-02
                   2.2247640e-01  -3.0137230e-02
                   1.8262500e-01  -2.9733840e-02
                   1.4331680e-01  -2.8781700e-02
                   1.0469320e-01  -2.6989690e-02
                   6.8153420e-02  -2.4205460e-02
                   3.9375470e-02  -2.0669700e-02
                   2.2635530e-02  -1.7421540e-02
                   1.3150860e-02  -1.4479590e-02
                   7.1869100e-03  -1.1479570e-02
                   3.4316200e-03  -8.1609190e-03
                   1.2664350e-03  -4.7221460e-03
                   1.8223130e-04  -1.4792530e-03
                   1.1376170e-04   1.5617770e-03  ];
               
        h5create(h5file, '/esp/taps', size(taps));
        h5write(h5file, '/esp/taps', taps);
        h5writeatt(h5file, '/esp', 'filename', espfile); 
        
        clear a;
    end
    
    if ~isempty(tsvfile)
        
        [pos, t, events] = import_tsv(tsvfile);
        
        evtime = zeros(numel(events),1);
        for ki = 1:numel(events)
            evtime(ki) = events(ki).time;
        end
        
        % shift time so that trigger event is at t = 0
        if numel(evtime) > 0
            t = t - evtime(1);
        end
        
        evtime = evtime - evtime(1);
        h5create(h5file, '/qsys/evtime', size(evtime));
        h5write(h5file, '/qsys/evtime', evtime);
        
        for ki = 1:numel(events)
            h5writeatt(h5file, '/qsys', ['event_name' num2str(ki)], ...
                       events(ki).name);
        end
        
        [nt, ncol] = size(pos);
        h5create(h5file, '/qsys/time', [nt,1], ...
                 'ChunkSize', [nt,1], 'Deflate', 1);
        h5write(h5file, '/qsys/time', t);
        h5create(h5file, '/qsys/pos', [nt, ncol],  ...
                 'ChunkSize', [nt, 1], 'Deflate', 1);
        h5write(h5file, '/qsys/pos', pos);
        h5writeatt(h5file, '/qsys', 'filename', tsvfile);
       
    end

end