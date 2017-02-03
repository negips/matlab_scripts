function [segments] = split_segments(hfile, uoo, deltacase, dt)
% [segments] = split_segments(hfile, uoo, deltacase, dt)
% 
% Read a HDF5 file with pressure and position measurements, split out 
% segments marked by pairs of events and postprocess. 
%
% hfile:     file name for HDF5 file containing measurements
% uoo:       corresponding free-stream velocity, ususally part of the file name
% deltacase: nominal flap deflection (-4,0,8,14,24)
%            actual measured flap deflection is reported in each segment
% dt:        time delay between pressure and optical measurments [optional]
%
% segments:  struct array with one entry per time interval where frequency 
%            is kept constant. 

    if nargin > 3
        tshift = dt;
    else
        tshift = 0;
    end

    b = 0.5*0.5;
    qoo = 0.5*1.225*uoo^2;
    evt = h5read(hfile, '/qsys/evtime');
    if ~isempty(strfind(hfile, 'afsteps'))
        rpms = [60, 150:150:2400];
    elseif ~isempty(strfind(hfile, 'hfsteps'))
        rpms = [600, 1200, 1500:150:2400];
    else
        rpms = [60, 150:150:1500];
    end
    
    % first event is trigger
    nev = numel(evt);
    nseg = fix((nev - 1)/2);
    if nev ~= 2*numel(rpms) + 1
        isfsteps = 0;
    else
        isfsteps = 1;
    end
    
    [alpha, delta, tq] = marker2angles(hfile, deltacase);
    
    % split motion segments out 
    segments = struct([]);
    pos = h5read(hfile, '/qsys/pos');
    for ki = 1:nseg
        t1 = evt(2*ki);
        t2 = evt(2*ki+1);
        idx = (tq >= t1) & (tq <= t2);
        if isfsteps
            segments(ki).rpm = rpms(ki);
            segments(ki).rfreq = 2*pi*rpms(ki)/300. * b/uoo; 
        end
        segments(ki).qtime = tq(idx);
        segments(ki).pos = pos(idx, :);
        segments(ki).alpha = alpha(idx);
        segments(ki).delta = delta(idx);
    end
    
    % mask out the pressure tap covered by hinge tape
    xytaps = h5read(hfile, '/esp/taps');
    itaps = true(62,1);
    if deltacase == 14 || deltacase == 24
        itaps(27) = false;
    end
    xtaps = xytaps(itaps,1);
    ytaps = xytaps(itaps,2);
    
    % extract pressure segments
    tp = h5read(hfile, '/esp/time');
    p = h5read(hfile, '/esp/pressure');
    
    % DTC clock skew correction 
    tp = 1.002491251607135*tp;
    
    for ki = 1:nseg
        it1 = find((tp+tshift)>=(evt(2*ki)), 1, 'first');
        it2 = find((tp+tshift)<=(evt(2*ki+1)), 1, 'last');
        segments(ki).ptime = tp(it1:it2) + tshift;
        segments(ki).pressure = p(it1:it2, :);
        segments(ki).Cz = trapz(xtaps, p(it1:it2, itaps),2) / qoo;
        dxt = ones(it2-it1+1,1) * (0.25 - xtaps)';
        segments(ki).Cm = trapz(xtaps, (p(it1:it2, itaps).*dxt), 2) / qoo;
        segments(ki).xtaps = xtaps;
        segments(ki).ytaps = ytaps;
        segments(ki).itaps = itaps;
    end
    
end
