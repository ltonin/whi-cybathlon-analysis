function C = whi_processing_sliding_covariance(data, wlenght, wshift)


    nsamples  = size(data, 1);
    nchannels = size(data, 2);
    
    wstart = 1:wshift:(nsamples-wlenght);
    wstop  = wstart + wlenght - 1;
    
    nwins = length(wstart);
    
    C = nan(nwins, nchannels, nchannels);
    for wId = 1:nwins
        whi_util_progress_display(wId, nwins, '    ');
        cstart = wstart(wId);
        cstop  = wstop(wId);
        C(wId, :, :) = cov(data(cstart:cstop, :));
    end

end
