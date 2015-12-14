%{
vangogh.RF (computed) # RF mapping for van Gogh noise
-> rf.Sync
-> vangogh.RFMethod
-> rf.Trace
-----
nbins              : smallint                      # temporal bins
bin_size           : float                         # (ms) temporal bin size
degrees_x          : float                         # degrees along x
degrees_y          : float                         # degrees along y
map                : longblob                      # receptive field map
%}

classdef RF < dj.Relvar & dj.AutoPopulate

	properties
		popRel = rf.Sync * rf.Segment * vangogh.RFMethod & (psy.Session * psy.Trial * psy.VanGogh)
	end

	methods (Access = protected)

		function makeTuples(self, key)
            
            % temporal binning
            nBins = 8;
            binSize = 0.1;  %s
            
            disp 'loading movies ...'
            caTimes = fetch1(rf.Sync & key,'frame_times');
            dt = median(diff(caTimes));
            trials = rf.Sync * psy.Trial * psy.VanGogh * pro(psy.VanGoghLookup, 'movie');
            trials = trials & key & 'trial_idx between first_trial and last_trial';
            [stimTimes, movie] = trials.fetchn('flip_times', 'movie', 'ORDER BY trial_idx');
            % compute physical dimensions
            sess = fetch(rf.Sync * psy.Session & key, 'resolution_x', 'resolution_y', 'monitor_distance', 'monitor_size');
            rect = [sess.resolution_x sess.resolution_y];
            degPerPix = 180 / pi * sess.monitor_size * 2.54 / norm(rect(1 : 2)) / sess.monitor_distance;
            degSize = degPerPix * rect;
            
            disp 'concatenating stimulus movies...'
            stimTimes = cat(2, stimTimes{:});
            movie = double(cat(3, movie{:})) / 127 - 1;
            
            disp 'interpolation...'
            % clip stimulus movie to fit within the calcium recording to avoid extrapolation
            ix = stimTimes > caTimes(1) & stimTimes < caTimes(end) - nBins * binSize;
            movie = movie(:, :, ix);
            stimTimes = stimTimes(ix);
            
            t0 = max(caTimes(1), stimTimes(1) + (nBins - 1) * binSize) + 0.1;   % start time for calcium traces
            t1 = t0 - (nBins - 1) * binSize;                                 % start time for stimulus
            t2 = min(caTimes(end), stimTimes(end));                     % end time for both
            movie = permute(interp1(stimTimes', permute(movie,[3 1 2]), (t1 : binSize : t2)', 'linear'), [2 3 1]);
            
            % normalize movie
            movie = bsxfun(@minus, movie, mean(movie, 3));
            movie = bsxfun(@rdivide, movie, std(movie, [], 3));
            sz = size(movie);

            method = fetch1(vangogh.RFMethod & key, 'method');
            
            disp 'computing RF...'
            [traces, traceKey] = fetchn(rf.Trace & key, 'ca_trace');
            for iTrace = 1 : length(traces)
                fprintf('trace %d\n', traceKey(iTrace).trace_id)
                tuple = dj.struct.join(key, traceKey(iTrace));
                
                % highpass filter and deconvolve
                cutoff = 0.03;
                k = hamming(round(1 / dt / cutoff) * 2 + 1);
                k = k / sum(k);
                trace = double(traces{iTrace});
                trace = (trace - ne7.dsp.convmirr(double(trace), k)) / mean(trace);
                trace = fast_oopsi(trace, struct('dt', dt), struct('lambda', 0.3));
                
                % interpolate to common time bins
                k = hamming(2 * round(binSize / dt) + 1);
                k = k / sum(k);
                trace = ne7.dsp.convmirr(trace, k);
                trace = interp1(caTimes, trace ,(t0 : binSize : t2)', 'linear');
                trace = trace / sum(trace);
                
                disp 'computing RF...'
                switch method
                    case 'STA'                        
                        map = reshape(conv2(fliplr(reshape(movie, sz(1) * sz(2), sz(3))), trace', 'valid'), sz(1), sz(2), []);
                    otherwise
                        error('The "%s" method is not implemented yet', method)
                end
                
                disp 'saving..'
                
                tuple.nbins = nBins;
                tuple.bin_size = binSize * 1000;
                tuple.degrees_x = degSize(1);
                tuple.degrees_y = degSize(2);
                tuple.map = single(map);
                
                imagesc(map(:, :, 2), [-1 1] * 0.05), axis image
                colormap(ne7.vis.doppler)
                drawnow
                self.insert(tuple)
            end
            disp done
		end
	end

end