% args:
%   x: the seq to be filtered;
%   fs: sampling freq of x
%   fn: number of freq pts between [fmin, fmax]
%   [fmin, fmax]: min and max freq 
% returns:
%   spec: the spetral, in power; to db: 10*log10(spec) or pow2db(spec);
%   fsamples: the freq pts.
function [spec, fsamples] = get_spectral(x, fs, fn, fmin, fmax)
if(~exist('fmin','var'))
    fmin = 1;
end
if(~exist('fmax','var'))
    fmax = 40;
end
fvec = (linspace(fmin,fmax, fn)).'; % freq range in Hz
%hbp=firls(100,[0 0.2*fmin 0.9*fmin fmax-2 fmax+5 100]*2/fs, [0 0 1 1 0 0]); % for detrending, a bandpass filter
hbp=firls(10,[0 0.2*fmin 0.9*fmin fmax-2 fmax+5 200]*2/fs, [0 0 1 1 0 0]); % you may use other filters
q = double(x);
q = filter(hbp,1,q);
[q, fsamples] = pmtm(q, 3,fvec, fs); % fsamples is actual fequencies in Hz returned by pburg

lpf = [1 2 5 2 1]; 
lpf = lpf/sum(lpf(:));
spec = conv(q, lpf, 'same');

end