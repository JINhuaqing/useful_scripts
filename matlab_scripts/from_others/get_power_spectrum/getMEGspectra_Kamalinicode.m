clear all;close all;clc
load ('timecourse_ucsfCONT_group.mat')
% cd /Users/kamaliniranasinghe/Dropbox/RESEARCH/SPECGRAPH/matfiles
%AD
for jj=1:length(dk10(:,1,1))
  dk1=squeeze(dk10(jj,:,:));
  MEGdata=dk1;
  fmin = 1; %Hz
  fmax = 40;
  fvec = (linspace(fmin,fmax, 40)).'; % freq range in Hz
  fs=600;
  hbp=firls(100,[0 0.2*fmin 0.9*fmin fmax-2 fmax+5 100]*2/fs, [0 0 1 1 0 0]); % for detrending, a bandpass filter
  lpf = [1 2 5 2 1]; lpf = lpf/sum(lpf(:));
  nroiMEG=68;
  figure;hold
  for i = 1:nroiMEG
    q = double(MEGdata(i,:));
    q = filter(hbp,1,q);
    [FMEGdata(i,:), fsamples] = pmtm(q, 3,fvec, fs); % fsamples is actual fequencies in Hz returned by pburg
    q = FMEGdata(i,:);
    FMEGdata(i,:) = conv(q, lpf, 'same');
    plot(FMEGdata(i,:));
    title(char(radid(jj,:)));
    fmegall(jj,i,:)=FMEGdata(i,:);
  end
  pause
  % eval([‘print -dtiff -r600 ’ char(dktimecourse.radid(jj,:)) ‘spectrum’])
  close all;
end
save fmegallcont_tinittus fmegall
