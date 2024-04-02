%%%%%%%%%
% This is "tfrscalo.m", tabulations, variable names, etc. slightly modified by 
% St\'ephane Rossignol in order to make the program more readable
%%%%%%%%%

function [tfr,t,f,wt] = lagis_mytfrscalo(X, time, wave, fmin, fmax, NNN, trace);
%TFRSCALO Scalogram, for Morlet or Mexican hat wavelet.
%	[TFR,T,F,WT]=TFRSCALO(X,T,WAVE,FMIN,FMAX,NNN,TRACE) computes 
%	the scalogram (squared magnitude of a continuous wavelet
%	transform). 
%
%	X : signal (in time) to be analyzed (Nx=length(X)). Its
%	    analytic version is used (zzz=hilbert(real(X))).  
%	T : time instant(s) on which the TFR is evaluated 
%	     					(default : 1:Nx).
%	WAVE : half length of the Morlet analyzing wavelet at coarsest 
% 	    scale. If WAVE = 0, the Mexican hat is used. WAVE can also be
%           a vector containing the time samples of any bandpass
%           function, at any scale.        	(default : sqrt(Nx)). 
%	FMIN,FMAX : respectively lower and upper frequency bounds of 
%	    the analyzed signal. These parameters fix the equivalent
%	    frequency bandwidth (expressed in Hz). When unspecified, you
%	    have to enter them at the command line from the plot of the
%	    spectrum. FMIN and FMAX must be >0 and <=0.5.
%	NNN : number of analyzed voices.
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 	(default : 0).
%	TFR : time-frequency matrix containing the coefficients of the
%	    decomposition (abscissa correspond to uniformly sampled time,
%	    and ordinates correspond to a geometrically sampled
%	    frequency). First row of TFR corresponds to the lowest 
%	    frequency. When called without output arguments, TFRSCALO
%	    runs TFRQVIEW.
%	F : vector of normalized frequencies (geometrically sampled 
%	    from FMIN to FMAX).
%	WT : Complex matrix containing the corresponding wavelet
%	    transform. The scalogram TFR is the square modulus of WT.
%
%	Example :    
%	 sig=altes(64,0.1,0.45); tfrscalo(sig);  
%
%	See also all the time-frequency representations listed in
%	the file CONTENTS (TFR*)

%	P. Goncalves, October 1995 - O. Lemoine, June 1996. 
%	Copyright (c) 1995 Rice University - CNRS 1996.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

if (nargin == 0),
  error('At least one parameter required');
end;

[xrow,xcol] = size(X);
if nargin<=6, trace=0; end

if (nargin == 1),
  time=1:xrow; wave=sqrt(xrow);
elseif (nargin == 2),
  wave=sqrt(xrow);
elseif (nargin==4),
  disp('FMIN will not be taken into account. Determine it with FMAX');
  disp('     from the following plot of the spectrum.'); 
elseif nargin==5,
  NNN=xrow;
end;

[trow,tcol] = size(time);
if (xcol==0)|(xcol>2),
  error('X must have one or two columns');
elseif (trow~=1),
  error('TIME must only have one row'); 
elseif wave<0,
  error('WAVE must be positive');
end; 

sss = (real(X) - mean(real(X)))';  
zzz = hilbert(sss);

if trace, disp('Scalogram distribution'); end;

if nargin<=4		        % fmin,fmax,NNN unspecified
 STF = fft(fftshift(zzz(min(time):max(time)))); Nstf=length(STF);
 sp = (abs(STF(1:round(Nstf/2)))).^2; Maxsp=max(sp);
 fff = linspace(0,0.5,round(Nstf/2)+1);
 fff = fff(1:round(Nstf/2));
 plot(f,sp);
 grid;
 xlabel('Normalized frequency');
 title('Analyzed signal energy spectrum');
 indmin=min(find(sp>Maxsp/100));
 indmax=max(find(sp>Maxsp/100));
 fmindflt=max([0.01 0.05*fix(f(indmin)/0.05)]);
 fmaxdflt=0.05*ceil(f(indmax)/0.05);
 txtmin=['Lower frequency bound [',num2str(fmindflt),'] : '];
 txtmax=['Upper frequency bound [',num2str(fmaxdflt),'] : '];
 fmin = input(txtmin); fmax = input(txtmax);
 if isempty(fmin), fmin=fmindflt; end
 if isempty(fmax), fmax=fmaxdflt; end
 txt=['Number of frequency samples [',num2str(2^nextpow2(xrow)),'] : ']; 
 NNN=input(txt); 
 if isempty(NNN), NNN=2^nextpow2(xrow); end
end

fmin_s=num2str(fmin); fmax_s=num2str(fmax); 
N_s=num2str(NNN);

if fmin >= fmax
 error('FMAX must be greater or equal to FMIN');
elseif fmin<=0.0 | fmin>0.5,
 error('FMIN must be > 0 and <= 0.5');
elseif fmax<=0.0 | fmax>0.5,
 error('FMAX must be > 0 and <= 0.5');
end
if trace,
 disp(['Frequency runs from ',fmin_s,' to ',fmax_s,' with ',N_s,' points']);
end

whichspace=2;
switch whichspace
  case 1
    %%% logarithmic
    %%% => original way
    fff = logspace(log10(fmin),log10(fmax),NNN);
    aaa = logspace(log10(fmax/fmin),log10(1),NNN);

    %figure(1);
    %plot(fff);
    %figure(2);
    %plot(aaa);
    figure(3);
    plot(fff.*aaa);
    fprintf(1,"| %f %f | %f %f |\n", min(fff), max(fff), min(aaa), max(aaa));
    fflush(1);

  case 2
    %%% linear
    fff = linspace(fmin,fmax,NNN);
    aaa = max(fff)./fff;
    %aaa = linspace(fmax/fmin,1,NNN); %%% => 30/04/07; does not work: the product aaa.*fff needs to be constant

    %figure(1);
    %plot(fff);
    %figure(2);
    %plot(aaa);
    %figure(3);
    %plot(fff.*aaa);
    fprintf(1,"| %f %f | %f %f |\n", min(fff), max(fff), min(aaa), max(aaa));
    fflush(1);

end;


wt  = zeros(NNN,tcol);
tfr = zeros(NNN,tcol);

if wave > 0
 if trace, disp(['using a Morlet wavelet']) ; end
 for ptr=1:NNN
  if trace, disprog(ptr,NNN,10); end
  nha = wave*aaa(ptr);
  tha = -round(nha) : round(nha);
  ha  = exp(-(2*log(10)/nha^2)*tha.^2).*exp(i*2*pi*fff(ptr)*tha);
  %fprintf(1,'aa=%f  ff=%f\n',aaa(ptr),fff(ptr));pause
  %if ptr==NNN
  %  size(ha)
  %  pause;
  %end;
  %figure(1);plot(abs(ha));figure(2);plot(real(ha));pause;
  detail = conv(zzz,ha)./sqrt(aaa(ptr));
  detail = detail(round(nha)+1:length(detail)-round(nha)) ;
  wt(ptr,:)  = detail(time) ;
  tfr(ptr,:) = detail(time).*conj(detail(time)) ;
 end
elseif wave == 0
 if trace, disp(['using a Mexican hat wavelet']) ; end
 for ptr = 1:NNN
  if trace, disprog(ptr,NNN,10); end
  % St\'ephane Rossignol added this parameter "zoo" the 24th of April 2007
  % => but it does not work, thus put "zoo=1"
  zoo=1;
  ha  = mexhat(fff(ptr), zoo) ;
  nha = (length(ha)-1)/2 ;
  %plot(ha);pause;
  detail = conv(zzz,ha)./sqrt(aaa(ptr));
  detail = detail(round(nha)+1:length(detail)-round(nha)) ;
  wt(ptr,:)  = detail(time);
  tfr(ptr,:) = detail(time).*conj(detail(time)) ;
 end;
elseif length(wave) > 1
 [rwav,cwav]=size(wave);
 if cwav>rwav, wave=wave.'; end
 wavef = fft(wave) ;
 nwave = length(wave) ;
 f0 = find(abs(wavef(1:nwave/2)) == max(abs(wavef(1:nwave/2))));
 f0 = mean((f0-1).*(1/nwave));
 if trace, disp(['mother wavelet centered at f0 = ',num2str(f0)]); end
 aaa = logspace(log10(f0/fmin),log10(f0/fmax),NNN);
 B = 0.99;
 R = B/((1.001)/2); 
 nscale = max(128,round((B*nwave*(1+2/R)*log((1+R/2)/(1-R/2)))/2));
 if trace, disp('Scale computation :'); end
 wts = scale(wave,aaa,fmin,fmax,nscale,trace);
 for ptr = 1:NNN
  clear detail
  if trace, disprog(ptr,NNN,10); end
  ha = wts(ptr,:);
  nha = length(ha)/2;
  detail = conv(zzz,ha)./sqrt(aaa(ptr));
  detail = detail(fix(nha):length(detail)-round(nha));
  wt(ptr,:) = detail(time);
  tfr(ptr,:) = detail(time).*conj(detail(time));
 end
end


ttt = time;
fff = fff';

% Normalization
SP = fft(zzz); 
indmin = 1+round(fmin*(xrow-2));
indmax = 1+round(fmax*(xrow-2));
SPana=SP(indmin:indmax);
tfr=tfr*norm(SPana)^2/integ2d(tfr, ttt, fff)/NNN;

if (nargout==0),
 tfrqview(tfr, hilbert(real(X)), ttt, 'tfrscalo', wave, NNN, fff);
end;

