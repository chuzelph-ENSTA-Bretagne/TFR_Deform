# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 22:23:48 2015

@author: Philippe
"""

'''
function [tfr,t,f] = tfrspwv(x,t,N,g,h,trace);
%TFRSPWV Smoothed Pseudo Wigner-Ville time-frequency distribution.
%[TFR,T,F]=TFRSPWV(X,T,N,G,H,TRACE) computes the Smoothed Pseudo
%Wigner-Ville distribution of a discrete-time signal X, or the
%cross Smoothed Pseudo Wigner-Ville representation between two
%signals.
%
%X     : signal if auto-SPWV, or [X1,X2] if cross-SPWV.
%T     : time instant(s)          (default : 1:length(X)).
%N     : number of frequency bins (default : length(X)).
%G     : time smoothing window, G(0) being forced to 1.
%                                 (default : Hamming(N/10)).
%H     : frequency smoothing window in the time-domain,
%H(0) being forced to 1   (default : Hamming(N/4)).
%TRACE : if nonzero, the progression of the algorithm is shown
%                                 (default : 0).
%TFR   : time-frequency representation. When called without
%        output arguments, TFRSPWV runs TFRQVIEW.
%F     : vector of normalized frequencies.
%
%Example :
% sig=fmlin(128,0.05,0.15)+fmlin(128,0.3,0.4);  
% g=window(15,'Kaiser'); h=window(63,'Kaiser');
% tfrspwv(sig,1:128,64,g,h,1);
%
%See also all the time-frequency representations listed in
% the file CONTENTS (TFR*)
 
%F. Auger, May-August 1994, July 1995.
%Copyright (c) 1996 by CNRS (France).
%
%------------------- CONFIDENTIAL PROGRAM --------------------
%This program can not be used without the authorization of its
%author(s). For any comment or bug report, please send e-mail to
%f.auger@ieee.org
 
if (nargin == 0),
 error('At least 1 parameter required');
end;
[xrow,xcol] = size(x);
if (xcol==0)|(xcol>2),
 error('X must have one or two columns');
end
 
if (nargin <= 2),
 N=xrow;
elseif (N<0),
 error('N must be greater than zero');
elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');
end;
 
hlength=floor(N/4); hlength=hlength+1-rem(hlength,2);
glength=floor(N/10);glength=glength+1-rem(glength,2);
 
if (nargin == 1),
 t=1:xrow; g = window(glength); h = window(hlength); trace = 0;
elseif (nargin == 2)|(nargin == 3),
 g = window(glength); h = window(hlength); trace = 0;
elseif (nargin == 4),
 h = window(hlength); trace = 0;
elseif (nargin == 5),
 trace = 0;
end;
 
[trow,tcol] = size(t);
if (trow~=1),
 error('T must only have one row');
end;
 
[grow,gcol]=size(g); Lg=(grow-1)/2; % g=g/sum(g);
if (gcol~=1)|(rem(grow,2)==0),
 error('G must be a smoothing window with odd length');
end;
 
[hrow,hcol]=size(h); Lh=(hrow-1)/2; h=h/h(Lh+1);
if (hcol~=1)|(rem(hrow,2)==0),
  error('H must be a smoothing window with odd length');
end;
 
tfr= zeros (N,tcol) ; 
if trace, disp('Smoothed pseudo Wigner-Ville distribution'); end;
for icol=1:tcol,
 ti= t(icol); taumax=min([ti+Lg-1,xrow-ti+Lg,round(N/2)-1,Lh]);
 if trace, disprog(icol,tcol,10); end;
 points= -min([Lg,xrow-ti]):min([Lg,ti-1]);
 g2=g(Lg+1+points); g2=g2/sum(g2);
 tfr(1,icol)= sum(g2 .* x(ti-points,1) .* conj(x(ti-points,xcol)));
 for tau=1:taumax,
  points= -min([Lg,xrow-ti-tau]):min([Lg,ti-tau-1]);
  g2=g(Lg+1+points); g2=g2/sum(g2);
  R=sum(g2 .* x(ti+tau-points,1) .* conj(x(ti-tau-points,xcol)));
  tfr(  1+tau,icol)=h(Lh+tau+1)*R;
  R=sum(g2 .* x(ti-tau-points,1) .* conj(x(ti+tau-points,xcol)));
  tfr(N+1-tau,icol)=h(Lh-tau+1)*R;
 end;
 tau=round(N/2);
 if (ti<=xrow-tau)&(ti>=tau+1)&(tau<=Lh),
  points= -min([Lg,xrow-ti-tau]):min([Lg,ti-tau-1]);
  g2=g(Lg+1+points); g2=g2/sum(g2);
  tfr(tau+1,icol) = 0.5 * ...
   (h(Lh+tau+1)*sum(g2 .* x(ti+tau-points,1) .* conj(x(ti-tau-points,xcol)))+...
    h(Lh-tau+1)*sum(g2 .* x(ti-tau-points,1) .* conj(x(ti+tau-points,xcol))));
 end;
end;
if trace, fprintf('\n'); end;
tfr= fft(tfr);
if (xcol==1), tfr=real(tfr); end ;
 
if (nargout==0),
 tfrqview(tfr,x,t,'tfrspwv',g,h);
elseif (nargout==3),
 f=(0.5*(0:N-1)/N)';
end;

'''