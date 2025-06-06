%  Function to calculate the dft of an input sequence
%
%  Usage is f = dft(x,w)
%
%  where
%
%  x = input sequence
%  w = vector of normalised frequencies to calculate dft at
%  f = dft
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function f = dft(x,w);

ww = exp(-j*w); ww = ww(:).'; x = x(:).';
f = polyval(fliplr(x),ww);