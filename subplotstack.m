function p = subplotstack(m,ROW)
%SUBPLOTSTACK
%
%   p = SUBPLOTSTACK(m,r) returns a position array that can be passed to
%   SUBPLOT. 
%
%   subplot('position',subplotstack(2,1));plot(rand(10,1)); set(gca,
%   'Xticklabel', 'off')
%   subplot('position',subplotstack(2,2));plot(rand(10,1));
%   
%   See also SUBPLOT.
  
N     = m;
dlef  = 0.07;        % Space on left
drig  = 0.07;        % Space on right
w     = 1-dlef-drig; % Width of plot area
gap   = 0.02;       % Gap between plot areas
dtop  = 0.06;        % Space on top
dbot  = 0.09;        % Space on bottom
He    = 1 - ( gap*(N-1)+dtop+dbot ); 
h     = He/N;

t = 1-dtop-h;
l = dlef;

for i = 1:ROW
  tlast = t;
  p = [l t w h];
  t = tlast - (h + gap);
end
