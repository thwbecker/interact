function r = roots(c)
% local replacement for roots() replicating MATLAB's ordering:
% roots at zero (from trailing zero coefficients) come FIRST,
% followed by eigenvalues of the companion matrix.
% Octave's own roots() puts zero roots last, which breaks codes
% that index the output assuming MATLAB's convention.
c = c(:).';
inz = find(c ~= 0);
if isempty(inz)
  r = zeros(0,1);
  return
end
nlead = inz(1) - 1;      %#ok<NASGU> leading zeros discarded
ntrail = length(c) - inz(end);
c = c(inz(1):inz(end));
r = zeros(ntrail,1);     % roots at zero, first
n = length(c);
if n > 1
  a = diag(ones(1,n-2),-1);
  a(1,:) = -c(2:n)./c(1);
  r = [r; eig(a)];
end
