%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function [ci_min, ci_max] = support_get_confidence_limit(y, ci)

N = length(y);
nMin = min(y);
nMax = max(y);

if nMin == nMax
  ci_min = nMin;
  ci_max = nMax;
  return
end

B = 1000; % 10000 max precision 99.99% (1000 for 99.9%)
b = (nMin:((nMax - nMin) / (B - 1)):nMax)';
p = zeros(B, 1);
for n = 1:B 
  p(n) = sum(y <= b(n)) / N;
  if n == B
    p(n) = sum(y < b(n)) / N;
  end
end
q = nMin + p * (nMax - nMin);

ci_min_i = find(p <= (1 - ci));
ci_max_i = find(p >= ci);

if isempty(ci_min_i) || isempty(ci_max_i)
  ci_min = nMin;
  ci_max = nMax;
else
  ci_min = q(ci_min_i(end));
  ci_max = q(ci_max_i(1));
end

end % end

%-------------------------------------------------------------------------------