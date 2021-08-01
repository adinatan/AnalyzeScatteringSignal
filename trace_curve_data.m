function out = trace_curve_data(t)

[yy, xx]=find(t);
 
% sort x,y based on proximity to the curve, themselves....
% clear data pd N result v
data = [xx  , yy ];
pd = pdist2(data,data);

N = size(data,1);
result = NaN(1,N);
result(1) = 1; % first point is first row in data matrix

     
for ii=2:N
    pd(:,result(ii-1)) = Inf;
    [~, closest_idx] = min(pd(result(ii-1),:));
    result(ii) = closest_idx;
   
%    plot( xx(result(ii)), yy(result(ii)),'ro')
    v(ii)=G0(yy(result(ii)),  xx(result(ii)),n );     
    
end

out.v=v;
out.res=result;
