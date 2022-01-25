function [header out] = parse_csv(data);
lines = regexp(data,'\n','split');
[a,n]  = size(lines);
header = regexp(lines(1),',','split');
[b,v]  = size(header);
if isempty(lines{n})
	n = n-1;
end
out = cell(n-1,v);
for i=2:n
	ln = regexp(lines(i),',','split');
	for j=1:v
		out(i-1,j) = ln(j);
	end
end
