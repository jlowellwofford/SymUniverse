function MOV = mkmov( filename )
%MKMOV Create a movie from a universe csv file
%   input: filename
%   output: MOV object

M = dlmread(filename,',',1,0);
[ts,ia,ic] = unique(M(:,1));
clear ic;
ia = [ia; length(M(:,1))+1];
ax = [ 
    min(M(ia(2):length(M(:,1)),13)) max(M(ia(2):length(M(:,1)),13)) ...
    min(M(ia(2):length(M(:,1)),14)) max(M(ia(2):length(M(:,1)),14)) ...
    min(M(ia(2):length(M(:,1)),15)) max(M(ia(2):length(M(:,1)),15)) ];
for i = 1:length(ts)
    N = (ia(i+1) - ia(i));
    scatter3(...
        M(ia(i):ia(i+1)-1,13),...
        M(ia(i):ia(i+1)-1,14),...
        M(ia(i):ia(i+1)-1,15),4, [transpose(1:N)/N zeros(N,2) ], 'filled');
    axis manual;
    axis(ax);
    MOV(i) = getframe(gcf);
end

end

