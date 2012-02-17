function meshplot(X, Y,  varargin)


nx = size(X, 2);
ny = size(X, 1);

% newplot;`

% axis([gridx(1),gridx(end)],[gridy(1),gridy(end)]);
for i = 1:nx
    hold on;
    line( X(:,i),  Y(:,i),  varargin{:});
    hold off;
end

for i = 1:ny
    hold on;
    line( X(i,:), Y(i,:),  varargin{:});
    hold off;
end