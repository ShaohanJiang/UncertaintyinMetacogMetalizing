function res = splitroibeta(file, tasknsub)
% if ~isfile(file) % matlab 2018
%     error('Please input a txt file!')
% end

if nargin<2
    tasknsub = [28 30 30 28 28];
end

if ~ischar(file)
    tmpbeta = mean(file,2);
else
    tmpbeta = load(file);
end

% nsub = length(tmpbeta);

res{1} = tmpbeta(1 : tasknsub(1));
res{2} = tmpbeta(tasknsub(1)+1 : tasknsub(1)+tasknsub(2));
res{3} = tmpbeta(tasknsub(1)+tasknsub(2)+1 : tasknsub(1)+tasknsub(2)+tasknsub(3));
res{4} = tmpbeta(tasknsub(1)+tasknsub(2)+tasknsub(3)+1 : tasknsub(1)+tasknsub(2)+tasknsub(3)+tasknsub(4));
res{5} = tmpbeta(tasknsub(1)+tasknsub(2)+tasknsub(3)+tasknsub(4)+1 : end);

    
end
