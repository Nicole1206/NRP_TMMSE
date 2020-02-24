function nummarkers(h,hnum)
% NUMMARKERS takes a vector of line handles in h
% and reduces the number of plot markers on the lines
% to num. This is useful for closely sampled data.
%
% example:
% t = 0:0.01:pi;
% p = plot(t,sin(t),'-*',t,cos(t),'r-o');
% nummarkers(p,10);
% legend('sin(t)','cos(t)')
%

% Created: Magnus Sundberg Feb 08, 2001
% Modified: [email]felonwan@gmail.com[/email] May 25, 2013
lh=length(h);
lhn=length(hnum);
if (lhn~=lh && lhn~=1)
    error('The number of markers should be a scalar or a vector with equal length to the number of lines!')
end
for n = 1:lh
    if lhn==1
        num=hnum;
    elseif lhn==lh
        num=hnum(n);
    end
    if strcmp(get(h(n),'type'),'line')
        axes(get(h(n),'parent'));
        x = get(h(n),'xdata');
        y = get(h(n),'ydata');
        lx=length(x);
    elseif(lx<2*num)
        disp('Warning: Data points are not so many. Not necessary to use this function!')
    end
    xd=diff(get(gca,'xlim'));
    yd=diff(get(gca,'ylim'));
    s = [0 cumsum(sqrt(diff(x/xd).^2+diff(y/yd).^2))];
    si = (0:num-1)*s(end)/(num-1);
    ti=zeros(num,1);
    ti(1)=1;
    ti(num)=lx;
    for i=2:num-1
        [~,ti(i)]=min(abs(s-si(i)));
    end
    xi = x(ti);
    yi = y(ti);
    linewidth = get(h(n),'linewidth');
    marker = get(h(n),'marker');
    markersize = get(h(n),'markersize');
    color = get(h(n),'color');
    style = get(h(n),'linestyle');
    % make a line with just the markers
    set(line(xi,yi),'marker',marker,'markersize',markersize,'linestyle','none','color',color,'linewidth',linewidth);
    % make a copy of the old line with no markers
    set(line(x,y),'marker','none','linestyle',style,'color',color,'linewidth',linewidth);
    % set the x- and ydata of the old line to [], this tricks legend to
    %  keep on working
    set(h(n),'xdata',[],'ydata',[]);
end
end