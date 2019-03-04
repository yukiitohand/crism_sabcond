function [logt_new] = dragpoints_logt(wv,logt,varargin)

logt_ref = [];
ref_axis = 'left';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for n=1:2:(length(varargin)-1)
        switch upper(varargin{n}) 
            case 'REF'
                logt_ref = varargin{n+1};
            case 'REF_AXIS'
                ref_axis = varargin{n+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{n} '''']);
        end
    end
end



fig = figure;
x = wv;
y = logt;

% ax = axes('xlimmode','manual','ylimmode','manual');
% ax.XLim = [xLower xUpper];
% ax.YLim = [yLower yUpper];

%can change the marker size or marker type to make it more visible.
%Currently is set to small points at a size of 2 so is not very visible.
clrs = lines(7);
l = line(x,y,'marker','.','markersize',6,'LineStyle','-','Color',clrs(1,:));
if ~isempty(logt_ref)
    yyaxis(ref_axis);
    lref = line(x,logt_ref,'marker','.','markersize',6,'LineStyle','-',...
                'Color',clrs(5,:));
    yyaxis left;
end
lnew = line(x,y,'marker','.','markersize',6,'LineStyle','-','Color',clrs(2,:),...
            'hittest','on','buttondownfcn',@clickmarker);

b = uicontrol(fig,'Style','pushbutton','String','OK',...
    'Units','normalized','Position', [0.02 0.9 0.1 0.05],...
    'Callback','uiresume(gcbf)');

uiwait(fig);

logt_new = lnew.YData;

close(fig);

end

function fin_update(src,ev)
    set(gcbo,'UserData',1);
end

function clickmarker(src,ev)
%get current axes and coords
h1=gca;
coords=get(h1,'currentpoint');

%get all x and y data 
x=src.XData;
%y=h1.Children.YData;

%check which data point has the smallest distance to the dragged point
x_diff=abs(x-coords(1,1,1));

[value,index]=min(x_diff);

set(ancestor(src,'figure'),'windowbuttonmotionfcn',{@dragmarker,src,index})
set(ancestor(src,'figure'),'windowbuttonupfcn',@stopdragging)

end

function dragmarker(fig,ev,src,index)

%get current axes and coords
h1=gca;
coords=get(h1,'currentpoint');

%get all x and y data 
x=src.XData;
y=src.YData;

%check which data point has the smallest distance to the dragged point
%x_diff=abs(x-coords(1,1,1));
%y_diff = 0;
%y_diff=abs(y-coords(1,2,1));
%[value index]=min(x_diff+y_diff);

%create new x and y data and exchange coords for the dragged point
%x_new=x;
y_new=y;
%x_new(index)=coords(1,1,1);
y_new(index)=coords(1,2,1);

%update plot
set(src,'xdata',x,'ydata',y_new);

end

function stopdragging(fig,ev)
set(fig,'windowbuttonmotionfcn','')
set(fig,'windowbuttonupfcn','')

end


