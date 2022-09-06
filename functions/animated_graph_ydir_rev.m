
function animated_graph(trials_x,trials_y,x_name,y_name,title_name,save_name,frame_rate,colormap)
% function plots x and y data from a number of trials one point at a time
% then save to AVI video file. File will save in current directory.
% INPUT:
%   - trials_x: x axis data. Individual trials span rows, number of rows is
%               the number of trials.
%   - trials_y: y axis data. Individual trials span rows, number of rows is
%               the number of trials.
%   - x_name:  x axis label
%   - y_name: y axis label.
%   - title: title of graph in str format
%   - save_name: name of the saved video.
%   - frame_rate: frames per second that the video will be written to. 
%   - colormap:  a cbrewer colormap to plot the lines accoring to
    
    % write frames (they will write to the current path)
    writerObj = VideoWriter(save_name,'MPEG-4');
    writerObj.FrameRate = frame_rate;
    
    % open the video writer
    open(writerObj);

    % loop through each step to plot
    for stepnum = 1:size(trials_x,2)
        
        hold on
        for trial = 1:size(trials_x,1)
            plot(trials_x(trial,1:stepnum),trials_y(trial,1:stepnum),'Color',colormap(trial,:))
        end
        
        % format axes 
        tx = reshape(trials_x,[numel(trials_x),1]);
        ty = reshape(trials_y,[numel(trials_y),1]);
        xlim([min(tx) max(tx)])
        ylim([min(ty) max(ty)])
        xlabel(x_name)
        ylabel(y_name)
        title(title_name)
        set(gca, 'YDir','reverse');
        
        % get frame
        F = getframe(gcf);
        
        % convert the image to a frame
        frame = F.cdata ;    
        
        % write to video
        writeVideo(writerObj, frame);
        
        cla
    end
    % close the writer object
    close(writerObj);
    
end
