% It is a code to calculate velocity and position from accelerometer
% acc_data: accelerometer data
% sf: sampling frequency ex) sf=600.615;
% Copyright (c) 2023 Hong Gi Yeom

function [velocity, position]=acc2pos(acc_data,sf)

d_ord=10; % detrend oder
eve=length(acc_data);
for i=1:size(acc_data,2)
    events_size(i)=size(acc_data{i},3);
end
%--------------------initialization-------------------------
for i=1:eve
    for j=1:events_size(i)
        for m=1:3
            fil_acc_data{i}(m,:,j)=detrend(acc_data{i}(m,:,j),'linear',d_ord);
        end
    end
end

%---------------calculation of velocity--------------------
for i=1:eve
    for j=1:events_size(i)
        v{i}(1,:,j)=Integral(fil_acc_data{i}(1,:,j),1/sf);
        v{i}(2,:,j)=Integral(fil_acc_data{i}(2,:,j),1/sf);
        v{i}(3,:,j)=Integral(fil_acc_data{i}(3,:,j),1/sf);
    end
end

%---------------detrend velocity---------------------------
for i=1:eve
    for j=1:events_size(i)
        v{i}(1,:,j)=detrend(v{i}(1,:,j),'linear',d_ord); 
        v{i}(2,:,j)=detrend(v{i}(2,:,j),'linear',d_ord); 
        v{i}(3,:,j)=detrend(v{i}(3,:,j),'linear',d_ord);
    end
end

%--------------calculation of position---------------------
for i=1:eve
    for j=1:events_size(i)
        p{i}(1,:,j)=Integral(v{i}(1,:,j),1/sf);
        p{i}(2,:,j)=Integral(v{i}(2,:,j),1/sf);
        p{i}(3,:,j)=Integral(v{i}(3,:,j),1/sf);        
    end
end

%--------------detrend position----------------------------
for i=1:eve
    for j=1:events_size(i)
        p{i}(1,:,j)=detrend(p{i}(1,:,j),'linear',d_ord); p{i}(2,:,j)=detrend(p{i}(2,:,j),'linear',d_ord); p{i}(3,:,j)=detrend(p{i}(3,:,j),'linear',d_ord);
    end
end

%--------------setting initial point to zero----------------------------
for i=1:eve
    for j=1:events_size(i)
        v{i}(1,:,j)=v{i}(1,:,j)-v{i}(1,1,j); v{i}(2,:,j)=v{i}(2,:,j)-v{i}(2,1,j); v{i}(3,:,j)=v{i}(3,:,j)-v{i}(3,1,j); 
        p{i}(1,:,j)=p{i}(1,:,j)-p{i}(1,1,j); p{i}(2,:,j)=p{i}(2,:,j)-p{i}(2,1,j); p{i}(3,:,j)=p{i}(3,:,j)-p{i}(3,1,j); 
    end
end

velocity=v; position=p;