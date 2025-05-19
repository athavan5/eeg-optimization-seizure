tt = edfread('chb01_01.edf');
info = edfinfo('chb01_01.edf');

fps = 256;

second = 3600; % enter second (usually 1-3600)
channel = 23; % enter channel (1-23)
t = (0:fps-1)/fps;
y = tt.(channel){second};

plot(t,y)
legend(strcat("Second of File: ",int2str(second), ...
    ", Channel: ",info.SignalLabels(channel)))