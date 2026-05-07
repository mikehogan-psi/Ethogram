function[] = chessboard_video()
device = serialport("COM4",9600);
T = 0.5; %for experimentL 0.05
Nimg = 300;
pause(4*T);
for n = 1:Nimg
    write(device,'t',"char");
    pause(T);
end
clear device;