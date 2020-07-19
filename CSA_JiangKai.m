%% 正侧视cs成像算法  蒋恺 学号：19011110166
clear;
close all
%% 数据主要参数 参看课本P157页P波段合成孔径雷达数据
c=3e8;              %光速3e8
lambda=0.4;     %波长0.4m
f_c=c/lambda;        %频率=光速/波长
T_p=10e-6;           %脉冲宽度10us
PRF=1000;            %脉冲重复频率1000Hz
B=50e6;             %带宽50MHz
gamma=B/T_p;          %调频率
F_s=64e6;            %距离向采样率64MHz
%R_s=15e3;             %平台离场景中心的距离15km
R_s=4e3;             
V=120;              %飞机速度120m/s
D_a=1.2;             %天线方位向孔径长度1.2m
W_r = 9e3;           %条带宽度9km
theta=lambda/D_a;        %波束宽度=波长/天线方位向孔径长度
T_a=theta*R_s/V;       %方位向合成孔径时间=波束宽度*航线到场景中心的距离
num_range=floor(8*T_p*F_s/2)*2;       %距离向采样点数
num_azimuth=floor(PRF*(T_a*1)/2)*2;      %方位向采样点数
t_hat=(-num_range/2:num_range/2-1)'/F_s+2*R_s/c; %距离向采样起始与结束时间，快时间
t_m=(-num_azimuth/2:num_azimuth/2-1)/PRF;%方位向采样起始与结束时间，慢时间
%% 点目标的坐标位置
x_target_position=[0];    %点目标的横坐标，方位向
y_target_position=[0]+R_s;%点目标的纵坐标，距离向
figure;
% subplot(221);
plot(x_target_position,y_target_position,'rp');
axis([-30 30 3000 5000]);
title('目标位置');
xlabel('方位向');ylabel('距离向');
target_num=length(x_target_position);%目标数量
coordiante_azimuth=t_m*V;%方位向坐标=距离向坐标*飞机速度
signal_of_echo=zeros(num_range,num_azimuth);%用于存储回波信号
y_target=y_target_position;
%% 产生回波信号 目标回波的基频信号在距离快时间—方位慢时间域的公式 参考课本P132页公式(5.9)
for k=1:num_azimuth
    %以目标所在的位置作为横坐标原点
    x_target=x_target_position-coordiante_azimuth(k);
    for point_target_index=1:target_num
        %判断目标是否在当前孔径内，加方位窗，方位窗函数a_a(t_m)
        if (atan(x_target(point_target_index)/y_target(point_target_index))<theta/2)&&(atan(x_target(point_target_index)/y_target(point_target_index))>-theta/2)
            %任意时刻雷达天线相位中心至目标的斜距为R=R(t_m;;R_B)
            R=sqrt(x_target(point_target_index)^2+y_target(point_target_index)^2); 
            %加距离窗，距离窗函数a_r(t_hat-2*R/c)
            t0=2*R/c;
            a_r=((t_hat-t0)<T_p/2)&((t_hat-t0)>-T_p/2);
            %目标的回波的基频信号
            %课本对应公式：s(t_hat,t_m;R_B)=a_r*a_a*exp(j*pi*gamma*(t_hat-2*R/c)^2)*exp(-j*4*pi*R/lambda)
            target_echo=a_r.*exp(1i*pi*gamma*(t_hat-t0).^2)* exp(-1i*4*pi*R/lambda);
            %信号线性叠加
            signal_of_echo(:,k)=signal_of_echo(:,k)+target_echo;
        end
    end
end

% 作图
% 原始回波数据
figure;
subplot(2,2,1);
imagesc(real(signal_of_echo));
title('（a）实部');
xlabel('距离时域');
ylabel('方位时域');

subplot(2,2,2);
imagesc(imag(signal_of_echo));
title('（b）虚部');
xlabel('距离时域');
ylabel('方位时域');

subplot(2,2,3);
imagesc(abs(signal_of_echo));
title('（c）幅度');
xlabel('距离时域');
ylabel('方位时域');

subplot(2,2,4);
imagesc(angle(signal_of_echo));
title('（d）相位');
xlabel('距离时域');
ylabel('方位时域');

figure,
subplot(1,2,1);
imagesc(abs((fft2(signal_of_echo))));
title('二维频谱幅度');
subplot(1,2,2);
imagesc(angle((fft2(signal_of_echo))));
title('二维频谱相位');
%% Chirp Scaling 算法  

%点目标回波的多普勒
f_a=-PRF/2:PRF/num_azimuth:(PRF/2-PRF/num_azimuth);      %方位向频率（多普勒频率）
%位于载机正前方点目标的回波的多普勒，即最大多普勒
f_aM=2*V/lambda;
%斜视角
sin_theta=f_a/f_aM;
cos_theta=sqrt(1-sin_theta.^2);
%回波沿着快时间方向做线性调频的调频率gamma_e，课本P140页(5.38)
%课本对应公式：gamma_e = 1/(1/gamma-R_B*2*lambda*sin_theta*sin_theta/(c^2*cos_theta*cos_theta*cos_theta))
gamma_e=1./(1/gamma-2*R_s*lambda*sin_theta.^2./(cos_theta.^3*c^2)); 
%cs因子 a_f_a=1/sqrt(1-(f_a/f_aM)^2)-1=1/cos_theta-1
a_f_a=1./cos_theta-1;
%课本对应公式：R(f_a;R_B)=R_B*a_f_a+R_B
R_fa_Rs=R_s+R_s*a_f_a;
%P155页用于改变线调频率尺度的chirp scaling二次相位函数
%课本对应公式：H_1=exp(j*pi*gamma_e*a_f_a*(t_hat-2*R/c)^2)
H_1=exp(1i*pi*(ones(num_range,1)*(gamma_e.*a_f_a)).*(t_hat*ones(1,num_azimuth)-2*ones(num_range,1)*R_fa_Rs/c).^2);
H_1=fftshift(H_1,2);
%% 1. 方位向FFT，将信号变换到快时间-方位频率域
for i=1:num_range
    signal_of_echo(i,:)=fft(signal_of_echo(i,:));
end
figure;
imagesc(abs(signal_of_echo));
title('原始数据变换到距离多普勒域，幅度');
figure;
subplot(1,2,1);
imagesc(abs((signal_of_echo)));
title('距离多普勒域 幅度谱');
subplot(1,2,2);
imagesc(angle((signal_of_echo)));
title('距离多普勒域 相位谱');
%% 2. 将快时间-方位多普勒域的信号与chirp scaling二次相位函数H_1相乘
signal_of_echo=signal_of_echo.*H_1;
figure;
subplot(1,2,1);
imagesc(abs(signal_of_echo));
title('距离多普勒域，chirp scaling操作后，幅度谱');
subplot(1,2,2);
imagesc(angle(signal_of_echo));
title('距离多普勒域，chirp scaling操作后，相位谱');
%% 3. 距离向FFT，将chirp scaling操作后的信号变换到距离频率-方位频率域
for i=1:num_azimuth
    signal_of_echo(:,i)=fft(signal_of_echo(:,i));
end
figure;
subplot(1,2,1);
imagesc(abs(signal_of_echo));
title('将chirp scaling操作后的信号变换到二维频域，幅度谱');
subplot(1,2,2);
imagesc(angle(signal_of_echo));
title('将chirp scaling操作后的信号变换到二维频域，相位谱');
%距离向频率坐标
f_r=(-F_s/2:F_s/num_range:(F_s/2-F_s/num_range))';      
%用于距离压缩，距离徙动校正的相位函数H_2
%课本对应公式：H_2=exp(j*pi*f_r^2/(gamma_e(1+a_f_a)))exp(j*4*pi*R_s*a_f_a*f_r/c)
H_2=exp(1i*pi*f_r.^2./(gamma_e.*(1+a_f_a))).*exp(1i*4*pi*R_s*f_r*a_f_a/c);
H_2=fftshift(H_2);
%% 4. 将距离频率-方位频率域的信号与相位函数H_2相乘，对信号进行距离压缩和距离徙动校正
signal_of_echo=signal_of_echo.*H_2;

s_rc=ifft2(signal_of_echo);
figure;
subplot(1,2,1);
imagesc(abs(s_rc));title('距离压缩和距离徙动校正后的信号，幅度谱');xlabel('方位向');ylabel('距离向');
subplot(1,2,2);
imagesc(angle(s_rc));title('距离压缩和距离徙动校正后的信号，相位谱');xlabel('方位向');ylabel('距离向');
%% 5. 对进行距离向逆傅里叶变换，变换到快时间-方位频率域
for i=1:num_azimuth
    signal_of_echo(:,i)=ifft(signal_of_echo(:,i));
end
figure;
subplot(1,2,1);
imagesc(abs(signal_of_echo));title('距离压缩和距离徙动校正后的信号，幅度谱');xlabel('方位向');ylabel('距离向');
subplot(1,2,2);
imagesc(angle(signal_of_echo));title('距离压缩和距离徙动校正后的信号，相位谱');xlabel('方位向');ylabel('距离向');
%% 6. 用于方位压缩处理和补偿由chirp scaling引起的剩余相位函数H_3
%课本对应公式：H_3=exp(j*2*pi*R_B*sqrt(f_aM^2-f_a^2)/V)exp(j*Theta_delta)
R__B=c/2/F_s*(-num_range/2:num_range/2-1)'+R_s;  %距离向采样
%由于chirp scaling操作引起的剩余相位Theta_delta(f_a;R_B)
Theta_delta=4*pi*(R__B-R_s).^2*(gamma_e.*a_f_a.*(1+a_f_a))/c^2;
H_3=exp(1j*2*pi/V*R__B*sqrt(f_aM.^2-f_a.^2)).*exp(1j*Theta_delta);
H_3=fftshift(H_3,2);
%% 7. 将快时间-方位频率域的信号与H_3相乘，进行方位压缩处理，补偿由chirp scaling引起的剩余相位
signal_of_echo=signal_of_echo.*H_3;
figure;
subplot(1,2,1);
imagesc(abs(signal_of_echo));title('方位脉压和相位补偿后快时间-方位频率域的信号，幅度谱');xlabel('方位向');ylabel('距离向');
subplot(1,2,2);
imagesc(angle(signal_of_echo));title('方位脉压和相位补偿后快时间-方位频率域的信号，相位谱');xlabel('方位向');ylabel('距离向');
%% 8. 将信号进行方位向逆傅里叶变换
for i=1:num_range
    signal_of_echo(i,:)=ifft(signal_of_echo(i,:));
end
figure;
subplot(1,2,1);
imagesc(abs(signal_of_echo));title('方位脉压和相位补偿后最终的成像结果');xlabel('方位向');ylabel('距离向');
subplot(1,2,2);
imagesc(angle(signal_of_echo));title('方位脉压和相位补偿后最终的成像结果');xlabel('方位向');ylabel('距离向');
%% 分辨率分析
Imaging_Result = signal_of_echo;
%插值点数
Interp_num = 8;
num_row = 32;
num_col = 100;
%目标的位置，横纵坐标
target_position_row = num_range/2+1;
traget_position_col = num_azimuth/2+1;
data_0 = Two_Dim_Interpolate(Imaging_Result(target_position_row-num_row/2+1:target_position_row+num_row/2,traget_position_col-num_col/2+1:traget_position_col+num_col/2),Interp_num);
%绘制成像结果的等高线
figure,contour(abs(data_0),30)
title('CS算法最终的成像结果');xlabel('方位向');ylabel('距离向');
%% 对距离向进行切片，绘制sinc函数波形，计算sinc函数3db带宽和距离向分辨率，进行分辨率分析
Interp_num = 8;
num_row = 64;
num_col = 200;
Img_Result_range_ori = Imaging_Result(target_position_row-num_row/2+1:target_position_row+num_row/2,traget_position_col);
%补零插值
Img_Res_range_irp = zeros(num_row*Interp_num,1);
Img_Res_range_irp(num_row*Interp_num/2-num_row/2+1:num_row*Interp_num/2+num_row/2,1) = fftshift(fft(fftshift(Img_Result_range_ori)));
Img_Res_range_irp = fftshift(ifft(fftshift(Img_Res_range_irp)));
figure,plot(20*log10(abs(Img_Res_range_irp)./max(abs(Img_Res_range_irp))))
title('CS算法成像的距离向sinc函数波形');xlabel('距离向');ylabel('幅度（dB）');
disp('PSLR：')
PSLR(Img_Res_range_irp)
disp('ISLR：')
ISLR(Img_Res_range_irp)
disp('距离向分辨率理论值： （m）')
c/2/B
disp('距离向分辨率实际值： （m）')
1.2*IRW(Img_Res_range_irp)/Interp_num*(c/2/F_s)     
%% 对方位向进行切片，绘制sinc函数波形，计算sinc函数3db带宽和方位向分辨率，进行分辨率分析
Img_Result_azimuth_ori = Imaging_Result(target_position_row,traget_position_col-num_col/2+1:traget_position_col+num_col/2);
%补零插值
Img_Res_azimuth_irp = zeros(1,num_col*Interp_num);
Img_Res_azimuth_irp(1,num_col*Interp_num/2-num_col/2+1:num_col*Interp_num/2+num_col/2) = fftshift(fft(fftshift(Img_Result_azimuth_ori)));
Img_Res_azimuth_irp = fftshift(ifft(fftshift(Img_Res_azimuth_irp)));
figure,plot(20*log10(abs(Img_Res_azimuth_irp)./max(abs(Img_Res_azimuth_irp))))
title('CS算法成像的方位向sinc函数波形');xlabel('方位向');ylabel('幅度（dB）');
disp('PSLR：')
PSLR(Img_Res_azimuth_irp)
disp('ISLR：')
ISLR(Img_Res_azimuth_irp)
disp('方位向分辨率理论值： （m）')
%lambda/(2*(V*num_azimuth/(PRF*R_s)))              % lambda/(2*theta_BW)
D_a/2
disp('方位向分辨率实际值： （m）')
1.2*IRW(Img_Res_azimuth_irp)/Interp_num*(V/PRF)