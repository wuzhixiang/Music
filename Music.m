clear all
close all

derad = pi / 180;
radeg = 180 / pi;
twpi = 2 * pi;

kelm = 8;                                          %阵元个数
dd = 0.5;                                          %阵元间距
d = 0: dd: (kelm - 1) * dd;

iwave = 3;                                         %待估计信号源数
theta = [10 30 60];                                %待估计角度为10 30 60
snr = 10;                                          %信噪比
n = 500;                                           %快拍数

A = exp(-j * twpi * d.' * sin(theta * derad));
S = randn(iwave, n);
X = A * S;
X1 = awgn(X, snr, 'measured');
Rxx = X1 * X1' / n;
InvS = inv(Rxx);

[EV, D] = eig(Rxx);                                %这里D为特征值的对角矩阵
EVA = diag(D)';                                    %这里为特征值所组成的向量
[EVA, I] = sort(EVA);
EVA = fliplr(EVA);                                 %特征值从大到小排列
EV = fliplr(EV(:, I));                             %特征向量依据特征值排列

for iang = 1:361                                   %对每个角度求出谱估计
	angle(iang) = (iang - 181) / 2;
	phim = derad * angle(iang);
	a = exp(-j * twpi * d * sin(phim)).';
	L = iwave;
	En = EV(:, L+1:kelm);                          %将噪声子空间的投影算子求出
	SP(iang) = (a' * a) / (a' * En * En' * a);
end
SP = abs(SP);
SPmax = max(SP);                                   %取谱估计的最大值
SP = 10 * log10(SP / SPmax);                       %dB = 10lg(p/po)     (p:功率，po：基准功率)

h = plot(angle, SP);
set(h, 'Linewidth', 2)
xlabel('angle (degree)')
ylabel('magnitude (dB)')
axis([-90 90 -60 0])
set(gca, 'XTick', [-90:30:90])
grid on
