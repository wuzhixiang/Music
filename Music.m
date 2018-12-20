clear all
close all

derad = pi / 180;
radeg = 180 / pi;
twpi = 2 * pi;

kelm = 8;                                          %��Ԫ����
dd = 0.5;                                          %��Ԫ���
d = 0: dd: (kelm - 1) * dd;

iwave = 3;                                         %�������ź�Դ��
theta = [10 30 60];                                %�����ƽǶ�Ϊ10 30 60
snr = 10;                                          %�����
n = 500;                                           %������

A = exp(-j * twpi * d.' * sin(theta * derad));
S = randn(iwave, n);
X = A * S;
X1 = awgn(X, snr, 'measured');
Rxx = X1 * X1' / n;
InvS = inv(Rxx);

[EV, D] = eig(Rxx);                                %����DΪ����ֵ�ĶԽǾ���
EVA = diag(D)';                                    %����Ϊ����ֵ����ɵ�����
[EVA, I] = sort(EVA);
EVA = fliplr(EVA);                                 %����ֵ�Ӵ�С����
EV = fliplr(EV(:, I));                             %����������������ֵ����

for iang = 1:361                                   %��ÿ���Ƕ�����׹���
	angle(iang) = (iang - 181) / 2;
	phim = derad * angle(iang);
	a = exp(-j * twpi * d * sin(phim)).';
	L = iwave;
	En = EV(:, L+1:kelm);                          %�������ӿռ��ͶӰ�������
	SP(iang) = (a' * a) / (a' * En * En' * a);
end
SP = abs(SP);
SPmax = max(SP);                                   %ȡ�׹��Ƶ����ֵ
SP = 10 * log10(SP / SPmax);                       %dB = 10lg(p/po)     (p:���ʣ�po����׼����)

h = plot(angle, SP);
set(h, 'Linewidth', 2)
xlabel('angle (degree)')
ylabel('magnitude (dB)')
axis([-90 90 -60 0])
set(gca, 'XTick', [-90:30:90])
grid on
