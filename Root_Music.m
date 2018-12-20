clear all;
close all;

K = 2;                                                                              %��Դ��
M = 8;                                                                              %��Ԫ��
L = 200;                                                                            %�źų���
w = [pi/4, pi/6].';                                                                 %�ź�Ƶ��
lamda = ((2 * pi * 3e8) / w(1) + (2 * pi * 3e8) / w(2)) / 2;                        %�źŲ���
d_lamda = 0.5;                                                                      %��Ԫ���
snr = 20;                                                                           %�����
theta1 = [49, 60];                                                                  %�ź������

for k = 1: K
    A(:, k) = exp(-1j * [0: M-1]' * d_lamda * 2 * pi * sin(theta1(k) * pi / 180));  %��������
end
for kk = 1: L
    s(:, kk) = sqrt(10.^((snr / 2) / 10)) * exp(1j * w * (kk - 1));                 %�����ź�
end
x = A * s + (1 / sqrt(2)) * (randn(M, L) + 1j * randn(M, L));                       %�����˹������
R = (x * x') / L;                                                                   %Э�������

%%%%%%��һ�ַ���%%%%%%%%%%
[V, D] = eig(R);                                                                    %��Э���������������ֽ�?����û����
Un = V(:, 1: M-K);                                                                  %ȡ�����ӿռ�
Gn = Un * Un';
a = zeros(2 * M - 1, 1)';                                                           %�ҳ�����ʽ��ϵ�������������Ӹ���������
for i = -(M - 1): (M - 1)
    a(i + M) = sum(diag(Gn, i));
end
a1 = roots(a);                                                                      %ʹ��ROOTS�����������ʽ�ĸ�
a2 = a1(abs(a1) < 1);                                                               %�ҳ��ڵ�λԲ������ӽ���λԲ��N����
[lamda, I] = sort(abs(abs(a2) -1 ));                                                %��ѡ����ӽ���λԲ��N����
f = a2(I(1: K));                                                                    %�����źŵ��﷽���
source_doa = [asin(angle(f(1)) / pi) * 180 / pi asin(angle(f(2)) / pi) * 180 / pi];
source_doa = sort(source_doa);
disp('source_doa');
disp(source_doa);

%%%%%%%�ڶ��ַ���%%%%%%%%%
% [V, D] = eig(R);
% Un = V(:, 1: M - K);                                                              %(4.5.7��
% Un1 = Un(1: M, :);
% Un2 = Un(K + 1: M, :);
% T = [1 0 0 0 0 0]';                                                               %����Ԫ��-��Դ����*1
% c = Un1 * inv(Un2) * T;                                                           % (K*(M-K))*((M-K)*(M-K))*((M-K)*1)=K*1
% c = [1, c(2, 1), c(1, 1)];                                                        %��4.5.8��
% f = roots(c);
% source_doa = [asin(angle(f(1)) / pi) * 180 / pi, asin(angle(f(2)) / pi) * 180 / pi];
% source_doa = sort(source_doa);
% disp('source_doa');
% disp(source_doa);
