clc;                % 清除命令窗口的内容
clear all;                % 清除了所有的变量
close all;                % 关闭所有窗口

fs=48000;                % ??采样率，1秒保存的信号数据点数

Tw = 20;                % analysis frame duration (ms)  毫秒，到MFCC中再转换为帧长+
Ts = 10;                % analysis frame shift (ms)  毫秒，到MFCC中再转换为帧移+
alpha = 0.97;           % preemphasis coefficient  预加重系数+
M = 20;                 % number of filterbank channels   Mel滤波器（filterbank channels）个数一般取20~26个
C = 12;                 % number of cepstral coefficients 倒频谱系数
L = 22;                 % cepstral sine lifter parameter  倒谱正弦提升参数
LF = 2400;               % lower frequency limit (Hz)  
HF = 24000;              % upper frequency limit (Hz)

%训练用数据
[train_signals train_labels] = load_audio_from_folder('train_car');%读取语音文件、生成标签数据（1 * 样本数）
%disp(size(train_labels));
%转MFCC参数
MFCC_feature={};               % 存储训练特征的广义矩阵，整体是列向量，一行是一个样本
min1 = 1000000;               % 存储训练特征矩阵最小的维度
min2 = 1000000;               % 存储测试特征矩阵最小的维度    ???????

% 所有样本和标签进行提取特征训练模型,输入信号(1 * 样本数)，得到的MFCC特征矩阵(样本数 * 1)
for speech=train_signals
    [x1 x2] = vad(speech{1});      % 端点检测
    signal  = speech{1}(round(x1):round(x2),:);
    %disp(x2 - x1); % ------------------------------------------------------------------------------
    [ MFCCs, FBEs, frames ] = mfcc( signal, fs, Tw, Ts, alpha, @hamming, [LF HF], M, C+1, L );    
    MFCC_feature(end+1,1)={MFCCs};
    %disp(size(MFCCs)); % ---------------------------------------------------------------------------
    if min1 > size(MFCCs,2)
        min1 = size(MFCCs,2);
    end;
end;
for i = 1: size(MFCC_feature)      % 自己写的，把输入特征全搞成一样的大小
    data1 = cell2mat(MFCC_feature(i,1));
    MFCC_feature(i,1) = {data1(:, 1:min1)};
end;
%测试用数据
[test_signals fs] = audioread('D:\软件安全下载目录\4c3222f5a0.wav');        %读取语音文件，生成测试数据和标签
%转MFCC参数
test_feature={};
[x1 x2] = lianvad(test_signals);    % 端点检测
for k=1:length(x1)
    %display(sprintf('x1chanfdu %d %d',length(x1),k));
    signal  = test_signals(round(x1(1,k)):round(x2(1,k)),:);
    [ MFCCs, FBEs, frames ] = mfcc( signal, fs, Tw, Ts, alpha, @hamming, [LF HF], M, C+1, L );    
    test_feature(end+1,1)={MFCCs};
    if min2 > size(MFCCs,2)
        min2 = size(MFCCs,2);
    end;
end;

for i = 1: size(test_feature)       % 自己写的，把输入特征全搞成一样的大小
    data2 = cell2mat(test_feature(i,1));
    test_feature(i,1) = {data2(:, 1:min2)};
end;
%% HMM模型建立和训练
d = C-3;    % 一个观察符号的维度，MFCC特征系数
nstates = 3;    % 状态个数 -------------------------
nmix    = 3;    % 混合高斯分布个数------------------------

%训练的HMM
% pi0=[1  0  0  0  0]; %初始状态概率
% trans0=[1/2  1/2   0   0     0  ; %转移概率 
%         0    1/2  1/2   0    0  ; 
%         0    0    1/2  1/2   0  ; 
%         0    0     0   1/2  1/2 ; 
%         0    0     0    0    1  ];
 pi0=[1  0  0]; %初始状态概率
trans0=[1/2  1/2   0    ; %转移概率 
        0    1/2  1/2   ; 
        0    0    1   ; 
       ];
    
models={};
[unique_train_labels, ~, indices] = unique(train_labels); % 去掉重复，只取唯一值，~为忽略输出参数,unique_train_labels（1 * 分类数）
% indices（1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2）
% train_labels（large large ... small small)
for i = 1:length(unique_train_labels)
    display(sprintf('modeling %s...', char(unique_train_labels(i))))     % 打印hmm模型名称

    %发射概率初始值计算
    stackedData = cell2mat(MFCC_feature(indices == i)')'; % 训练用数据----------------
    mu = zeros(d, nmix); % 均值
    Sigma = zeros(d, d, nmix); % 方差
    for k = 1:nmix
        XX             = stackedData + randn(size(stackedData));
        mu(:, k)       = colvec(mean(XX)); % 所有训练数据的均值
        Sigma(:, :, k) = cov(XX); % 所有训练数据的方差
    end
    M = normalize(rand(nstates, nmix), 1);%混合高斯系数
    emission = condMixGaussTiedCpdCreate(mu, Sigma, M);%发射状态描述结构体获得
    
    %训练一个个HMM
    %verbose：是否打印信息
    %nRandomRestarts:随机训练几次
    %maxiter：最大迭代次数
    %nmix：混合高斯个数
    %pi0：HMM初始状态概率
    %trans0：HMM初始转移概率
    %emission0：初始发射状态描述
    %piPrior：防止除数为0的情况出现   a/b=a+piPrior/b+piPrior  变成这种形式
    %transPrior：防止除数为0的情况出现   a/b=a+piPrior/b+piPrior  变成这种形式
    models(end+1,1) = {hmmFit(MFCC_feature(indices == i), nstates, 'mixGaussTied', 'verbose', false, ...
                        'nRandomRestarts', 3, 'maxiter', 50, 'nmix', nmix,...
                        'pi0',pi0,'trans0',trans0,'emission0',emission, ...
                        'piPrior',pi0,'transPrior',trans0.*10)}; 
end;

%训练样本识别率计算
errorcount=0;
for j=1:length(MFCC_feature)
    p=zeros(length(unique_train_labels),1);
    for i = 1:length(unique_train_labels)
        p(i)=hmmLogprob(models{i}, MFCC_feature(j)); % 计算每个HMM模型的概率值
    end;
    [~, i]=max(p);%取最大概率的模型作为识别，返回第i个模型
    display(sprintf('"%s" is recognized as "%s"', train_labels{j},char(unique_train_labels(i))))
    if ~strcmp(train_labels{j},char(unique_train_labels(i)))
        errorcount=errorcount+1;%错误累计vv
    end;
end;
display(sprintf('train accuracy is %0.2f', (length(MFCC_feature)-errorcount)*100/length(MFCC_feature)));

%% 测试样本识别率计算
errorcount=0;
for j=1:length(test_feature)
    p=zeros(length(unique_train_labels),1);
    for i = 1:length(unique_train_labels)
        p(i)=hmmLogprob(models{i}, test_feature(j)); % 计算概率值
    end;
    [~, i]=max(p); % 取最大概率的模型作为识别
    display(sprintf('x%d is recognized as "%s"', j,char(unique_train_labels(i))))
    %if ~strcmp(test_labels{j},char(unique_train_labels(i)))
     %   errorcount=errorcount+1; % 错误累计
    %end;
end;
%display(sprintf('test accuracy is %0.2f', (length(test_labels)-errorcount)*100/length(test_labels)));

