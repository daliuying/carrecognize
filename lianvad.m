function [v_Begin, v_End] =lianvad(x)  
 %================================i=========================  
 % 端点检测  
 % Input:音频数据x,采样率fs  
 % Output：经过端点检测提取的语音信号  
 %=========================================================  
  
%幅度归一化到[-1,1] 
%[x fs]=audioread('D:\软件安全下载目录\两辆车.wav');
x = x / max(abs(x));%幅度归一化到[-1,1] 
  
%常数设置  
FrameLen = 960;%帧长为256点  
FrameInc = 480;%帧移为80点  
amp1 = 14000;%初始短时能量高门限  
amp2 = 12000;%初始短时能量低门限  
zcr1 = 14000;%初始短时过零率高门限  
zcr2 = 12000;%初始短时过零率低门限  
maxsilence = 10;  % 8*10ms  = 80ms  
  
%语音段中允许的最大静音长度，如果语音段中的静音帧数未超过此值，则认为语音还没结束；如果超过了  
%该值，则对语音段长度count进行判断，若count<minlen，则认为前面的语音段为噪音，舍弃，跳到静音  
%状态0；若count>minlen，则认为语音段结束；  
  
minlen  = 15;    % 15*10ms = 150ms  
%语音段的最短长度，若语音段长度小于此值，则认为其为一段噪音  
  
  
status  = 0;     %初始状态为静音状态  
count   = 0;     %初始语音段长度为0  
silence = 0;     %初始静音段长度为0  
  
%计算过零率  
x1=x(1:end-1);  
x2=x(2:end);  
%分帧  
tmp1=enframe(x1,FrameLen,FrameInc);  
tmp2=enframe(x2,FrameLen,FrameInc);  
signs = (tmp1.*tmp2)<0;  
diffs = (tmp1 -tmp2)>0.02;  
zcr   = sum(signs.*diffs, 2);%一帧一个值，每帧一个过零率，存到zcr数组里
  
%计算短时能量  
%一帧一个值  
%amp = sum(abs(enframe(filter([1 -0.9375], 1, x), FrameLen, FrameInc)), 2);  
amp = sum((abs(enframe( x, FrameLen, FrameInc))).^2, 2);  
  
%调整能量门限  
  
amp1 = min(amp1, max(amp)/4);  
amp2 = min(amp2, max(amp)/8);  
   
  
%开始端点检测  
%For循环，整个信号各帧比较  
%根据各帧能量判断帧所处的阶段  
x1 = 0;  
x2 = 0;  
v_num=0;%记录语音段数  
v_Begin=[];%记录所有语音段的起点  
v_End=[];%记录所有语音段的终点  
  
%length(zcr)即为帧数  
for n=1:length(zcr)  
   goto = 0;  
   switch status  
   case {0,1}                   % 0 = 静音, 1 = 可能开始  
      if amp(n) > amp1          % 确信进入语音段  
         x1 = max(n-count-1,1);  
%          '打印每个x1*FrameInc'  
%          x1*FrameInc  
         status  = 2;  
         silence = 0;  
         count   = count + 1;  
      elseif amp(n) > amp2 | ... % 可能处于语音段  
             zcr(n) > zcr2  
         status = 1;  
         count  = count + 1;  
      else                       % 静音状态  
         status  = 0;  
         count   = 0;  
      end  
   case 2,                       % 2 = 语音段  
      if amp(n) > amp2 | ...     % 保持在语音段  
         zcr(n) > zcr2  
         count = count + 1;  
      else                       % 语音将结束  
         silence = silence+1;  
         if silence < maxsilence % 静音还不够长，尚未结束  
            count  = count + 1;  
         elseif count < minlen   % 语音长度太短，认为是噪声  
            status  = 0;  
            silence = 0;  
            count   = 0;  
         else                    % 语音结束  
            status  = 3;  
         end  
      end  
   case 3,  
      %break;  
      %记录当前语音段数据  
      v_num=v_num+1;   %语音段个数加一  
      count = count-silence/2;  
      x2 = x1 + count -1;  
      v_Begin(1,v_num)=x1*FrameInc;   
      v_End(1,v_num)=x2*FrameInc;  
      %不跳出 数据归零继续往下查找下一段语音  
      status  = 0;     %初始状态为静音状态  
      count   = 0;     %初始语音段长度为0  
      silence = 0;     %初始静音段长度为0  
  
  
   end  
end    
  
if length(v_End)==0  
    x2 = x1 + count -1;  
    v_Begin(1,1)=x1*FrameInc;   
    v_End(1,1)=x2*FrameInc;  
end  
  
 subplot(311)    %subplot(3,1,1)表示将图排成3行1列，最后的一个1表示下面要画第1幅图  
 plot(x)  
 axis([1 length(x) -1 1])    %函数中的四个参数分别表示xmin,xmax,ymin,ymax，即轴的范围  
 ylabel('Speech');  
 for k=1:length(v_End)  
 line([v_Begin(1,k) v_Begin(1,k)], [-1 1], 'Color', 'red');  
  
%这里作用为用直线画出语音段的起点和终点，看起来更直观。第一个[]中的两个参数为线起止点的横坐标，  
  
%第二个[]中的两个参数为线起止点的纵坐标。最后两个参数设置了线的颜色。  
 line([v_End(1,k) v_End(1,k)], [-1 1], 'Color', 'red');  
 end  
%   
% subplot(312)     
% plot(amp);  
% axis([1 length(amp) 0 max(amp)])  
% ylabel('Energy');  
% line([x1 x1], [min(amp),max(amp)], 'Color', 'red');  
% line([x2 x2], [min(amp),max(amp)], 'Color', 'red');  
%   
% subplot(313)  
% plot(zcr);  
% axis([1 length(zcr) 0 max(zcr)])  
% ylabel('ZCR');  
% line([x1 x1], [min(zcr),max(zcr)], 'Color', 'red');  
% line([x2 x2], [min(zcr),max(zcr)], 'Color', 'red');    
