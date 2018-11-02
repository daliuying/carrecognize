function [audio_signals, word_labels] = load_audio_from_folder(audio_folder)
    audio_signals = {}; % 一行是一个分类的所有样本
    word_labels = {};
    
    % audio_folder为train_car,dir列出所有文件夹，包括. 和..；word_folder为四个文件夹
    for word_folder = struct2cell(dir(audio_folder))
        
        for word_file = struct2cell(dir(sprintf('%s/%s/*.wav', audio_folder, char(word_folder(1))))) % 得到该文件夹下所有的wav文件名
            file_path = sprintf('%s/%s/%s', audio_folder, char(word_folder(1)), char(word_file(1))); % 得到一个语音文件的路径
            [x fs]=audioread(file_path); % 把一个样本语音信号读入x
            %disp(fs);
            audio_signals(end + 1) = {x(:,1)}; %#ok<AGROW> 1 * 30
            word_labels(end + 1) = word_folder(1); %#ok<AGROW> % 大车或者小车的文件夹名，即大车或小车的标签，1 * 30
        end
    end
end
