function CHAT_demo

% CHAT repository
CHAT_path=which('CHAT_demo.m');

% Sample data folder
fs=filesep;
DirMain=sprintf('%s%s%s%s',fileparts(CHAT_path),fs,'Sample Data Folder',fs);

% Pass DirMain to 'CHAT_ProcessImageData.m'
CHAT_ProcessImageData(DirMain)
