#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fcntl.h>    // 用於檔案控制選項
#define PI 3.14159265358979

typedef struct AudioHeader
{
    char riffID[4];
    int fileSize;
    char waveFormat[4];
    char fmtID[4];
    int fmtChunkSize;
    short audioFormat;
    short numChannels;
    int sampleRate;
    int byteRate;
    short blockAlign;
    short bitsPerSample;
    char dataID[4];
    int dataSize;
} AudioHeader;
double calcSQNR8(int sampleCount, unsigned char *samples8Bit, char *waveType, double *originalSignal, int channels);
double calcSQNR16(int sampleCount, short *samples16Bit, char *waveType, double *originalSignal, int channels);
double calcSQNR32(int sampleCount, int *samples32Bit, char *waveType, double *originalSignal, int channels);
int main(int argc, char *argv[]){

    AudioHeader header; // 音訊檔案的標頭資訊
    char *endPtr;
    int i;
    int j;
    int sampleRate = atoi(argv[1]);
    int bitDepth = atoi(argv[2]);// 位元深度
    int channels = atoi(argv[3]);//聲道數
    char *waveType = argv[4];  // 波形類型
    int frequency = atoi(argv[5]);
    double amplitude = strtod(argv[6], &endPtr); // 振幅
    int duration = atoi(argv[7]); // 持續時間
    double f = frequency; // 頻率的雙精度浮點表示
    size_t sampleCount = (size_t)(duration * sampleRate); // 總取樣數
    int N= sampleCount ; // 總取樣數（整數表示）

    double sqnr = 0.0;// 信號與噪聲比
 // 動態分配記憶體，用於儲存取樣資料和原始信號
    unsigned char *samples8Bit = (unsigned char *)malloc(sizeof(unsigned char) * sampleCount * channels);
    short *samples16Bit = (short *)malloc(sizeof(short) * sampleCount * channels);
    int *samples32Bit = (int *)malloc(sizeof(int) * sampleCount * channels);
    double *originalSignal = (double *)malloc(sizeof(double) * sampleCount * channels);
    if (bitDepth == 8)
    {
        amplitude *= __SCHAR_MAX__;
        if (strcmp(waveType, "sine") == 0)  //m=8 sine
        {
            for (int j = 0; j < N; j++)// N 是樣本數量，循環處理每個樣本
            { // 加上 __SCHAR_MAX__ 是為了將正弦波的範圍平移至 unsigned char 的範圍 (0-255)
                double tmpValue = amplitude * sin(2 * PI * frequency * j / sampleRate) + __SCHAR_MAX__;
                for (int ch = 0; ch < channels; ch++) {   // 對於每個聲道，將計算出的值轉換為無符號 8 位元整數並儲存
                // 使用 floor() 函數將浮點值取整，並將其轉換為 unsigned char 型態
            // 調整了加了 0.5 的偏移量來進行取整操作（這樣做是為了控制四捨五入的效果）
                    samples8Bit[j* channels + ch] = (unsigned char)floor(tmpValue + 0.5);
                    originalSignal[j* channels + ch] = tmpValue;
                }
            }
        }
        else if (strcmp(waveType, "square") == 0)
        { 
            for (int i = 0; i < N; i++)
            {
                double time = (double)i / sampleRate;
                double tmpValue = amplitude * ((time * f - floor(0.5 + time * f) < 0.0) ? -1.0 : 1.0) + __SCHAR_MAX__;
                // 加上 __SCHAR_MAX__ 是為了將方波的範圍平移至 unsigned char 的範圍 (0-255)
                // floor(0.5 + time * f) 是將時間轉換為一個周期內的離散值
                for (int ch = 0; ch < channels; ch++) { 
                    samples8Bit[i* channels + ch] = (unsigned char)floor(tmpValue + 0.25 + 0.25);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        else if (strcmp(waveType, "triangle") == 0)
        {
            for (int i = 0; i < N; i++)
            {
                double time = (double)i / sampleRate;
                // fmod(time * f + 0.75, 1.0) 用來將時間映射到周期內的範圍 (0 到 1)，
                // 然後減去 0.5，將範圍從 [-0.5, 0.5] 轉換到三角波的範圍 [0, 1]。
                // 公式 1.0 - 4.0 * fabs(...) 用來生成三角波的形狀
                // amplitude 是振幅，將三角波的範圍映射到指定的振幅範圍。
                double tmpValue = amplitude * (1.0 - 4.0 * fabs(fmod(time * f + 0.75, 1.0) - 0.5)) + __SCHAR_MAX__;
                for (int ch = 0; ch < channels; ch++) {    
                    samples8Bit[i* channels + ch] = (unsigned char)floor(tmpValue + 0.25 +0.25); // 使用 floor() 函數將浮點值取整，並將其轉換為 unsigned char 型態
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        
        else if (strcmp(waveType, "sawtooth") == 0)
        {
            for (int i = 0; i < N; i++)
            {
                double time = (double)i / sampleRate;
        // floor(0.5 + time * f) 用來將時間映射到周期範圍內（將時間轉為整數周期位置）
        // 2.0 * (...) 則產生鋸齒波的形狀，將其值放大到所需的範圍。
        // __SCHAR_MAX__ 是無符號 char 類型的最大值，將其加到計算出的值上來調整範圍。
                double tmpValue = amplitude * (2.0 * (time * f - floor(0.5 + time * f))) + __SCHAR_MAX__;
                for (int ch = 0; ch < channels; ch++) {
                    samples8Bit[i* channels + ch] = (unsigned char)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
    sqnr = calcSQNR8(sampleCount, samples8Bit, waveType, originalSignal, channels);
    }
    else if (bitDepth == 16)
    {
        amplitude *= __INT16_MAX__;
        if (strcmp(waveType, "sine") == 0)
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double tmpValue = amplitude * sin(2 * PI * frequency * i / sampleRate);
                for (int ch = 0; ch < channels; ch++) {
                    samples16Bit[i* channels + ch] = (short)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        else if (strcmp(waveType, "square") == 0)
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double time = (double)i / sampleRate;
                double tmpValue = amplitude * ((time * f - floor(0.5 + time * f) < 0.0) ? -1.0 : 1.0);
                for (int ch = 0; ch < channels; ch++) {
                    samples16Bit[i* channels + ch] = (short)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        else if (strcmp(waveType, "triangle") == 0)
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double time = (double)i / sampleRate;
                double tmpValue = amplitude * (1.0 - 4.0 * fabs(fmod(time * f + 0.75, 1.0) - 0.5));
                for (int ch = 0; ch < channels; ch++) {
                    samples16Bit[i* channels + ch] = (short)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        
        else if (strcmp(waveType, "sawtooth") == 0)
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double time = (double)i / sampleRate;
                double tmpValue = amplitude * (2.0 * (time * f - floor(0.5 + time * f)));
                for (int ch = 0; ch < channels; ch++) {
                    samples16Bit[i* channels + ch] = (short)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
    sqnr = calcSQNR16(sampleCount, samples16Bit, waveType, originalSignal, channels);
    }
    
    else if (bitDepth == 32)
    {
        amplitude *= __INT32_MAX__;  // 調整振幅為32位元的最大值
        if (strcmp(waveType, "sine") == 0)             // 生成正弦波
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double tmpValue = amplitude * sin(2 * PI * frequency * i / sampleRate);
                for (int ch = 0; ch < channels; ch++) {
                    samples32Bit[i* channels + ch] = (int)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        else if (strcmp(waveType, "square") == 0)             // 產生方波
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double time = (double)i / sampleRate;
                double tmpValue = amplitude * ((time * f - floor(0.5 + time * f) < 0.0) ? -1.0 : 1.0);
                for (int ch = 0; ch < channels; ch++) {
                    samples32Bit[i* channels + ch] = (int)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        else if (strcmp(waveType, "triangle") == 0)             // 生成三角波
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double time = (double)i / sampleRate;
                double tmpValue = amplitude * (1.0 - 4.0 * fabs(fmod(time * f + 0.75, 1.0) - 0.5));
                for (int ch = 0; ch < channels; ch++) {
                    samples32Bit[i* channels + ch] = (int)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        else if (strcmp(waveType, "sawtooth") == 0)             // 產生鋸齒波
        {
            for (int i = 0; i < sampleCount; i++)
            {
                double time = (double)i / sampleRate;
                double tmpValue = amplitude * (2.0 * (time * f - floor(0.5 + time * f)));
                for (int ch = 0; ch < channels; ch++) {
                    samples32Bit[i* channels + ch] = (int)floor(tmpValue + 0.5);
                    originalSignal[i* channels + ch] = tmpValue;
                }
            }
        }
        
    sqnr = calcSQNR32(sampleCount, samples32Bit, waveType, originalSignal, channels);
    }
    
    // 寫入音訊檔頭資訊
    header.riffID[0] = 'R'; header.riffID[1] = 'I'; header.riffID[2] = 'F'; header.riffID[3] = 'F'; // 設定 "RIFF" 標誌，表明檔案是 RIFF 格式的音訊檔案
    header.fileSize = sampleCount * channels * (bitDepth / 8) + 36;// 計算檔案的總大小 (檔案大小 = 樣本數 * 聲道數 * (位元深度/8) + 36)
    header.waveFormat[0] = 'W'; header.waveFormat[1] = 'A'; header.waveFormat[2] = 'V'; header.waveFormat[3] = 'E'; // 設定 "WAVE" 標誌，表明檔案是 WAV 格式
    header.fmtID[0] = 'f'; header.fmtID[1] = 'm'; header.fmtID[2] = 't'; header.fmtID[3] = ' '; // 設定 "fmt " 標誌，表示格式區塊。
    header.fmtChunkSize = 16; // 設定格式區塊的大小，固定為 16。
    header.audioFormat = 1; // 設定音訊格式為 PCM（1 表示未壓縮的 PCM 格式）。
    header.numChannels = channels; // 設定聲道數。
    header.sampleRate = sampleRate;// 設定取樣率。
    header.bitsPerSample = bitDepth; // 設定每個樣本的位元深度（如 8、16 或 32 位元）。
    header.byteRate = header.sampleRate * header.numChannels * (header.bitsPerSample / 8); // 計算每秒的位元數率 (取樣率 * 聲道數 * (位元深度 / 8))。
    header.blockAlign = (header.numChannels * header.bitsPerSample) / 8; // 計算每個音訊幀的字節數 (每個聲道的字節數)。
    header.dataID[0] = 'd'; header.dataID[1] = 'a'; header.dataID[2] = 't'; header.dataID[3] = 'a'; // 設定 "data" 標誌，表示音訊數據區塊的開始。
    header.dataSize = sampleCount * channels *(bitDepth / 8); // 計算音訊數據區塊的大小 (樣本數 * 聲道數 * (位元深度 / 8))。

    // 寫入檔案頭資訊到標準輸出（stdout）
    fwrite(&header.riffID, sizeof(char), 4, stdout);
    fwrite(&header.fileSize, sizeof(int), 1, stdout);
    fwrite(&header.waveFormat, sizeof(char), 4, stdout);
    fwrite(&header.fmtID, sizeof(char), 4, stdout);
    fwrite(&header.fmtChunkSize, sizeof(int), 1, stdout);
    fwrite(&header.audioFormat, sizeof(short), 1, stdout);
    fwrite(&header.numChannels, sizeof(short), 1, stdout);
    fwrite(&header.sampleRate, sizeof(int), 1, stdout);
    fwrite(&header.byteRate, sizeof(int), 1, stdout);
    fwrite(&header.blockAlign, sizeof(short), 1, stdout);
    fwrite(&header.bitsPerSample, sizeof(short), 1, stdout);
    fwrite(&header.dataID, sizeof(char), 4, stdout);
    fwrite(&header.dataSize, sizeof(int), 1, stdout);
    
    // 根據不同的位元深度選擇對應的數據類型進行寫入。
    if (bitDepth == 8)
    {
        fwrite(samples8Bit, sizeof(unsigned char), sampleCount * channels, stdout);
    }
    else if (bitDepth == 16)
    {
        fwrite(samples16Bit, sizeof(short), sampleCount * channels, stdout);
    }
    else if (bitDepth == 32)
    {
        fwrite(samples32Bit, sizeof(int), sampleCount * channels, stdout);
    }

    
    free(samples8Bit);
    free(samples16Bit);
    free(samples32Bit);
    free(originalSignal);

    return 0;
}

// 計算8位樣本的信噪比
double calcSQNR8(int sampleCount, unsigned char *samples8Bit, char *waveType, double *originalSignal, int channels)
{
    double error = 0.0, signalPower = 0.0, noisePower = 0.0, sqnr = 0.0;  // 初始化誤差、信號功率、噪聲功率和SQNR
    for (int i = 0; i < sampleCount* channels; i++)  // 計算誤差並累加信號功率和噪聲功率
    {
        error = originalSignal[i] - samples8Bit[i];  // 計算每個樣本的誤差（即原始信號與經過編碼後的樣本之間的差異）
        signalPower = signalPower + (originalSignal[i] - 128) * (originalSignal[i] - 128); // 計算信號的功率。這裡假設原始信號的範圍是 -128 到 127，並將其偏移 128 (中心化)
        noisePower += error * error; // 計算噪聲的功率，即誤差的平方和
    }// 計算平均信號功率和噪聲功率
    signalPower /= (sampleCount * channels);
    noisePower /= (sampleCount * channels);
    sqnr = 10 * log10(signalPower / noisePower);  // 計算 SQNR，單位是分貝 (dB)
    fprintf(stderr, "%.15f\n", sqnr); // 輸出 SQNR 的值到stderr
    return sqnr;
}
// 計算16位樣本的信噪比
double calcSQNR16(int sampleCount, short *samples16Bit, char *waveType, double *originalSignal, int channels)
{
    double error = 0.0, signalPower = 0.0, noisePower = 0.0, sqnr = 0.0;

    for (int i = 0; i < sampleCount* channels; i++)
    {
        error = originalSignal[i] - samples16Bit[i]; // 計算每個樣本的誤差（即原始信號與經過編碼後的樣本之間的差異）          
        signalPower = signalPower + originalSignal[i] * originalSignal[i];  // 計算信號的功率            
        noisePower = noisePower + error * error;  // 計算噪聲的功率，即誤差的平方和            
    }

    signalPower = signalPower / (sampleCount* channels);   
    noisePower = noisePower / (sampleCount* channels); 
    sqnr = 10 * log10(signalPower / noisePower);//轉換為分貝 dB)可以量化信號的質量
   
    fprintf(stderr, "%.15f\n", sqnr); 
    return sqnr;
} 
// 計算32位樣本的信噪比
double calcSQNR32(int sampleCount, int *samples32Bit, char *waveType, double *originalSignal, int channels)
{
    double error = 0.0, signalPower = 0.0, noisePower = 0.0, sqnr = 0.0;
    
    for (int i = 0; i < sampleCount* channels; i++)    
    {
        error = originalSignal[i] - samples32Bit[i]; 
        noisePower = noisePower + error * error;         
        signalPower = signalPower + originalSignal[i] * originalSignal[i];         
    }

    signalPower = signalPower / (sampleCount* channels);   
    noisePower /= (sampleCount * channels); 
    sqnr = 10 * log10(signalPower / noisePower);
    
    fprintf(stderr, "%.15f\n", sqnr); 
    return sqnr;
}
