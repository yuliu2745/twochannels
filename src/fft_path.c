#include "../include/setting.h"
#include "../include/fft_path.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>




//查找互相关峰值
float find_nbest_maximums(int *delays, float *values, int amount_max, float *xcorr_value, int margin, int fftwindow_size)
{
    (void)fftwindow_size; // 未使用，防止警告

    // 屏蔽范围
    int masking = 5;

    // 最大峰值数量：最多 2*margin 个
    int max_count = 0;
    float max_values[2048];  // 足够大的数组
    int max_pos[2048];
    float xcorr_sum = xcorr_value[0];

    // 查找所有峰值（跳过 i=0）
    for (int i = 1; i < (2 * margin) - 1; i++)
    {
        xcorr_sum += xcorr_value[i];
        if (xcorr_value[i] > xcorr_value[i-1] && xcorr_value[i] > xcorr_value[i+1])
        {
            max_values[max_count] = xcorr_value[i];
            max_pos[max_count] = i;
            max_count++;
        }
    }

    // 如果没有找到任何峰值，手动添加一个
    if (max_count == 0)
    {
        max_values[0] = 0.0001f;
        max_pos[0] = margin;
        max_count = 1;
        printf("WARNING!!! no maximum points found in crosscorrelation!!!\n");
    }

    // 提取前 N 个最大峰值
    float max_value;
    int max_idx, max_idx_here, good_max;

    for (int i = 0; i < amount_max; i++)
    {
        max_value = -1.0f;
        max_idx = -1;
        max_idx_here = -1;

        // 找当前最大值
        for (int count = 0; count < max_count; count++)
        {
            if (max_values[count] > max_value)
            {
                max_value = max_values[count];
                max_idx = max_pos[count];
                max_idx_here = count;
            }
        }

        if (max_value != -1.0f)
        {
            good_max = 1;

            // 检查是否与已选峰值太近
            for (int j = 0; j < i; j++)
            {
                int delay_pos = delays[j] + margin + 1;
                if (max_idx - masking < delay_pos && max_idx + masking > delay_pos)
                {
                    good_max = 0;
                    max_values[max_idx_here] = -1.0f;
                }
            }

            if (!good_max)
            {
                i--; // 重新找一个
            }
            else
            {
                // 计算最终延时
                delays[i] = max_idx - margin;
                values[i] = max_value;
                max_values[max_idx_here] = -1.0f;
            }
        }
        else
        {
            // 没有更多峰值，复制第一个
            delays[i] = delays[0];
            values[i] = values[0];
        }
    }

    if (max_idx == -1 || max_idx_here == -1)
    {
        printf("ERROR: Maximum index not found\n");
        exit(1);
    }

    return xcorr_sum;
}


// FFT-PHAT实现GCC_PATH  找到的延时值     互相关峰值           峰值数量     通道2数据    通道1数据（作为参考数据） 搜索延迟的范围    窗口大小   FFT大小
float FFT_Real_Gcc_Path(int *delays, float *peak_values, int *peak_num, float *channel2, float *channel1, int margin, int window, int fft_size)
{
    //变量声明、内存分配、创建FFT计划
    int i = 0;                                                                      //循环计数，用于加窗、FFT、数据移动等
    int fft_real_size = 2 * fft_size;                                               //实数FFT数组大小
    float *xcorr_value;                                                             //指向存放互相关函数的最终结果的地址,在+margin到-margin之间的值
    xcorr_value = (float *)malloc((2 * margin + 1) * sizeof(float));                //申请连续的内存空间存放互相关函数的最终结果
    int complex_size = fft_size;                                                //复数数组大小,使用实际 FFT 大小为 2*fft_size

    //分配输入实数缓冲区（长度 fft_real_size）和输出复数缓冲区（长度 complex_size）
    float *in1 = (float*)fftwf_malloc(sizeof(float) * fft_real_size);
    float *in2 = (float*)fftwf_malloc(sizeof(float) * fft_real_size);
    fftwf_complex *out1 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * complex_size);               //参考信号channel1的FFT结果
    fftwf_complex *out2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * complex_size);               //通道信号channel2的FFT结果

    // 分配互功率谱复数缓冲区（用于逆变换）
    fftwf_complex *cross_spectrum = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * complex_size);
    // 分配逆变换结果实数缓冲区
    float *xcorr = (float*)fftwf_malloc(sizeof(float) * fft_real_size);

    if (!in1 || !in2 || !out1 || !out2 || !cross_spectrum || !xcorr) 
    {
        // 内存分配失败，清理并返回 0
        if (in1) fftwf_free(in1);
        if (in2) fftwf_free(in2);
        if (out1) fftwf_free(out1);
        if (out2) fftwf_free(out2);
        if (cross_spectrum) fftwf_free(cross_spectrum);
        if (xcorr) fftwf_free(xcorr);
        *peak_num = 0;
        return 0.0f;
    }

    // 创建 FFTW 计划（实数→复数，长度 fft_real_size)
    fftwf_plan plan_r2c1 = fftwf_plan_dft_r2c_1d(fft_real_size, in1, out1, FFTW_ESTIMATE);
    fftwf_plan plan_r2c2 = fftwf_plan_dft_r2c_1d(fft_real_size, in2, out2, FFTW_ESTIMATE);
    // 创建逆变换计划（复数→实数）
    fftwf_plan plan_c2r = fftwf_plan_dft_c2r_1d(fft_real_size, cross_spectrum, xcorr, FFTW_ESTIMATE);

    if (!plan_r2c1 || !plan_r2c2 || !plan_c2r) 
    {
        // 计划创建失败，清理资源
        fftwf_destroy_plan(plan_r2c1);
        fftwf_destroy_plan(plan_r2c2);
        fftwf_destroy_plan(plan_c2r);
        fftwf_free(in1); fftwf_free(in2); fftwf_free(out1); fftwf_free(out2);
        fftwf_free(cross_spectrum); fftwf_free(xcorr);
        *peak_num = 0;
        return 0.0f;
    }

    float *channel_in, *refer_in, *result, *channelFFT, *refFFT;
    result = (float*)malloc(sizeof(float) * 2*fft_size);                                  //存放逆 FFT 的结果，即时域互相关序列。后续会从中提取出 [-margin, margin] 部分
    channel_in = (float*)malloc(sizeof(float) * 2*fft_size);                              //存放经过加窗和补零后的通道信号（channel2），作为 FFT 的输入
    refer_in = (float*)malloc(sizeof(float) * 2*fft_size);                                //存放经过加窗和补零后的参考信号（channel1），作为 FFT 的输入
    channelFFT = (float*)malloc(sizeof(float) * 2*fft_size);                              //存放 channel_in 的 FFT 结果
    refFFT = (float*)malloc(sizeof(float) * 2*fft_size);                                  //存放 refer_in 的 FFT 结果
    
    float hamm_val;                                                                             //用于存储汉明窗的系数值，在加窗循环中使用
    
    /*对输入信号加汉明窗*/
    for(i=0; i<(window); i++)
    {
        hamm_val = 0.54 - 0.46*cos(6.283185307*i/(window-1));
        channel_in[i] = channel2[i] * hamm_val;
        refer_in[i]   = channel1[i] * hamm_val;
    }
    //0填充到FFT长度
    for(i=window; i<(2*fft_size); i++)
    {
        channel_in[i] = 0;
        refer_in[i] = 0;
    }

    /*执行FFT，时域信号转频域*/
    memcpy(in1, refer_in,  sizeof(float)*fft_real_size);
    memcpy(in2, channel_in,sizeof(float)*fft_real_size);

    fftwf_execute(plan_r2c1); // 参考信号 FFT
    fftwf_execute(plan_r2c2); // 通道信号 FFT

    channelFFT[0] = out2[0][0];  // DC分量
    refFFT[0] = out1[0][0];
    
    for(int i = 1; i < fft_size; i++) {
        channelFFT[i] = out2[i][0];        // 实部
        channelFFT[fft_size + i] = out2[i][1]; // 虚部
        refFFT[i] = out1[i][0];           // 实部  
        refFFT[fft_size + i] = out1[i][1]; // 虚部
    }
    
    channelFFT[fft_size] = out2[fft_size][0]; // Nyquist频率
    refFFT[fft_size] = out1[fft_size][0];

    /*GCC_PATH核心计算*/ 
    float* tmpData = (float*)malloc(sizeof(float) * 2*fft_size);
    float abs_value;
    for(i=0; i < fft_size; i++) {
        //复数乘法
        if(i == 0 || i == fft_size) 
        {
            tmpData[i] = channelFFT[i] * refFFT[i];
            abs_value = fabsf(tmpData[i]);
        } 
        else //复数频率分量
        {
            tmpData[i] = channelFFT[i] * refFFT[i] + channelFFT[complex_size+i] * refFFT[complex_size+i];                       //实部：ac + bd
            tmpData[complex_size + i] = channelFFT[complex_size+i] * refFFT[i] - channelFFT[i] * refFFT[complex_size+i];        //虚部：ad - bc
            abs_value = sqrtf(tmpData[i]*tmpData[i] + tmpData[complex_size+i]*tmpData[complex_size+i]);
        }

        //PATH变换：除以幅度，只保留相位信息
        if(abs_value == 0.0f) abs_value = 1.0f;
        tmpData[i] /= (abs_value * fft_real_size);//归一化
        if(i != 0 && i != fft_size) 
        {
            tmpData[fft_size + i] /= (abs_value * fft_real_size);
        }
    }

    /*逆FFT*/
    //将计算结果复制到cross_spectrum用于逆FFT
    for(int i = 0; i <= fft_size; i++) {
        cross_spectrum[i][0] = tmpData[i];  // 实部
        if(i != 0 && i != fft_size) {
            cross_spectrum[i][1] = tmpData[fft_size + i]; // 虚部
        } else {
            cross_spectrum[i][1] = 0; // DC和Nyquist频率虚部为0
        }
    }
    fftwf_execute(plan_c2r);      // 执行逆FFT
    free(tmpData);                // 释放 tmpData

    /*互相关结果移位：正负延迟对齐*/
    for(i = 0; i < margin; i++)// 正延迟部分 positive part
    {
        xcorr_value[i + (margin + 1)] = xcorr[i];
    }

    for(i = 2 * fft_size - (margin + 1); i < 2 * fft_size; i++)// 负延迟部分 negative part
    {
        int target_index = i - (2 * fft_size - (margin + 1));
        if (target_index >= 0 && target_index < (2 * margin + 1)) {
            xcorr_value[target_index] = xcorr[i];
        }
    }

    /*查找峰值*/
    float xcorr_sum;
    xcorr_sum = find_nbest_maximums(delays, peak_values, *peak_num, xcorr_value, margin, fft_size);


    /*释放所有内存*/

    // 销毁 FFTW 计划
    fftwf_destroy_plan(plan_r2c1);
    fftwf_destroy_plan(plan_r2c2);
    fftwf_destroy_plan(plan_c2r);

    // 释放 FFTW 申请的内存
    fftwf_free(in1);
    fftwf_free(in2);
    fftwf_free(out1);
    fftwf_free(out2);
    fftwf_free(cross_spectrum);
    fftwf_free(xcorr);

    // 释放普通 malloc 内存
    free(channel_in);
    free(refer_in);
    free(xcorr_value);
    free(result);        
    free(channelFFT);    
    free(refFFT);        

    return xcorr_sum;


}
