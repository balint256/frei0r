/*
 * Copyright (C) 2010-2011 Simon Andreas Eugster (simon.eu@gmail.com)
 * This file is not a Frei0r plugin but a collection of ideas.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "frei0r.hpp"

// Other includes used for the examples below.
// Can be removed on copy/paste.

// Limits (min/max values) of various data types
#include <limits>
// For the CHAR_BIT constant
#include <climits>
// pow() and other mathematical functions
#include <cmath>

//#include <string>
using namespace std;
#include <stdio.h>
#include <string.h>

/**
This is a sample filter for easy copy/pasting.

The CMakeLists.txt needs to be adjusted as well (both in this directory and in the parent directory).
Also, don't forget src/Makefile.am.
*/

class time0r : public frei0r::filter
{
private:
    //f0r_param_double m_barSize;
    //f0r_param_bool m_pointerMethod;
    std::string m_paramOp;
public:
    enum Operation
    {
        OP_NULL             = 0,    // -
        OP_MIN              = 1,    // -
        OP_MAX              = 2,    // -
        OP_AVE              = 3,    // <window length>
        OP_MINMAX           = 4,    // -
        OP_MINMAXAVE        = 5,    // <ave history (background) factor>:<min/max threshold re-adj factor>:<min/max colour buffer re-adj factor>:<min/max split weight>:<min/max threshold weight>
        OP_AVETHRESH        = 6,    // <window length>:<ave diff threshold>:<skip>
        OP_AVETHRESHDEBUG   = 7,    //
        OP_MINMAXAVEDEBUG   = 8,    //
        OP_HSV              = 9,
    };
private:
    //std::vector<uint8_t> lookupTable;
    //std::vector<uint8_t> additionTable;
    uint8_t *m_pBuffer,
        *m_pBufferMonochrome, *m_pBufferMonochromeMin, *m_pBufferMonochromeMax,
        *m_pBufferMonochromeAve,
        *m_pBufferMin, *m_pBufferMax,
        *m_pBufferHistory;
    uint32_t *m_pBufferAve;
    Operation m_op;
    double m_fOp;
    double m_dAveragePeriod;
    int m_iAveragePeriod;
    int64_t m_iFrameCount;
    int m_iAverageBufferIndex;
    double m_dAverageDiffThreshold;
    float *m_fMonoMin, *m_fMonoMax, *m_fHistory, *m_fMin, *m_fMax, *m_fAve;
    double m_dAverageSmoothing, m_dMinMaxSplitWeight, m_dMinMaxThresholdWeight;
public:
    time0r(unsigned int width, unsigned int height)
        : m_pBuffer(NULL)
        //, m_strOp("null")
        , m_paramOp("null")
        , m_fOp(OP_NULL)
        , m_pBufferMonochrome(NULL)
        , m_pBufferMonochromeAve(NULL)
        , m_pBufferMonochromeMin(NULL)
        , m_pBufferMonochromeMax(NULL)
        , m_pBufferMin(NULL)
        , m_pBufferMax(NULL)
        , m_dAveragePeriod(0)
        , m_iAveragePeriod(0)
        , m_iFrameCount(-1)
        , m_iAverageBufferIndex(0)
        , m_dAverageDiffThreshold(0)
        , m_pBufferAve(NULL)
        , m_pBufferHistory(NULL)
        , m_fMonoMin(NULL)
        , m_fMonoMax(NULL)
        , m_fHistory(NULL)
        , m_dAverageSmoothing(0)
        , m_fMin(NULL)
        , m_fMax(NULL)
        , m_dMinMaxSplitWeight(0.5)
        , m_dMinMaxThresholdWeight(0.5)
        , m_fAve(NULL)
    {
        //register_param(m_barSize, "barSize", "Size of the black bar");
        //register_param(m_pointerMethod, "pointerMethod", "Pointer Method (internal)");
        //register_param(m_paramOp, "op", "Operation");
        register_param(m_fOp, "op", "Operation");
        register_param(m_dAveragePeriod, "ave_len", "Averaging period");
        register_param(m_dAverageDiffThreshold, "ave_diff_thresh", "Average difference threshold");
        register_param(m_dAverageSmoothing, "ave_smooth", "Average smoothing");
        register_param(m_dMinMaxSplitWeight, "split_weight", "Min/max split weight");
        register_param(m_dMinMaxThresholdWeight, "threshold_weight", "Min/max threshold weight");

        m_pBuffer = new uint8_t[width * height * 4];
        memset(m_pBuffer, 0x00, width * height * 4);
        
/*        m_barSize = 0.1;

        // Create the lookup table (see update()). Use std::vector instead of an array here, since:
        // http://stackoverflow.com/questions/381621/using-arrays-or-stdvectors-in-c-whats-the-performance-gap
        // Using std::numeric_limits just for fun here.
        lookupTable = std::vector<uint8_t>(std::numeric_limits<uint8_t>::max()+1, 0);
        std::cout << lookupTable.size() << " elements in the lookup table." << std::endl;

        // Calculate the entries in the lookup table. Applied on the R value in this example.
        // If the R value of an input pixel is r, then the output R value will be mapped to lookupTable[r].
        // We'll use a calculation that looks expensive here.
        float f, h;
        int tempVal;
        for (int i = 0; i < lookupTable.size(); i++) {
            f = i / (float)std::numeric_limits<uint8_t>::max(); // Normalize to [0,1]
            h = f/5;
            f = 2*f - 1; // Stretch to [-1,1]
            f = pow(f, 2); // Parabola
            f = f*.8 + h; // Modification to the parabola

            // Since we might get negative values above, directly putting f into an unsigned data type
            // may lead to interesting effects. To avoid this, use a signed integer and clamp it to valid ranges.
            tempVal = f * std::numeric_limits<uint8_t>::max();
            if (tempVal < 0) {
                tempVal = 0;
            }
            if (tempVal > std::numeric_limits<uint8_t>::max()) {
                tempVal = std::numeric_limits<uint8_t>::max();
            }
            lookupTable[i] = tempVal;
        }

        // This is a second lookup table for addition of two uint8 numbers. The result of {0,255}+{0,255}
        // ranges in {0,511}. The usual way to clamp the values to the {0,255} range is to use if/else
        // (if (k > 255) { k = 255; }), however using a lookup table is slightly faster:
        // http://stackoverflow.com/questions/4783674/lookup-table-vs-if-else
        additionTable = std::vector<uint8_t>( 2*std::numeric_limits<uint8_t>::max() + 1, 0 );
        for (int i = 0; i < 2*std::numeric_limits<uint8_t>::max(); i++) {
            if (i <= std::numeric_limits<uint8_t>::max()) {
                additionTable[i] = i;
            } else {
                additionTable[i] = std::numeric_limits<uint8_t>::max();
            }
        }
 */
    }

    ~time0r()
    {
        if (m_pBufferMonochromeAve)
            delete [] m_pBufferMonochromeAve;
        if (m_pBufferHistory)
            delete [] m_pBufferHistory;
        if (m_pBufferAve)
            delete [] m_pBufferAve;
        
        if (m_pBuffer)
            delete [] m_pBuffer;
        if (m_pBufferMonochrome)
            delete [] m_pBufferMonochrome;
        if (m_pBufferMonochromeMin)
            delete [] m_pBufferMonochromeMin;
        if (m_pBufferMonochromeMax)
            delete [] m_pBufferMonochromeMax;
        if (m_pBufferMin)
            delete [] m_pBufferMin;
        if (m_pBufferMax)
            delete [] m_pBufferMax;
        
        if (m_fMonoMin)
            delete [] m_fMonoMin;
        if (m_fMonoMax)
            delete [] m_fMonoMax;
        if (m_fHistory)
            delete [] m_fHistory;
        if (m_fMin)
            delete [] m_fMin;
        if (m_fMax)
            delete [] m_fMax;
        if (m_fAve)
            delete [] m_fAve;
    }
    
    static inline double to_grey(uint8_t* _in)
    {
        float in[] = {_in[0], _in[1], _in[2]};
        return to_grey(in);
    }
    
    static inline double to_grey(float* _in)
    {
        return RGBToV(_in);
    }
    static inline double RGBToGrey(float* _in)
    {
        return
            0.30/*0.34*//*0.00*/ * _in[0] +
            0.59/*0.66*//*0.00*/ * _in[1] +
            0.11/*0.00*//*1.00*/ * _in[2];
    }
    static inline double RGBToGrey(uint8_t* _in)
    {
        float in[] = {_in[0], _in[1], _in[2]};
        return RGBToGrey(in);
    }
    
    static inline uint8_t clamp8(float f)
    {
        return (uint8_t)min(255, max(0, (int)rint(f)));
    }
    
    static inline uint8_t to_grey8(uint8_t* _in)
    {
        return clamp8(to_grey(_in));
    }
    static inline float RGBToV(const /*double*/float R, const /*double*/float G, const /*double*/float B)
    {
        if ((B > G) && (B > R))
            return B;
        else if(G > R)
            return G;
        return R;
    }
    static inline float RGBToV(uint8_t* p)
    {
        return RGBToV((float)p[0], (float)p[1], (float)p[2]);
    }
    static inline float RGBToV(float* p)
    {
        return RGBToV(p[0], p[1], p[2]);
    }
    // 0 <= S <= 1 and 0 <= V <= 1, 0 <= H < 360
    static inline void RGBToHSV(const /*double*/float R, const /*double*/float G, const /*double*/float B, /*double*/float& H, /*double*/float& S, /*double*/float& V)
    {
        if((B > G) && (B > R))
        {
            V = B;
            if(V != 0)
            {
                double min;
                if(R > G) min = G;
                else      min = R;
                const double delta = V - min;
                if(delta != 0)
                { S = (delta/V); H = 4 + (R - G) / delta; }
                else
                { S = 0;         H = 4 + (R - G); }
                H *=   60; if(H < 0) H += 360;
                /*if(!NORM) */V =  (V/255);
                //else      S *= (100);
            }
            else
            { S = 0; H = 0;}
        }
        else if(G > R)
        {
            V = G;
            if(V != 0)
            {
                double min;
                if(R > B) min = B;
                else      min = R;
                const double delta = V - min;
                if(delta != 0)
                { S = (delta/V); H = 2 + (B - R) / delta; }
                else
                { S = 0;         H = 2 + (B - R); }
                H *=   60; if(H < 0) H += 360;
                /*if(!NORM) */V =  (V/255);
                //else      S *= (100);
            }
            else
            { S = 0; H = 0;}
        }
        else
        {
            V = R;
            if(V != 0)
            {
                double min;
                if(G > B) min = B;
                else      min = G;
                const double delta = V - min;
                if(delta != 0)
                { S = (delta/V); H = (G - B) / delta; }
                else
                { S = 0;         H = (G - B); }
                H *=   60; if(H < 0) H += 360;
                /*if(!NORM) */V =  (V/255);
                //else      S *= (100);
            }
            else
            { S = 0; H = 0;}
        }
    }
    
    // HSV values = [0,1]
    static inline void HSVToRGB(float H, float S, float V, /*int*/float& R, /*int*/float& G, /*int*/float &B)
    {
        if (S == 0)
        {
            R = (int)(V * 255.0f);
            G = (int)(V * 255.0f);
            B = (int)(V * 255.0f);
        }
        else
        {
            float var_h = H * 6;
            if (var_h == 6)
                var_h = 0;	// H must be < 1
    
            float var_i = floor(var_h);
            float var_1 = V * (1 - S);
            float var_2 = V * (1 - S * (var_h - var_i));
            float var_3 = V * (1 - S * (1 - (var_h - var_i)));
    
            float var_r, var_g, var_b;
            if      ( var_i == 0 ) { var_r = V     ; var_g = var_3 ; var_b = var_1; }
            else if ( var_i == 1 ) { var_r = var_2 ; var_g = V     ; var_b = var_1; }
            else if ( var_i == 2 ) { var_r = var_1 ; var_g = V     ; var_b = var_3; }
            else if ( var_i == 3 ) { var_r = var_1 ; var_g = var_2 ; var_b = V    ; }
            else if ( var_i == 4 ) { var_r = var_3 ; var_g = var_1 ; var_b = V    ; }
            else                   { var_r = V     ; var_g = var_1 ; var_b = var_2; }
    
            R = /*(int)*/(var_r * 255.0f);	// RGB results [0, 255]
            G = /*(int)*/(var_g * 255.0f);
            B = /*(int)*/(var_b * 255.0f);
        }
    }

    virtual void update()   // in, out, time, width, height, size
    {
        ++m_iFrameCount;
        
        uint8_t* _in = (uint8_t*)in;
        uint8_t* pOutput = m_pBuffer;
        
        if (m_iFrameCount == 0)
        {
            m_op = (Operation)((int)m_fOp);
            /*char* strOp;
            get_param_value(&strOp, 0);
            m_paramOp = strOp;
            
            if (m_paramOp == "min")
            {
                m_op = OP_MIN;
            }
            else
                fprintf(stderr, "time0r: unknown op: %s\n", m_paramOp.c_str());
            */
            
            m_iAveragePeriod = (int)m_dAveragePeriod;
        }
        
        if (m_op == OP_NULL)
        {
            memcpy(m_pBuffer, _in, size * 4);
        }
        else if ((m_op == OP_MIN) || (m_op == OP_MAX))   // 1, 2
        {
            if (m_iFrameCount == 0)
            {
                memcpy(m_pBuffer, in, size * 4);
                
                m_pBufferMonochrome = new uint8_t[width * height];
                
                if (m_op == OP_MIN)
                    memset(m_pBufferMonochrome, 0xFF, size);
                else if (m_op == OP_MAX)
                    memset(m_pBufferMonochrome, 0x00, size);
            }
            else
            {
                uint8_t* _in = (uint8_t*)in;
                
                for (int i = 0; i < size; i++)
                {
                    /*double dGrey =
                        0.30 * (double)_in[i*4 + 0] +
                        0.59 * (double)_in[i*4 + 1] +
                        0.11 * (double)_in[i*4 + 2];*/
                    uint8_t chGrey = to_grey8(_in + i*4);//(uint8_t)min(255, max(0, (int)floor(dGrey)));
                    
                    if ((((m_op == OP_MIN)) && (chGrey < m_pBufferMonochrome[i])) ||
                        (((m_op == OP_MAX)) && (chGrey > m_pBufferMonochrome[i])))
                    {
                        memcpy(m_pBuffer + i*4, &in[i], 4);
                        m_pBufferMonochrome[i] = chGrey;
                    }
                    
                    /*m_pBuffer[i*4 + 0] = min(m_pBuffer[i*4 + 0], _in[i*4 + 0]);
                    m_pBuffer[i*4 + 1] = min(m_pBuffer[i*4 + 1], _in[i*4 + 1]);
                    m_pBuffer[i*4 + 2] = min(m_pBuffer[i*4 + 2], _in[i*4 + 2]);
                    m_pBuffer[i*4 + 3] = _in[i*4 + 3];*/
                }
            }
        }
        else if ((m_op == OP_MINMAX) || (m_op == OP_MINMAXAVE) || (m_op == OP_MINMAXAVEDEBUG))  // 4, 5, 8
        {
            if (m_iFrameCount == 0)
            {
                //memcpy(m_pBuffer, in, size * 4);
                
                //memset(m_pBufferMonochrome, 0xFF/2, size);
                
                //m_pBufferMin = new uint8_t[width * height * 4];
                //m_pBufferMax = new uint8_t[width * height * 4];
                m_fMin = new float[size * 3];
                m_fMax = new float[size * 3];
                
                //m_pBufferMonochromeMin = new uint8_t[width * height];
                //m_pBufferMonochromeMax = new uint8_t[width * height];
                m_fMonoMin = new float[size];
                m_fMonoMax = new float[size];
                
                for (int i = 0; i < size; ++i)
                {
                    //m_fMonoMin[i] = 255.0f;
                    //m_fMonoMax[i] = 0.0f;
                    m_fMonoMin[i] = m_fMonoMax[i] = to_grey(_in + i*4);
                }
                
                //memset(m_pBufferMonochromeMin, 0xFF, size);
                //memset(m_pBufferMonochromeMax, 0x00, size);
                
                /*for (int i = 0; i < size; i++)
                {
                    double dGrey = 0.3 * (double)_in[i*4 + 0] + 0.59 * (double)_in[i*4 + 1] + 0.11 * (double)_in[i*4 + 2];
                    uint8_t chGrey = (uint8_t)min(255, max(0, (int)floor(dGrey)));
                    
                    m_pBufferMonochrome[i] = chGrey;
                }
                memcpy(m_pBufferMonochromeMin, m_pBufferMonochrome, size);
                memcpy(m_pBufferMonochromeMax, m_pBufferMonochrome, size);*/
                
                if ((m_op == OP_MINMAXAVE) || (m_op == OP_MINMAXAVEDEBUG))
                {
                    //m_pBufferHistory = new uint8_t[size * 4];
                    //memcpy(m_pBufferHistory, in, size * 4);
                    
                    m_fHistory = new float[size * 3];
                    for (int i = 0; i < size; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                            m_fHistory[i*3 + j] = _in[i*4 + j];
                    }
                }
            }
            
            const float fAveWeight = /*0.1f*/m_dAveragePeriod;
            const float fLimitWeight = /*0.1f*/m_dAverageDiffThreshold;
            
            pOutput = (uint8_t*)out;
//fprintf(stderr, "%i \n", (int)m_iFrameCount);
            for (int i = 0; i < size; i++)
            {
                if (((m_op == OP_MINMAXAVE) || (m_op == OP_MINMAXAVEDEBUG)) && (m_iFrameCount > 0))
                {
                    for (int j = 0; j < 3; j++)
                    {
                        //m_pBufferHistory[i*4 + j] = clamp8(((float)m_pBufferHistory[i*4 + j] * fAveWeight) + ((float)_in[i*4 + j] * (1.0f - fAveWeight)));
                        m_fHistory[i*3 + j] = (m_fHistory[i*3 + j] * fAveWeight) + ((float)_in[i*4 + j] * (1.0f - fAveWeight));
                    }
                }
                
                float fGrey = to_grey(_in + i*4);
                bool bCopy = false;
                //float r=0, g=0, b=0;
                //float fRGB[3] = {0.f, 0.f, 0.f};
                
                ///////////////////////////////////////////
                
                if (fGrey <= /*m_pBufferMonochromeMin[i]*/m_fMonoMin[i])
                {
                    //memcpy(m_pBuffer + i*4, &in[i], 4);
                    
                    //m_pBufferMonochromeMin[i] = chGrey;
                    m_fMonoMin[i] = fGrey;
                    
                    //if (chGrey < (0xFF - m_pBufferMonochromeMax[i]))
                    //    bCopy = true;
                    
                    //m_pBuffer[i*4 + 0] = chGrey;
                    
                    //memcpy(m_pBufferMin + i*4, _in + i*4, 4);
                    for (int j = 0; j < 3; ++j)
                        m_fMin[i*3 + j] = _in[i*4 + j];
                }
                
                //r += m_pBufferMin[i*4 + 0];
                //g += m_pBufferMin[i*4 + 1];
                //b += m_pBufferMin[i*4 + 2];
                
                ///////////////////////////////////////////
                
                if (fGrey >= /*m_pBufferMonochromeMax[i]*/m_fMonoMax[i])
                {
                    //memcpy(m_pBuffer + i*4, &in[i], 4);
                    
                    //m_pBufferMonochromeMax[i] = chGrey;
                    m_fMonoMax[i] = fGrey;
                    
                    //if ((0xFF - chGrey) < m_pBufferMonochromeMin[i])
                    //    bCopy = true;
                    
                    //m_pBuffer[i*4 + 2] = chGrey;
                    
                    //memcpy(m_pBufferMax + i*4, _in + i*4, 4);
                    for (int j = 0; j < 3; ++j)
                        m_fMax[i*3 + j] = _in[i*4 + j];
                }
                
                //r += m_pBufferMax[i*4 + 0];
                //g += m_pBufferMax[i*4 + 1];
                //b += m_pBufferMax[i*4 + 2];
                
                ///////////////////////////////////////////
                
                /*m_pBuffer*///pOutput[i*4 + 0] = clamp8(r / 2.0);
                /*m_pBuffer*///pOutput[i*4 + 1] = clamp8(g / 2.0);
                /*m_pBuffer*///pOutput[i*4 + 2] = clamp8(b / 2.0);
                
                float fAveGrey, fMinGrey, fMaxGrey;
                
                float fA[3], fB[3];
                if ((m_op == OP_MINMAXAVE) || (m_op == OP_MINMAXAVEDEBUG))
                {
                    fAveGrey = to_grey(m_fHistory + i*3);
                    
                    fMinGrey = to_grey(m_fMin + i*3);
                    fMaxGrey = to_grey(m_fMax + i*3);
                    
                    //if (abs(fMaxGrey - fAveGrey) >= abs(fAveGrey - fMinGrey))
                    if (abs(fMaxGrey - fAveGrey)*(m_dMinMaxThresholdWeight) >= abs(fAveGrey - fMinGrey)*(1.0f - m_dMinMaxThresholdWeight))
                    {
                        memcpy(fA, m_fMax + i*3, sizeof(float) * 3);
                        memcpy(fB, m_fMin + i*3, sizeof(float) * 3);
                    }
                    else
                    {
                        memcpy(fA, m_fMin + i*3, sizeof(float) * 3);
                        memcpy(fB, m_fMax + i*3, sizeof(float) * 3);
                    }
                }
                
                if (m_op != OP_MINMAXAVEDEBUG)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        if (m_op == OP_MINMAX)
                        {
                            pOutput[i*4 + j] = clamp8((m_fMin[i*3 + j] + m_fMax[i*3 + j]) / 2.0);
                        }
                        else if (m_op == OP_MINMAXAVE)
                        {
                            pOutput[i*4 + j] = clamp8((fA[j] * m_dMinMaxSplitWeight) + (fB[j] * (1.0f - m_dMinMaxSplitWeight)));
                        }
                        else if (m_op == OP_MINMAXAVEDEBUG)
                        {
                            //pOutput[i*4 + j] = clamp8(fMaxGrey - fMinGrey);
                            //pOutput[i*4 + j] = clamp8(m_fMonoMax[i] - m_fMonoMin[i]);
                            pOutput[i*4 + 0] = clamp8(m_fMonoMin[i]);
                            pOutput[i*4 + 1] = 0;
                            pOutput[i*4 + 2] = clamp8(m_fMonoMax[i]);
                            //pOutput[i*4 + 2] = clamp8(m_fMonoMax[i] - m_fMonoMin[i]);
                        }
                    }
                }
                else
                {
                    //pOutput[i*4 + 0] = clamp8(m_fMax - m_fMin);
                    //pOutput[i*4 + 1] = clamp8(m_fMax - m_fMin);
                    //pOutput[i*4 + 2] = clamp8(fAveGrey);
                    
                    pOutput[i*4 + 0] = clamp8(m_fMonoMin[i]);
                    pOutput[i*4 + 1] = clamp8(to_grey(m_fMax + i*3));
                    pOutput[i*4 + 2] = clamp8(m_fMonoMax[i]);
                }
                
                if (bCopy)
                {
                    //memcpy(m_pBuffer + i*4, &in[i], 4);
                }
                
                if ((m_op == OP_MINMAXAVE) || (m_op == OP_MINMAXAVEDEBUG))
                {
                    //uint8_t chAveGrey = to_grey(m_pBufferHistory + i*4);
                    
                    //m_pBufferMonochromeMin[i] = clamp8(((float)m_pBufferMonochromeMin[i] * (1.0f - fLimitWeight)) + ((float)chAveGrey * fLimitWeight));
                    //m_pBufferMonochromeMax[i] = clamp8(((float)m_pBufferMonochromeMax[i] * (1.0f - fLimitWeight)) + ((float)chAveGrey * fLimitWeight));
                    float fMin = (m_fMonoMin[i] * (1.0f - fLimitWeight)) + (fAveGrey * fLimitWeight);
                    if (fMin > m_fMonoMin[i])
                        m_fMonoMin[i] = fMin;
                    //else
                    {
                        for (int j = 0; j < 3; j++)
                            m_fMin[i*3 + j] = (m_fMin[i*3 + j] * (1.0 - m_dAverageSmoothing)) + (m_fHistory[i*3 + j] * m_dAverageSmoothing);
                    }
                    
                    float fMax = (m_fMonoMax[i] * (1.0f - fLimitWeight)) + (fAveGrey * fLimitWeight);
                    if (fMax < m_fMonoMax[i])
                        m_fMonoMax[i] = fMax;
                    //else
                    {
                        for (int j = 0; j < 3; j++)
                            m_fMax[i*3 + j] = (m_fMax[i*3 + j] * (1.0 - m_dAverageSmoothing)) + (m_fHistory[i*3 + j] * m_dAverageSmoothing);
                    }
                }
            }
        }
        else if ((m_op == OP_AVETHRESH) || (m_op == OP_AVETHRESHDEBUG)) // 6, 7
        {
            if (m_iFrameCount == 0)
            {
                m_iAveragePeriod = max(1, m_iAveragePeriod);
                //m_pBufferMonochromeAve = new uint8_t[size * m_iAveragePeriod];  // FIXME: Making float[]
                m_fHistory = new float[size * m_iAveragePeriod];
                memset(m_fHistory, 0x00, sizeof(float) * size * m_iAveragePeriod);
                m_fAve = new float[size];
                memset(m_fAve, 0x00, sizeof(float) * size);
            }
            
            for (int i = 0; i < size; i++)
            {
                /*double dGrey =
                    0.30 * (double)_in[i*4 + 0] +
                    0.59 * (double)_in[i*4 + 1] +
                    0.11 * (double)_in[i*4 + 2];*/
                //uint8_t chGrey = (uint8_t)min(255, max(0, (int)floor(dGrey)));
                float fGrey = /*to_grey*/RGBToGrey(_in + i*4);
                float fOrigGrey = fGrey;
                
                if (m_iFrameCount - (m_iAveragePeriod-1) > 0)
                {
                    //uint8_t chGrey = to_grey(_in + i*4);
                    int iSkip = (int)m_dAverageSmoothing;   // Default = 0
                    
                    float f = -1/*m_fAve[i]*/;
                    for (int j = ((m_iAverageBufferIndex/* + 1*/) % m_iAveragePeriod); (f == -1) || (j != m_iAverageBufferIndex); j = ((j + 1) % m_iAveragePeriod))
                    {
                        if (f == -1)
                            f = 0;
                        //f += m_pBufferMonochromeAve[i*m_iAveragePeriod + j];
                        f += m_fHistory[i*m_iAveragePeriod + j];
                    }
                    f /= (float)(m_iAveragePeriod/* - 1*/);
                    
                    if (iSkip != 0)
                    {
                        fGrey = -1;//m_fAve[i] - m_fHistory[i*m_iAveragePeriod + m_iAverageBufferIndex] + fGrey;
                        int iCount = 1;
                        for (int j = m_iAverageBufferIndex; (fGrey == -1) || (j != ((m_iAverageBufferIndex + /*(m_iAveragePeriod-1)*/iSkip)%m_iAveragePeriod)); j = ((j + 1) % m_iAveragePeriod))
                        {
                            if (fGrey == -1)
                                fGrey = /*0*/fOrigGrey;
                            //fGrey += m_pBufferMonochromeAve[i*m_iAveragePeriod + j];
                            fGrey += m_fHistory[i*m_iAveragePeriod + j];
                            ++iCount;
                        }
                        fGrey /= (float)(/*m_iAveragePeriod*//* - 1*/iCount);
                    }
                    
                    if (m_op == OP_AVETHRESH)
                    {
                        if (abs(fGrey - f) > m_dAverageDiffThreshold)
                        {
                            memcpy(m_pBuffer + i*4, &in[i], 4);
                        }
                    }
                    else if (m_op == OP_AVETHRESHDEBUG)
                    {
                        m_pBuffer[i*4 + 2] = clamp8(fGrey);
                        float fDiff = fGrey - f;
                        if (/*abs*/(fDiff) >= m_dAverageDiffThreshold)
                        {
                            m_pBuffer[i*4 + 0] = clamp8(/*fGrey*/fDiff * 32);
                            m_pBuffer[i*4 + 1] = 0;
                        }
                        else if (fDiff < -m_dAverageDiffThreshold)
                        {
                            m_pBuffer[i*4 + 0] = 0;
                            m_pBuffer[i*4 + 1] = clamp8(/*fGrey*/-fDiff * 32);
                        }
                        else
                        {
                            m_pBuffer[i*4 + 0] = 0;
                            m_pBuffer[i*4 + 1] = 0;
                        }
                        //if (fDiff > 0)
                        //    m_pBuffer[i*4 + 0] = clamp8(fDiff);
                        //else
                        //    m_pBuffer[i*4 + 2] = clamp8(-fDiff);
                    }
                }
                
                //m_pBufferMonochromeAve[i*m_iAveragePeriod + m_iAverageBufferIndex] = chGrey;
                m_fAve[i] -= m_fHistory[i*m_iAveragePeriod + m_iAverageBufferIndex];
                m_fHistory[i*m_iAveragePeriod + m_iAverageBufferIndex] = fOrigGrey;
                m_fAve[i] += fOrigGrey;
            }
            
            m_iAverageBufferIndex = (m_iAverageBufferIndex + 1) % m_iAveragePeriod;
        }
        else if (m_op == OP_AVE)    // 3
        {
            if (m_iFrameCount == 0)
            {
                m_iAveragePeriod = max(m_iAveragePeriod, 1);
                m_pBufferAve = new uint32_t[size * 3];
                memset(m_pBufferAve, 0x00, size * sizeof(uint32_t) * 3);
                m_pBufferHistory = new uint8_t[size * 4 * m_iAveragePeriod];
                memset(m_pBufferHistory, 0x00, size * 4 * m_iAveragePeriod);
            }
            
            pOutput = (uint8_t*)out;
            
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    uint8_t ch = m_pBufferAve[i*3 + j] / m_iAveragePeriod;  // Not bothering with floating-point
                    /*m_pBuffer*/pOutput[i*4 + j] = ch;
                    
                    m_pBufferAve[i*3 + j] -= m_pBufferHistory[i*4 + m_iAverageBufferIndex*size*4 + j];
                    m_pBufferAve[i*3 + j] += _in[i*4 + j];
                }
                
                //m_pBuffer[i*4 + 3] = ;
            }
            
            memcpy(m_pBufferHistory + size * 4 * m_iAverageBufferIndex, in, size * 4);
            
            m_iAverageBufferIndex = (m_iAverageBufferIndex + 1) % m_iAveragePeriod;
        }
        else if (m_op == OP_HSV)    // 9, V >= fGrey always
        {
            pOutput = (uint8_t*)out;
            
            for (int i = 0; i < size; ++i)
            {
                float fGrey = RGBToGrey(_in + i*4);
                //float V = RGBToV(_in + i*4);
                float H, S, V;
                RGBToHSV(_in[i*4 + 0], _in[i*4 + 1], _in[i*4 + 2], H, S, V);
                V *= 255;
                S *= 255;
                H *= (255.0 / 360.0);
                //for (int j = 0; j < 3; j++)
                //    pOutput[i*4 + j] = clamp8(abs(fGrey-V) * 2.0);
                /*float fDiff = fGrey - V;
                pOutput[i*4 + 1] = 0;
                if (fDiff < 0)
                {
                    pOutput[i*4 + 0] = clamp8(-fDiff * 8);
                    pOutput[i*4 + 2] = 0;
                }
                else
                {
                    pOutput[i*4 + 0] = 0;
                    pOutput[i*4 + 2] = clamp8(fDiff * 16);
                }*/
                pOutput[i*4 + 0] = clamp8(V);
                pOutput[i*4 + 1] = 0;//clamp8(H);   // Compressed even more
                pOutput[i*4 + 2] = clamp8(S * 2);   // Compressed quite heavily
            }
        }
        
        if ((uint8_t*)out != pOutput)
            memcpy(out, pOutput, size * 4);
        
        // Just copy input to output.
        // This is useful if ony few changes are made to the output.
        // If the whole image is processed, this makes no sense!
        //std::copy(in, in + width*height, out);

        // Fill the given amount of the top part with black.
        //std::fill(&out[0], &out[(int) (m_barSize*width*height)], 0);

        // Performance is important!
        // One way to improve the performance is to use a LOOKUP TABLE
        // instead of doing the same calculation several times. For example, the
        // SOP/Sat effect uses pow(), division, and multiplication for the Power parameter.
        // Since this only depends on the R/G/B value and not on previous frames, the target value
        // can be pre-computed; when applying the filter, only thing left to do is reading the value
        // in the lookup table.

        // This parameter allows to do simple benchmarking: Rendering a video with uint8_t pointers and with uint32_t pointers.
        // (Don't forget to substract the rendering time without this effect applied to avoid counting
        // encoding and decoding as well!)
/*        if (m_pointerMethod == 0) {
            uint8_t *in_pointer = (uint8_t *) in;
            uint8_t *out_pointer = (uint8_t *) out;
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {

                    // Apply the parabola to the R channel with a single lookup
                    *out_pointer++ = lookupTable[*in_pointer++];

                    // Add g+b and clamp with the second lookup table
                    *out_pointer = additionTable[*in_pointer + *(in_pointer+1)];
                    *out_pointer++;  *in_pointer++;

                    // Copy the other channels
                    *out_pointer++ = *in_pointer++;
                    *out_pointer++ = *in_pointer++;
                }
            }
        } else {
            // This method takes only 80% of the time if only processing the R channel,
            // and 90% of the time with additionally the G channel,
            // compared to the above solution using uint8_t pointers.
            for (int px = 0; px < width*height; px++) {
                out[px] =
                            // Parabola to Red channel
                            lookupTable[(in[px] & 0xFF)]
                            // g+b to Green channel
                          | (additionTable[((in[px] & 0xFF00) >> CHAR_BIT) + ((in[px] & 0xFF0000) >> 2*CHAR_BIT)] << CHAR_BIT)
                            // copy blue and alpha
                          | (in[px] & 0xFFFF0000);
            }
        }
 */
    }
};

frei0r::construct<time0r> plugin("time0r", "Time-based operations.",
                "Balint Seeber", 0,1,
                F0R_COLOR_MODEL_RGBA8888);
