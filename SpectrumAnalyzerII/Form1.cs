using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Net;
using SigProc;

/*************************************************************************/
/* This test application creates a power-of-two-sized buffer, fills it   */
/* with some generic signal data, runs an FFT on it, calculates          */
/* frequency-magnitude pairs and then plots them. Because there is no    */
/* imaginary input to the FFT, N/2 to N of the output array contains a   */
/* mirror of the spectra, which is why I only display points 0 - N/2.    */
/* Frequency of the bins is equal to (binIndex*samplingFrequency)/(N/2). */
/* It would be fairly simple to incrementally load audio data from a     */
/* file into the buffer, at which point it would make more sense to      */
/* display the magnitudes of the resulting FFT calculations on the       */
/* buffers on a waterfall plot instead of a simple graph.         - Eric */
/*************************************************************************/

namespace SpectrumAnalyzerII
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void Redraw(PaintEventArgs e) {

            Random rnd = new Random();

            // --- create signal buffer and vars --- //

            int bufMag = 10;
            int bufSize = 1 << bufMag;

            int samplingFrequency = 44100;

            double[] oscFreqs = { 10500, 15500 }; // frequencies of sine wave generators

            var sawtoothFreq = 1000; // frequency of sawtooth wave (infinite harmonics of amplitude 1/harmonic)

            double[] realAry = new double[bufSize];
            double[] imagAry = new double[bufSize];
            double[] magAry = new double[bufSize];

            // --- fill signal buffer with waveforms at specified frequencies --- //

            for (int i = 0; i < realAry.Length; i++)
            {

                // --- generate sawtooth wave at specified frequency --- //

                realAry[i] += ((double)(i % (samplingFrequency / sawtoothFreq))) / ((double)(samplingFrequency / sawtoothFreq));

                for (int j = 0; j < oscFreqs.Length; j++)
                {

                    // --- sine waves --- //

                    realAry[i] += Math.Sin(((double)i) * 2d * Math.PI * oscFreqs[j] / ((double)samplingFrequency)) / oscFreqs.Length;

                    imagAry[i] = 0;

                }

            }

            // --- create new FFT instance & initialize it --- //

            FFT myFFT = new FFT();

            myFFT.init((uint)bufMag);

            // --- analyze signal --- //

            myFFT.run(realAry, imagAry, false);

            // --- extract frequency magnitude information from sin/cos (real/imaginary) pairs --- //
            // --- (we don't care about phase information here, but if we did we could get it) --- //

            for (int i = 0; i < realAry.Length; i++)
            {

                magAry[i] = Math.Sqrt((realAry[i] * realAry[i]) + (imagAry[i] * imagAry[i]));

            }

            double magAryMax = magAry.Max();

            // --- normalize the array --- //
            // --- (max val is now 1)  --- //

            for (int i = 0; i < magAry.Length; i++)
            {

                magAry[i] = magAry[i] / magAryMax;

            }

            // --- //

            // Get the graphics object 
            Graphics gfx = e.Graphics;
            // Create a new pen that we shall use for drawing the line 
            Pen myPen = new Pen(Color.FromArgb(180, 255, 180));
            myPen.Width = 1;
            SolidBrush myBrush = new SolidBrush(Color.FromArgb(0, 0, 0));

            // --- draw a graph (frequency vs. magnitude) --- //

            int w = ClientRectangle.Width;
            int h = ClientRectangle.Height;

            gfx.FillRectangle(myBrush, Rectangle.FromLTRB(0, 0, w, h));

            for (int i = 0; i < (magAry.Length / 2) - 1; i++)
            {

                int xpos1, ypos1, xpos2, ypos2;

                xpos1 = w * i / (magAry.Length / 2);
                ypos1 = (int)(magAry[i] * ((double)h));

                xpos2 = w * (i + 1) / (magAry.Length / 2);
                ypos2 = (int)(magAry[i + 1] * ((double)h));

                ypos1 = h - ypos1;
                ypos2 = h - ypos2;

                gfx.DrawLine(myPen, xpos1, ypos1, xpos2, ypos2);

            }

        }

        private void Form1_Paint(object sender, PaintEventArgs e)
        {

            Redraw(e);

        } 
    }
}

// --- free, fast and small FFT lib --- //

namespace SigProc
{
    /**
     * Performs an in-place complex FFT.
     *
     * Released under the MIT License
     *
     * Copyright (c) 2010 Gerald T. Beauregard
     *
     * Permission is hereby granted, free of charge, to any person obtaining a copy
     * of this software and associated documentation files (the "Software"), to
     * deal in the Software without restriction, including without limitation the
     * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
     * sell copies of the Software, and to permit persons to whom the Software is
     * furnished to do so, subject to the following conditions:
     *
     * The above copyright notice and this permission notice shall be included in
     * all copies or substantial portions of the Software.
     *
     * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
     * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
     * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
     * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
     * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
     * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
     * IN THE SOFTWARE.
     */
    public class FFT
    {
        // Element for linked list in which we store the
        // input/output data. We use a linked list because
        // for sequential access it's faster than array index.
        class FFTElement
        {
            public double re = 0.0;		// Real component
            public double im = 0.0;		// Imaginary component
            public FFTElement next;		// Next element in linked list
            public uint revTgt;			// Target position post bit-reversal
        }

        private uint m_logN = 0;		// log2 of FFT size
        private uint m_N = 0;			// FFT size
        private FFTElement[] m_X;		// Vector of linked list elements

        /**
         *
         */
        public FFT()
        {
        }

        /**
         * Initialize class to perform FFT of specified size.
         *
         * @param	logN	Log2 of FFT length. e.g. for 512 pt FFT, logN = 9.
         */
        public void init(
            uint logN)
        {
            m_logN = logN;
            m_N = (uint)(1 << (int)m_logN);

            // Allocate elements for linked list of complex numbers.
            m_X = new FFTElement[m_N];
            for (uint k = 0; k < m_N; k++)
                m_X[k] = new FFTElement();

            // Set up "next" pointers.
            for (uint k = 0; k < m_N - 1; k++)
                m_X[k].next = m_X[k + 1];

            // Specify target for bit reversal re-ordering.
            for (uint k = 0; k < m_N; k++)
                m_X[k].revTgt = BitReverse(k, logN);
        }

        /**
         * Performs in-place complex FFT.
         *
         * @param	xRe		Real part of input/output
         * @param	xIm		Imaginary part of input/output
         * @param	inverse	If true, do an inverse FFT
         */
        public void run(
            double[] xRe,
            double[] xIm,
            bool inverse = false)
        {
            uint numFlies = m_N >> 1;	// Number of butterflies per sub-FFT
            uint span = m_N >> 1;		// Width of the butterfly
            uint spacing = m_N;			// Distance between start of sub-FFTs
            uint wIndexStep = 1; 		// Increment for twiddle table index

            // Copy data into linked complex number objects
            // If it's an IFFT, we divide by N while we're at it
            FFTElement x = m_X[0];
            uint k = 0;
            double scale = inverse ? 1.0 / m_N : 1.0;
            while (x != null)
            {
                x.re = scale * xRe[k];
                x.im = scale * xIm[k];
                x = x.next;
                k++;
            }

            // For each stage of the FFT
            for (uint stage = 0; stage < m_logN; stage++)
            {
                // Compute a multiplier factor for the "twiddle factors".
                // The twiddle factors are complex unit vectors spaced at
                // regular angular intervals. The angle by which the twiddle
                // factor advances depends on the FFT stage. In many FFT
                // implementations the twiddle factors are cached, but because
                // array lookup is relatively slow in C#, it's just
                // as fast to compute them on the fly.
                double wAngleInc = wIndexStep * 2.0 * Math.PI / m_N;
                if (inverse == false)
                    wAngleInc *= -1;
                double wMulRe = Math.Cos(wAngleInc);
                double wMulIm = Math.Sin(wAngleInc);

                for (uint start = 0; start < m_N; start += spacing)
                {
                    FFTElement xTop = m_X[start];
                    FFTElement xBot = m_X[start + span];

                    double wRe = 1.0;
                    double wIm = 0.0;

                    // For each butterfly in this stage
                    for (uint flyCount = 0; flyCount < numFlies; ++flyCount)
                    {
                        // Get the top & bottom values
                        double xTopRe = xTop.re;
                        double xTopIm = xTop.im;
                        double xBotRe = xBot.re;
                        double xBotIm = xBot.im;

                        // Top branch of butterfly has addition
                        xTop.re = xTopRe + xBotRe;
                        xTop.im = xTopIm + xBotIm;

                        // Bottom branch of butterly has subtraction,
                        // followed by multiplication by twiddle factor
                        xBotRe = xTopRe - xBotRe;
                        xBotIm = xTopIm - xBotIm;
                        xBot.re = xBotRe * wRe - xBotIm * wIm;
                        xBot.im = xBotRe * wIm + xBotIm * wRe;

                        // Advance butterfly to next top & bottom positions
                        xTop = xTop.next;
                        xBot = xBot.next;

                        // Update the twiddle factor, via complex multiply
                        // by unit vector with the appropriate angle
                        // (wRe + j wIm) = (wRe + j wIm) x (wMulRe + j wMulIm)
                        double tRe = wRe;
                        wRe = wRe * wMulRe - wIm * wMulIm;
                        wIm = tRe * wMulIm + wIm * wMulRe;
                    }
                }

                numFlies >>= 1; 	// Divide by 2 by right shift
                span >>= 1;
                spacing >>= 1;
                wIndexStep <<= 1;  	// Multiply by 2 by left shift
            }

            // The algorithm leaves the result in a scrambled order.
            // Unscramble while copying values from the complex
            // linked list elements back to the input/output vectors.
            x = m_X[0];
            while (x != null)
            {
                uint target = x.revTgt;
                xRe[target] = x.re;
                xIm[target] = x.im;
                x = x.next;
            }
        }

        /**
         * Do bit reversal of specified number of places of an int
         * For example, 1101 bit-reversed is 1011
         *
         * @param	x		Number to be bit-reverse.
         * @param	numBits	Number of bits in the number.
         */
        private uint BitReverse(
            uint x,
            uint numBits)
        {
            uint y = 0;
            for (uint i = 0; i < numBits; i++)
            {
                y <<= 1;
                y |= x & 0x0001;
                x >>= 1;
            }
            return y;
        }
    }
}