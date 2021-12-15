/**
 * An ImageJ plugin to filter all pixels of an image by padding the edges of the image with grey pixels.
 * 
 * Part of exercise 5.3 from "Principles of Digital Image Processing" by Wilhelm Burger and Mark J. Burge. Based on Program 5.3.
 * 
 * @author Isaac Woodard
 */

import ij.process.*;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class My_Wrap_Filter implements PlugInFilter {
    
    /**
     * Tells ImageJ the plugin works on 8-bit grayscale images.
     * 
     * @param arg optional arguments
     * @param im the image to be processed
     * @return works on 8-bit greyscale images
     */
    public int setup(String arg, ImagePlus im) {
        return DOES_8G;
    }   
    
    /**
     * Performs the filtering of the image.
     * 
     * @param orig the original image
     */
    public void run(ImageProcessor orig) {
        int M = orig.getWidth();
        int N = orig.getHeight();
        
        // filter matrix of size (2K + 1) × (2L + 1) (columns by rows)
        int[][] filter = {
            {0,0,1,1,1,0,0},
            {0,1,1,1,1,1,0},
            {1,1,1,1,1,1,1},
            {0,1,1,1,1,1,0},
            {0,0,1,1,1,0,0}
        };
        
        double s = 1.0/filterSum(filter);

        int K = filter[0].length/2;
        int L = filter.length/2;
        
        ImageProcessor copy = orig.duplicate();

        //Create padded image and copy orginal image in the center of the padding
        ImageProcessor padded = new ByteProcessor(M+2*K, N+2*L);
        padded.setValue(127); //middle intensity value (grey)
        padded.fill();
        for (int v = L; v <= N+L-1; v++) {
            for (int u = K; u <= M+K-1; u++) {
                int a = orig.getPixel(u-K, v-L);
                padded.putPixel(u, v, a);
            }
        }
        
        //Filtering process
        for (int v = 0; v < N; v++) {
            for (int u = 0; u < M; u++) {
                // compute filter result for position (u, v)
                int sum = 0;
                for (int j = -L; j <= L; j++) {
                    for (int i = -K; i <= K; i++) {
                        int p = padded.getPixel(u+K+i, v+L+j);
                        int c = filter[j+L][i+K];
                        sum = sum + c * p;
                    }
                }
                int q = (int) Math.round(s * sum);
                if (q < 0) q = 0;
                if (q > 255) q = 255;
                orig.putPixel(u, v, q);
            }
        }
    }

    /**
     * Helper method to return the sum of the entries of a 2x2 array of integers.
     * 
     * @param filter the filter matrix
     * @return the sum of the elements in the filter
     */
    int filterSum(int[][] filter) {
        int rows = filter.length; 
        int cols = filter[0].length;
        int sum = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sum = sum + Math.abs(filter[i][j]);
            }
        }
        return sum;
    }
}
