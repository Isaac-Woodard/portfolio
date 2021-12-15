/**
 * An ImageJ plugin to filter for the edges in a greyscale image using the Sobel edge operator. Generates 
 * an image for the edge magnitude and an image for the edge orientation.
 * 
 * Part of exercise 6.2 from "Principles of Digital Image Processing" by Wilhelm Burger and Mark J. Burge.
 * 
 * @author Isaac Woodard
 */

import ij.process.*;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.lang.Math;

public class My_Edge_Filter implements PlugInFilter {
    
    /**
     * Tells ImageJ the plugin works on 8-bit greyscale images and doesn't modify
     * the original image.
     * 
     * @param arg optional arguments
     * @param im the image to be processed
     * @return works on 8-bit greyscale images and doesn't modify the original image
     */
    public int setup(String arg, ImagePlus im) {
        return DOES_8G + NO_CHANGES;
    }   
    
    /**
     * Defines the filters and creates the images for the edge magnitude and edge orientation.
     * 
     * @param orig the image to be processed
     */
    public void run(ImageProcessor orig) {
        int M = orig.getWidth();
        int N = orig.getHeight();
        
        //Horizontal edge operator filter (2K + 1) × (2L + 1) (columns by rows)
        int[][] hFilter = {
            {-1,0,1},
            {-2,0,2},
            {-1,0,1}
        };

        //Vertical edge operator filter (2K + 1) × (2L + 1) (columns by rows)
        int[][] vFilter = {
            {-1,-2,-1},
            { 0, 0, 0},
            { 1, 2, 1}
        };

        //Make images for derivates
        ImageProcessor Dx = orig.duplicate();
        ImageProcessor Dy = orig.duplicate();

        //Apply derivative filters
        filter(hFilter, Dx);
        filter(vFilter, Dy);

        //Debugging
        // ImagePlus a = new ImagePlus("dx", Dx);
        // a.show();
        // ImagePlus b = new ImagePlus("dy", Dy);
        // b.show();

        //Create edge magnitude image
        ImageProcessor edgeMag = new ByteProcessor(M, N);
        for (int v = 0; v < N; v++) {
            for (int u = 0; u < M; u++) {
                int x = Dx.get(u, v);
                int y = Dy.get(u, v);
                int m = (int) (Math.sqrt(x*x + y*y) + 0.5);
                edgeMag.putPixel(u, v, m);
            }
        }
        
        ImagePlus EdgeMag = new ImagePlus("Edge Magnitude", edgeMag);
        EdgeMag.show();

        //Create edge orientation image
        ImageProcessor edgeOr = new ByteProcessor(M,N);
        for (int v = 0; v < N; v++) {
            for (int u = 0; u < M; u++) {
                int x = Dx.get(u, v);
                int y = Dy.get(u, v);
                int r = (int) (((255/Math.PI) * Math.atan2(y, x) + Math.PI) + 0.5); //normalize pi to intensity range
                edgeOr.putPixel(u, v, r);
            }
        }
        
        ImagePlus EdgeOr = new ImagePlus("Edge Orientation", edgeOr);
        EdgeOr.show();        
    }

    /**
     * Helper method to perform the filtering on an image. Calls paddedImage().
     * 
     * @param filter the filter matrix
     * @param orig   the original image
     */
    void filter(int[][] filter, ImageProcessor orig) {
        double s = 1.0 / filterSum(filter);
        int K = filter[0].length/2;
        int L = filter.length/2;
        
        int M = orig.getWidth();
        int N = orig.getHeight();
        ImageProcessor padded = paddedImage(orig, L, K);

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

    /**
     * Helper method to create a padded image from the original image. Padding is
     * based on rows L in filter and columns K in filter.
     * 
     * @param orig the original image
     * @param L the number of rows in the filter
     * @param K the number of columns in the filter
     * @return the original image with padding
     */
    ImageProcessor paddedImage(ImageProcessor orig, int L, int K) {
        int M = orig.getWidth();
        int N = orig.getHeight();
        ImageProcessor padded = new ByteProcessor(M+2*K, N+2*L);
        padded.setValue(127); //middle intensity value (grey)
        padded.fill();
        for (int v = L; v <= N+L-1; v++) {
            for (int u = K; u <= M+K-1; u++) {
                int a = orig.getPixel(u-K, v-L);
                padded.putPixel(u, v, a);
            }
        }
        return padded;
    }
}
