/**
 * Equalizes an image based on either the linear or square root cumulative histogram. Also displays
 * the histogram and cumulative histogram as separate images.
 * 
 * Based on Programs 3.3 and 4.2 from "Principles of Digital Image Processing" by Wilhelm Burger and Mark J. Burge.
 * 
 * Note: Square root equalization isn't behaving as expected.
 * 
 * @author Isaac Woodard
 */

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import java.lang.Math;
import java.awt.Label;

public class My_Equalization implements PlugInFilter {
    String title = null;

    /**
     * Tells ImageJ the plugin works on 8-bit grayscale images and stores the image's name.
     * 
     * @param arg optional arguments
     * @param im the image to be processed
     * @return works on 8-bit greyscale images
     */
    public int setup(String arg, ImagePlus im) {
        title = im.getTitle();
        return DOES_8G;
    }
    
    /**
     * Creates the histograms and equalizes the image.
     * 
     * @param ip the original image
     */
    public void run(ImageProcessor ip) {
        
        //User dialog options
        String[] choices = {"Linear","Square Root"};
        GenericDialog gd = new GenericDialog("My_Equalization Settings");
        gd.addChoice("Cumulative Histogram Type", choices, choices[0]);
        gd.showDialog();
        if (gd.wasCanceled()) { //Cancel plugin if dialog box is closed
            IJ.error("PlugIn canceled!");
            return;
        }
        String choice = gd.getNextChoice();

        //Make original histogram
        int[] h1 = ip.getHistogram();
        int max = arrayMax(h1); //find largest histogram value to normalize against
        drawHistogram(h1, max, "Histogram of " + title);

        //Make original cumulative histogram
        if (choice.equals(choices[1])) { //calculate square root of first entry
            for (int j = 1; j < h1.length; j++) {
                h1[j] = (int) (Math.sqrt(h1[j]) + 0.5); 
            }
        }
        int[] H = cumulativeHistogram(h1);
        drawHistogram(H, H[255], "Cumulative Histogram of " + title);

        //Equalize the image
        ImageProcessor newIp = equalize(ip, H);

        //Make equalized histogram
        int[] h2 = newIp.getHistogram();
        max = arrayMax(h2);
        drawHistogram(h2, max, "Equalized Histogram of " + title);

        //Make equalized cumulative histogram
        if (choice.equals(choices[1])) { //calculate square root of first entry
            for (int j = 1; j < h2.length; j++) {
                h2[j] = (int) (Math.sqrt(h2[j]) + 0.5); 
            }
        }
        int[] J = cumulativeHistogram(h2);
        drawHistogram(J, J[255], "Equalized Cumulative Histogram of " + title);
    }

    /**
     * Helper method to equalize the image. Uses the cumulative histogram H.
     * @param ip an ImageProcessor
     * @param H  cumulative histogram
     */
    ImageProcessor equalize(ImageProcessor ip, int H[]) {
        int M = H[255];
        int K = 256; //8-bit range
        for (int v = 0; v < h; v++) {
            for (int u = 0; u < w; u++) {
                int a = ip.get(u, v);
                int b = H[a] * (K-1) / M;
                ip.set(u, v, b);
            }
        }
        return ip;
    }

    /**
     * Helper method to draw a histogram on the given image using a max scaling factor. Displays the image once created.
     * @param H histogram or cumulative histogram
     * @param max max value of histogram
     * @param title name for histogram
     */
    void drawHistogram(int H[], int max, String title) {
        int w = 256;
        int h = 100;
        ImageProcessor histIp = new ByteProcessor(w, 100);
        histIp.setValue(255);
        histIp.fill();

        for(int i = 0; i < w; i++) {
            int scale = (int)( H[i] / ((double) max) * 90); //...casting abuse? factor of 90 instead of 100 to keep room at top of image.
            for(int j = 0; j < scale; j++) { //draw line from image bottom
                histIp.set(i, 99 - j, 0); 
            }
        }

        ImagePlus histIm = new ImagePlus(title, histIp);
        histIm.show();
    }

    /**
     * Helper method to take a histogram array and compute the cumulative histogram.
     * @param h histogram array
     * @return the cumulative histogram
     */
    int[] cumulativeHistogram(int[] h) {
        int[] H = h.clone();
        for (int j = 1; j < H.length; j++) {
            H[j] = H[j-1] + H[j];
        }
        return H;
    }

    /**
     * Helper method to find the max value of an integer array.
     * @param M array of integers
     * @return max value of M
     */
    int arrayMax(int[] M) {
        int max = 0;
        for (int j = 1; j < M.length; j++) {
            if (max < M[j]) {
               max = M[j]; 
            }
        }
        return max;
    }
}


