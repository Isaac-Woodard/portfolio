/**
 * An ImageJ plugin to remove spots up to a radius of 5 pixels using the shrink-grow technique. The foreground (black pixels)
 * is shrunk until the spots are gone and then grown back to its original size.
 * 
 * Part of exercise 7.2 from "Principles of Digital Image Processing" by Wilhelm Burger and Mark J. Burge.
 * 
 * @author Isaac Woodard
 */

import ij.*;
import ij.process.*;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import java.lang.Math;
import ij.process.Blitter;
import java.awt.Label;
import ij.gui.*;

public class My_Morphological_Filter implements PlugInFilter {

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
     * Defines the structuring element and performs the erosion and dilation process.
     * 
     * @param orig the original image
     */
    public void run(ImageProcessor orig) {
        int M = orig.getWidth();
        int N = orig.getHeight();
        
        //dilation/erosion filter (a.k.a. structuring element)
        //TODO: let the user define the structuring element
        int[][] struc = {
            {0,1,0},
            {1,1,1},
            {0,1,0}
        };

        //User dialog options
        GenericDialog gd = new GenericDialog("My_Equalization Settings");
        gd.addNumericField("Number of Iterations", 1, 2);
        gd.showDialog();
        if (gd.wasCanceled()) { //Cancel plugin if dialog box is closed
            IJ.error("PlugIn canceled!");
            return;
        }
        double choice = gd.getNextNumber();

        //Shrink then grow the image
        int repeats = (int) (choice + 0.5);
        for (int i = 0; i < repeats; i++) {
            erode(orig, struc);
        }
        for (int i = 0; i < repeats; i++) {
            dilate(orig, struc);
        }
    }

    /**
     * Dilates the image with the given structuring element. Implementation and comments copied from textbook (pg 180).
     * @param I the image
     * @param H the structuring element
     */
    void dilate(ImageProcessor I, int[][] H){
        //assume that the hot spot of H is at its center (ic,jc):
        int ic = (H[0].length-1)/2;
        int jc = (H.length-1)/2;
        
        //create a temporary (empty) image:
        ImageProcessor tmp = I.createProcessor(I.getWidth(), I.getHeight());
        for (int j = 0; j < H.length; j++){
            for (int i = 0; i < H[j].length; i++){
                if (H[j][i] > 0) { // this pixel is set
                    //copy image into position (i-ic,j-jc):
                    tmp.copyBits(I, i-ic, j-jc, Blitter.MAX);
                }
            }
        }
        //copy the temporary result back to original image
        I.copyBits(tmp, 0, 0, Blitter.COPY); //The text wrote "tmp" as "np". Probably a typo.
    }

    /**
     * Erodes the image with the given structuring element. Implementation and comments copied from textbook (pg 181).
     * Dilating the inverted image is equivalent to eroding the original image. (Just as efficient?)
     * @param I the image
     * @param H the structuring element
     */
    void erode(ImageProcessor I, int[][] H) {
        I.invert(); //The text wrote "I" as "ip". Probably a typo.
        dilate(I, reflect(H));
        I.invert();
    }

    /**
     * Helper method to reflect a 2x2 array of integers both horizontally and vertically.
     * @param H a matrix of integers
     * @return the mirrored matrix
     */
    int[][] reflect(int[][] H) {
        int rows = H.length; 
        int cols = H[0].length;

        //print for debugging
        // for (int u = 0; u < cols; u++) {
        //     for (int v = 0; v < rows; v++) {
        //         System.out.println(H[u][v]);
        //     }
        // } 
        // System.out.println("separator");

        //mirror horizontally
        for (int v = 0; v < (rows / 2); v++) {
            for (int u = 0; u < cols; u++) {
                int a = H[v][u];
                int b = H[rows-1-v][u];

                H[v][u] = b;
                H[rows-1-v][u] = a;
            }
        }

        //mirror vertically
        for (int u = 0; u < (cols / 2); u++) {
            for (int v = 0; v < rows; v++) {
                int a = H[v][u];
                int b = H[v][cols-1-u];

                H[v][u] = b;
                H[v][cols-1-u] = a;
            }
        }
        
        //print for debugging
        // for (int u = 0; u < cols; u++) {
        //     for (int v = 0; v < rows; v++) {
        //         System.out.println(H[u][v]);
        //     }
        // } 

        return H;
    }
}
