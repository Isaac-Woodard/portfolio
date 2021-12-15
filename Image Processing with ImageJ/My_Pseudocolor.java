/**
 * An ImageJ plugin to create a indexed color image from a grayscale image.
 * Part of exercise 8.5 from "Principles of Digital Image Processing" by Wilhelm Burger and Mark J. Burge.
 * 
 * @author Isaac Woodard
 */

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ColorProcessor;

public class My_Pseudocolor implements PlugInFilter {
    
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
     * Makes the indexed color image and displays it.
     * 
     * @param ip the original image
     */
    public void run(ImageProcessor ip) {
        int w = ip.getWidth();
        int h = ip.getHeight();

        //create color table
        //TODO: allow the user to specify the colors
        int[] color1 = {139, 0, 0};
        int[] color2 = {255, 255, 0};
        int[] color3 = {255, 255, 255};
        int[][] palette = makeColorTable(color1, color2, color3, 256);
        
        //iterate over dimensions replacing one pixel at a time
        ColorProcessor cp = new ColorProcessor(w, h);
        for (int v = 0; v < h; v++) {
            for (int u = 0; u < w; u++) {
                int k = ip.getPixel(u, v);
                int R = palette[0][k];
                int G = palette[1][k];
                int B = palette[2][k];               
                int[] RGB = {R, G, B};
                cp.putPixel(u, v, RGB);
            }
        }

        //display new image
        ImagePlus cimg = new ImagePlus("Indexed Color Image", cp);
        cimg.show();
    }

    /**
     * Makes a color table with a number of entries equal to size. Interpolates between colors 1 and 2 for the first half
     * of the entries and between colors 2 and 3 for the second half. Each color should have three entries for RGB values.
     * Limitations: Size must be even. 
     * 
     * @warning Assumes monotonically increasing RGB values from color 1 to color 2 to color 3.
     * @param color1 array of RGB values for color 1
     * @param color2 array of RGB values for color 2
     * @param color3 array of RGB values for color 3
     * @param size number of entries in the color table
     * @return the color table
     */
    int[][] makeColorTable(int[] color1, int[] color2, int[] color3, int size) {
        int[][] palette = new int[3][size];
        int[] dif12 = new int[3];
        int[] dif23 = new int[3];

        //calculate color differences
        for (int i = 0; i < 3; i++) {
            dif12[i] = color2[i] - color1[i];
            dif23[i] = color3[i] - color2[i];
        }        

        //first half of color entries
        for (int k = 0; k < size/2; k++) {
            for (int i = 0; i < 3; i++) {
                palette[i][k] = (int) (color1[i] + dif12[i] * ((double) k / (size/2 - 1)));
            }
        }

        //second half of color entries
        for (int k = size/2 ; k < size; k++) {
            for (int i = 0; i < 3; i++) {
                palette[i][k] = (int) (color2[i] + dif23[i] * (((double) k - size/2) / (size/2 - 1)));
            }
        }
        
        return palette;
    }
}