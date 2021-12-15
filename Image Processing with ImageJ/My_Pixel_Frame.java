/**
 * An ImageJ plugin to paint a 10 pixel wide white frame within an 8-bit grayscale image. Exercise 2.3 from "Principles of Digital Image Processing" 
 * by Wilhelm Burger and Mark J. Burge.
 * 
 * Successful tests with a 1x1 and 100 x 100 image.
 * 
 * @author Isaac Woodard
 */

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class My_Pixel_Frame implements PlugInFilter {
    
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
     * Draws rectangles around the image edge to create the frame.
     * 
     * @param ip the original image
     */
    public void run(ImageProcessor ip) {
        int w = ip.getWidth();
        int h = ip.getHeight();
        int intensity = 255;

        //if height or width is less than 21 pixels all pixels will be painted to white
        //otherwise "normal" behavior (i.e. paint a frame)
        if (w < 21 || h < 21) {
            this.rectangle(ip, 0, 0, w, h, intensity);
        } else {
            this.rectangle(ip, 0, 0, w, 10, intensity); //top
            this.rectangle(ip, 0, h-10, w, 10, intensity); //bottom
            this.rectangle(ip, 0, 10, 10, h-20, intensity); //left
            this.rectangle(ip, w-10, 10, 10, h-20, intensity); //right
        }
    }

    /**
     * Paints a rectangle with a width and height starting at pixel 1. Color is grayscale with values from 0-255.
     * 
     * @warning Assumes the rectangle is in bounds.
     * @param ip the original image
     */
    void rectangle(ImageProcessor ip, int u, int v, int w, int h, int intensity) {
        for (int i = 0; i < w; i++) {
            int vtemp = v;
            for (int j = 0; j < h; j++) {
                ip.putPixel(u, vtemp, intensity);
                vtemp++;
            }
            u++;
        }
    }
}
