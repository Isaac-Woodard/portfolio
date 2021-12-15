/**
 * An ImageJ plugin to reflect an image horizontally. Exercise 2.1 from "Principles of Digital Image Processing" by Wilhelm Burger and Mark J. Burge.
 * Successful tests with 1x1, 2x2, 3x3, 9x9 and 10x10 images.
 * 
 * @author Isaac Woodard
 */

 import ij.ImagePlus;
 import ij.plugin.filter.PlugInFilter;
 import ij.process.ImageProcessor;

 public class My_Mirror_Horizontal implements PlugInFilter {

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
     * Performs the mirroring for the image.
     * 
     * @param ip original image
     */
    public void run(ImageProcessor ip) {
        int w = ip.getWidth();
        int h = ip.getHeight();
        int pTop, pBot;

        //Integer division rules allow us to avoid defining separate odd and even cases.
        for (int v = 0; v < (h / 2); v++) {
            for (int u = 0; u < w; u++) {
                pTop = ip.getPixel(u, v);
                pBot = ip.getPixel(u, h-1-v);

                ip.putPixel(u, v, pBot);
                ip.putPixel(u, h-1-v, pTop);
            }
        }

        //iterate over top half of image
        //if image is odd middle line of pixels will be untouched
        //switch pixel values across axis
        //make sure to store pixel values before overwriting!
    }
}