/**
 * This ImageJ plugin creates a grayscale image of width and height 100 pixels. The pixel values are uniformally randomly distributed from 0 - 255.
 * 
 * Part of exercise 3.4 from "Principles of Digital Image Processing" by Wilhelm Burger and Mark J. Burge.
 * 
 * @author Isaac Woodard
 */

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;

public class My_Random_Image implements PlugIn {
    
    /**
     * Creates the image and displays it.
     * 
     * @param arg optional argument
     */
	public void run(String arg) {
        //TODO: Allow the user to specify the image size.
        int w = 100;
        int h = 100;
        
        //Make and fill the image.
        ImageProcessor randIp = new ByteProcessor(w, h);
        double temp;
        for (int v = 0; v < h; v++) {         //it is apparently more efficient to iterate row by row due to how images are stored in memory
            for (int u = 0; u < w; u++) {
                int value = (int) (256 * Math.random());
                randIp.putPixel(u, v, value);
            }
        }

        //Display the image and give it a title.
        String title = "Random Grayscale Image";
        ImagePlus randIm = new ImagePlus(title, randIp);
        randIm.show();
	}
  
}
