package ella;

import ij.ImageJ;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import java.util.Random;


public class MiniMeSphere {
	
	public static void drawSpheres(Img<FloatType> img, Point center, int radius) {
		HyperSphere<FloatType> mamaSphere = new HyperSphere <FloatType> (img, center,  radius);
		
		HyperSphereCursor< FloatType > c = mamaSphere.cursor();
		
		Random rand = new Random();
		
		while ( c.hasNext()) {
			c.fwd();
			if (rand.nextDouble() < 0.002) {
				HyperSphere<FloatType> babySphere = new HyperSphere <FloatType> (img, c ,  Math.round(rand.nextFloat()*10));
				
				float intensity = rand.nextFloat();
				for (FloatType f : babySphere) {
					if (intensity > f.get()) {
						f.set(intensity);
					}
				}
			}
		}
		
		
	}
	
	public static void main(String[] args) {
		new ImageJ();
		final Img<FloatType> img = ArrayImgs.floats(500, 500, 500);
		
		drawSpheres(img, new Point (250,250, 250), 200);
		
		ImageJFunctions.show(img);
		
	}
}