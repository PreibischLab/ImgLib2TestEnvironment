package ella;

import ij.ImageJ;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import net.imglib2.AbstractRealLocalizable;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;


public class Sphere {
	
	public static void main(String[] args) {
		new ImageJ();
		final Img<FloatType> img = ArrayImgs.floats(100, 100);
		
		HyperSphere<FloatType> mySphere = new HyperSphere <FloatType> (img, new Point(50,50) ,  30);
		
		HyperSphereCursor< FloatType > c = mySphere.cursor();
		c.fwd();
		
		int i = 0;
		while ( c.hasNext()) {
			c.fwd();
			c.get().set(i++);
		}
		
		
		ImageJFunctions.show(img);
		
	}
}
