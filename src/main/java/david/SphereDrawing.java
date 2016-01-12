package david;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

public class SphereDrawing {
	
	public static <T extends NumericType<T>> void drawSphereCenter(Img<T> img, long size){
		Point center = new Point(img.numDimensions());
		for (int i = 0; i < img.numDimensions(); i++){
			center.setPosition(img.dimension(i)/2, i);
		}
		HyperSphere<T> hs = new HyperSphere<T>(img, center, size);
		Cursor<T> c = hs.cursor();
		
		
		T value = img.firstElement().createVariable();
		value.setZero();
		T myOne = img.firstElement().createVariable();
		myOne.setOne();
		
		while (c.hasNext()){
			c.fwd();
			c.get().set(value);
			value.add(myOne);
		}
		
		
	}
	
	public static void main(String[] args) {
		Img<FloatType> img = ArrayImgs.floats(100, 100, 100);
		new ImageJ();
		drawSphereCenter(img, 20);
		ImageJFunctions.show(img);
	}

}
