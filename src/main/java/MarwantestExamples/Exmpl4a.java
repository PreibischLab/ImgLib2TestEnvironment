package MarwantestExamples;

import java.util.Random;

import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.Point;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.exception.ImgLibException;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.util.Util;

public class Exmpl4a {

	public Exmpl4a() {
		
		ImagePlusImg<UnsignedByteType, ?> img = new ImagePlusImgFactory<UnsignedByteType>().create(new long[] {256, 256,256}, new UnsignedByteType());
		drawSpheeres(img,0,255);
		
		try {
			ImagePlus imp = img.getImagePlus();
			imp.show();
		}catch(ImgLibException e) {
			System.out.println("this is image is not native");
			ImageJFunctions.show(img);
		}
	}
	private <T extends RealType<T>>void drawSpheeres(RandomAccessibleInterval<T> randomAccessible, double minValue, double maxValue) {
		int numDimensions = randomAccessible.numDimensions();
		Point center = new Point(randomAccessible.numDimensions());
		long minSize = randomAccessible.dimension(0);
		for (int d =0; d<numDimensions;++d) {
			long size = randomAccessible.dimension(d);
			center.setPosition(size/2,d);
		}
		
		int maxRadius = 5;
		long radiusLargeSphere = minSize /2 - maxRadius-1;
		Random rnd = new Random(System.currentTimeMillis());
		HyperSphere<T> hyperSphere = new HyperSphere<T>(randomAccessible, center, radiusLargeSphere);
		HyperSphereCursor<T> cursor = hyperSphere.cursor();
		while(cursor.hasNext()) {
			cursor.fwd();
			int radius = rnd.nextInt(maxRadius)+1;
			HyperSphere<T> smallSphere = new HyperSphere<T>(randomAccessible, center, radius);
			double randomValue = rnd.nextDouble();
			if(Math.round(randomValue*100)% Util.pow(4, numDimensions)==0) {
				randomValue = rnd.nextDouble()*(maxValue-minValue)+minValue;
				
				for (final T value : smallSphere)
					value.setReal(Math.max(randomValue, value.getRealDouble()));
			}
		
		}
		
		
	}
	public static void main(String[] args) {
		new ImageJ();
		new Exmpl4a();
		
	}
}
