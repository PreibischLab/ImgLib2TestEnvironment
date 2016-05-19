package klim;

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

public class Ex4a {
	
	public Ex4a(){
		ImagePlusImg<UnsignedByteType, ?> img = new ImagePlusImgFactory<UnsignedByteType>().create(new long[] {256, 256, 256}, new UnsignedByteType());	
		drawSpheres(img, 0, 255);
		try{
			ImagePlus imp = img.getImagePlus();
			imp.show();
		}
		catch (ImgLibException e){
			System.out.println("Erroneous image");
			ImageJFunctions.show(img);
		}
	}
	
	public <T extends RealType<T>> void drawSpheres(
		final RandomAccessibleInterval<T> randomAccesible, 
		final double minVal, final double maxVal){
		
		int numDimensions = randomAccesible.numDimensions();
		Point center = new Point (randomAccesible.numDimensions());
		long minSize = randomAccesible.dimension(0);
		
		for (int d = 0; d < numDimensions; ++d){
			long size = randomAccesible.dimension(d);
			center.setPosition(size/2, d);
			minSize = Math.min(minSize, size);
		}
		
		int maxRadius = 5; 
		long radiusLargeSphere = minSize/2 - maxRadius - 1; 
		
		Random rnd = new Random(System.currentTimeMillis() ); //
		
		HyperSphere<T> hyperSphere = new HyperSphere<T>(randomAccesible, center, radiusLargeSphere);
		HyperSphereCursor<T> cursor = hyperSphere.cursor();
		
		while (cursor.hasNext()){
			cursor.fwd();
			int radius = rnd.nextInt(maxRadius) + 1;
			HyperSphere<T> smallSphere = new HyperSphere<>(randomAccesible, cursor, radius);
			double randomValue = rnd.nextDouble();
			
			if (Math.round(randomValue*100)%Util.pow(4, numDimensions) == 0){
				randomValue = rnd.nextDouble()*(maxVal - minVal) + minVal;
				for (final T value : smallSphere)
					value.setReal(Math.max(randomValue, value.getRealDouble()));
			}
		}
		
	}
	
	public static void main(String[] args){
		new ImageJ();
		new Ex4a();
	}
}
