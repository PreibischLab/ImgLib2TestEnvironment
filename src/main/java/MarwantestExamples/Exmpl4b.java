package MarwantestExamples;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.localneighborhood.old.LocalNeighborhood;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

public class Exmpl4b {

	public <T extends RealType<T> & NativeType<T>> Exmpl4b() throws ImgIOException {
		
		Img<T> img =  (Img<T>) new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new ArrayImgFactory<FloatType>());
		
		Gauss.inDoubleInPlace(new double[] {1, 1}, img);
		
		Img<BitType> display = findAndDisplayMinima(img, new ArrayImgFactory<BitType>(), new BitType());
		
		ImageJFunctions.show(img);
		ImageJFunctions.show(display);
		
	}
	public static <T extends Comparable<T>, U extends RealType<U>> Img<U> findAndDisplayMinima(RandomAccessibleInterval<T> source, ImgFactory<U> imageFactory, U outputeType) {
		Img<U> output = imageFactory.create(source, outputeType);
		
		Interval interval = Intervals.expand(source, -1);
		source = Views.interval(source, interval);
		
		final Cursor<T> center = Views.iterable(source).cursor();
		
		final RectangleShape shape = new RectangleShape(1, true);
		
		for (final Neighborhood<T> localNeighborhood : shape.neighborhoods(source)) {
			final T centerValue = center.next();
			
			boolean isMinimum = true;
			
			for (final T value : localNeighborhood) {
				if (centerValue.compareTo(value)>=0) {
					isMinimum = false;
					break;
				}
			}
			if(isMinimum) {
				HyperSphere<U> hyperSphere = new HyperSphere<U>(output, center, 1);
				for(U value:hyperSphere)
					value.setOne();
			}
		}
		
		return output;
	}
	public static void main(String[] args) throws ImgIOException {
		new ImageJ();
		new Exmpl4b();
	}
}
