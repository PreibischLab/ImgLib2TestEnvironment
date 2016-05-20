package klim;

import java.io.File;
import java.util.Random;

import ij.ImageJ;
import net.imglib2.Interval;
import net.imglib2.IterableRealInterval;
import net.imglib2.KDTree;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.neighborsearch.InverseDistanceWeightingInterpolatorFactory;
import net.imglib2.interpolation.neighborsearch.NearestNeighborSearchInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.neighborsearch.KNearestNeighborSearch;
import net.imglib2.neighborsearch.KNearestNeighborSearchOnKDTree;
import net.imglib2.neighborsearch.NearestNeighborSearch;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex8b {
	public Ex8b(){
		Img <FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/LENNA.JPG"));
		ImageJFunctions.show(img).setTitle("initial pictire");	
		RealRandomAccessible<FloatType> realRandomAccessible = Views.interpolate(Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<FloatType>());
		for (int numPoints = 2; numPoints <= 32768; numPoints *= 4) {
			ImageJFunctions.show(randomSampling(realRandomAccessible, img, numPoints), numPoints + " points (NN)");
			ImageJFunctions.show(randomSamplingKNearest(realRandomAccessible, img, numPoints), numPoints + " points (kNN)");	
		}
	}
	
	
	public <T extends Type<T>> RandomAccessibleInterval<T> randomSampling(RealRandomAccessible<T> input, Interval interval, int numPoints){
		IterableRealInterval <T> realInterval = sampleRandomPoints(input, interval, numPoints);
		NearestNeighborSearch<T> search = new NearestNeighborSearchOnKDTree<T>(new KDTree<T>(realInterval));
		RealRandomAccessible<T> realRandomAccessible = Views.interpolate(search, new NearestNeighborSearchInterpolatorFactory<T>());
		RandomAccessible<T> randomAccessible = Views.raster(realRandomAccessible);
		return Views.interval(randomAccessible, interval);
	}
	
	public <T extends RealType<T>> RandomAccessibleInterval<T> randomSamplingKNearest(
			RealRandomAccessible<T> input, Interval interval, int numPoints){
		IterableRealInterval<T> realInterval = sampleRandomPoints(input, interval, numPoints);
		KNearestNeighborSearch<T> search = new KNearestNeighborSearchOnKDTree<>(new KDTree<T>(realInterval), Math.min(20, (int) realInterval.size()));
		RealRandomAccessible<T> realRandomAccessible = Views.interpolate(search, new InverseDistanceWeightingInterpolatorFactory<T>());
		RandomAccessible<T> randomAccessible = Views.raster(realRandomAccessible);
		return Views.interval(randomAccessible, interval);
	}
	
	public static <T extends Type<T>> RealPointSampleList<T> sampleRandomPoints(RealRandomAccessible<T> input, RealInterval interval, int numPoints){
		int numDimensions = interval.numDimensions();
		Random rnd = new Random(System.currentTimeMillis());
		RealPointSampleList<T> elements = new RealPointSampleList<T>(numDimensions);
		RealRandomAccess<T> realRandomAccess = input.realRandomAccess();
		for (int i = 0; i < numPoints; i++) {
			RealPoint point = new RealPoint(numDimensions);
			for (int d = 0; d < numDimensions; d++) {
				point.setPosition( rnd.nextDouble() *
						( interval.realMax( d ) - interval.realMin( d ) ) + interval.realMin( d ), d );
			}
			realRandomAccess.setPosition(point);
			elements.add(point, realRandomAccess.get().copy());
		}
		return elements;
	}
	
	public static void main(String [] args){
		new ImageJ();
		new Ex8b();
	}
}
