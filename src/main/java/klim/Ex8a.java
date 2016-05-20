package klim;

import java.util.Random;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.IterableRealInterval;
import net.imglib2.KDTree;
import net.imglib2.Point;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithm.gauss.Gauss;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.neighborsearch.NearestNeighborSearchInterpolatorFactory;
import net.imglib2.neighborsearch.NearestNeighborSearch;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Ex8a {
	public Ex8a(){
		FinalInterval interval = new FinalInterval(new long[] {375, 200});
		IterableRealInterval<FloatType> realInterval = createRandomPoints(interval, 250);
		NearestNeighborSearch<FloatType> search = 
				new NearestNeighborSearchOnKDTree<FloatType>(
						new KDTree<FloatType>(realInterval));
		RealRandomAccessible<FloatType> realRandomAccessible = 
				Views.interpolate(search, new NearestNeighborSearchInterpolatorFactory<FloatType>());
		RandomAccessible<FloatType> randomAccessible = Views.raster(realRandomAccessible);
		RandomAccessibleInterval<FloatType> view = Views.interval(randomAccessible, interval);
		ImageJFunctions.show(view);
		Img<FloatType> convolved = new ArrayImgFactory<FloatType>().create(interval, new FloatType());
		Gauss.inFloat(new double[] {3, 3}, view, interval, convolved, new Point(view.numDimensions()), convolved.factory());
		ImageJFunctions.show(convolved);
	}
	
	public static RealPointSampleList<FloatType> createRandomPoints(RealInterval interval, int numPoints){
		int numDimensions = interval.numDimensions();
		Random rnd = new Random(System.currentTimeMillis());
		RealPointSampleList<FloatType> elements = new RealPointSampleList<FloatType>(numDimensions);
		
		for (int i =0; i < numPoints; ++i){
			RealPoint point = new RealPoint(numDimensions);
			
			for (int d = 0; d < numDimensions; d++) {
				point.setPosition( rnd.nextDouble() *
						( interval.realMax( d ) - interval.realMin( d ) ) + interval.realMin( d ), d );
			}
			elements.add(point, new FloatType(rnd.nextFloat()));
		}
		
		return elements;
	}
	
	public static void main(String[] args){
		new ImageJ(); 
		new Ex8a();
	}
}
