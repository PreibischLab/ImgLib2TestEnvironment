package klim;

import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;

public class CentralCurveDetection {

	// detects the start and end points
	public static <T extends RealType<T>> void detectEnds(PointSampleList<T> endPoints, T intensity){
		// TODO: automate search
		// start 42 846
		// end 985 485
		
		endPoints.add(new Point(new long[]{42, 846}), intensity);
		endPoints.add(new Point(new long[]{985, 485}), intensity);
	}
	
	// detects central line 
	public static <T extends RealType<T>> void detectCentralCurve(RandomAccessibleInterval<T> img, PointSampleList<T> worm, PointSampleList<T> endPoints, RandomAccessible<T> out){
		RandomAccess<T> rImg = img.randomAccess();
		RandomAccess<T> rOut = out.randomAccess();
				
	}
	
	public static <T extends RealType<T>>void main(String[] args){
		
	}
	
}