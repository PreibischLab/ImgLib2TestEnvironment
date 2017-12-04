package marwan;

import java.util.Vector;

import net.imglib2.Cursor;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Helper {

	public static Boolean log ;
	public static Vector<RandomAccessibleInterval<FloatType>> splitImage(RandomAccessibleInterval<FloatType> input,
			int columns, int rows) {

		// TODO make it dynamic with 3d
		// input.numDimensions();
		long totalWidth = input.dimension(0);
		long totalHeight = input.dimension(1);

		// TODO get the rest / make dynamic 
		long width = totalWidth / columns-1;
		long height = totalHeight / rows-1;

		long lostWidth = totalWidth % columns;
		long lostHeight = totalHeight % rows;
		if (log) System.out.println("width: " + width +" | height: " + height +"columns: " + columns+ " | rows:" + rows);
		if (log) System.out.println("totalWidth: " + totalWidth +" | totalHeight: " + totalHeight +"lostWidth: " + lostWidth+ " | lostHeight:" + lostHeight);
		// TODO overlap
		// double overlap = 0.02;
		
		Vector<RandomAccessibleInterval<FloatType>> views = new Vector<RandomAccessibleInterval<FloatType>>();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				RandomAccessibleInterval<FloatType> view = Views.interval(input, new long[] { j * width, i * height },
						new long[] { (j + 1) * width, (i + 1) * height });
				if (log) System.out.println("start x: " + j * width +	"  y: " + i * height +"| arrive x: " + (j + 1) * width+ " | y:" + (i + 1) * height);
				views.add(view);
			}
		}
		if(log) System.out.println("split finished");
		return views;
	}

	public static float computeMiLocation(final Img<FloatType> input, int startPoint, int arrivePoint) {
		final Cursor<FloatType> cursor = input.cursor();

		cursor.jumpFwd(startPoint);
		float min = cursor.next().getRealFloat();
		for (long j = 0; j < arrivePoint; ++j) {
			final float v = cursor.next().getRealFloat();

			min = Math.min(min, v);
		}
		return min;
	}
}
