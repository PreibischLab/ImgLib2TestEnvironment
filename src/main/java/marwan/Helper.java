package marwan;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class Helper {

	public enum Task {
		Gaus, OriginalSize, Collect;
	}

	public static Boolean log;
	public static int count, sigma;

	public static ArrayList<Portion> splitImage(RandomAccessible<FloatType> input, FinalInterval interval, int columns,
			int rows) {

		// TODO make it dynamic with 3d
		// input.numDimensions();
		long totalWidth = interval.dimension(0) - (Helper.sigma * 2);
		long totalHeight = interval.dimension(1) - (Helper.sigma * 2);

		// TODO get the rest / make dynamic
		long width = totalWidth / columns - 1;
		long height = totalHeight / rows - 1;

		long lostWidth = totalWidth % columns;
		long lostHeight = totalHeight % rows;

		ArrayList<Portion> portions = new ArrayList<Portion>();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				RandomAccessibleInterval<FloatType> view = Views.interval(input,
						new long[] { j * width - Helper.sigma, i * height - Helper.sigma },
						new long[] { (j + 1) * width + Helper.sigma, (i + 1) * height + Helper.sigma });
				Rectangle shape = new Rectangle(new Point((int) (j * width), (int) (i * height)),
						new Dimension((int) (width), (int) (height)));
				portions.add(new Portion(view, shape));
			}
		}
		return portions;
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

	public static RandomAccessible<FloatType> extend(Img<FloatType> image, int sigma2) {
		RandomAccessible<FloatType> infiniteImg = Views.extendMirrorSingle(image);

		long[] min = new long[image.numDimensions()];
		long[] max = new long[image.numDimensions()];
		for (int d = 0; d < image.numDimensions(); ++d) {
			min[d] = -Helper.sigma;
			max[d] = image.dimension(d) + Helper.sigma;
		}
		FinalInterval interval = new FinalInterval(min, max);
		return Views.interval(infiniteImg, interval);

	}

	public static long[] getDimensions(Img<FloatType> image) {
		long[] dimensions = new long[image.numDimensions()];
		for (int d = 0; d < image.numDimensions(); ++d) {
			dimensions[d] = image.dimension(d);
		}

		return dimensions;
	}

	public static ArrayList<Portion> getMainSize(ArrayList<Portion> portions) {
		ArrayList<Portion> newPortions = new ArrayList<Portion>();
		for (Portion portion : portions) {
			RandomAccessibleInterval<FloatType> view = Views.interval(portion.getView(),
					new long[] { Helper.sigma, Helper.sigma },
					new long[] { portion.getView().dimension(0) - Helper.sigma,
							portion.getView().dimension(0) - Helper.sigma });
			newPortions.add(new Portion(view, portion.getShape()));
		}
		return newPortions;
	}

	public static RandomAccessibleInterval<FloatType> targetPositon(RandomAccessibleInterval<FloatType> targetView,
			Rectangle shape) {
		RandomAccessibleInterval<FloatType> view = Views.interval(targetView,
				new long[] { (long) shape.getX(), (long) shape.getY() },
				new long[] { (long) shape.getX() + (long) shape.getWidth(),
						(long) shape.getY() + (long) shape.getHeight() });

		return view;
	}

	public static FinalInterval getFinalInterval(Img<FloatType> image) {
		long[] min = new long[image.numDimensions()];
		long[] max = new long[image.numDimensions()];
		for (int d = 0; d < image.numDimensions(); ++d) {
			min[d] = -Helper.sigma;
			max[d] = image.dimension(d) + Helper.sigma;
		}
		FinalInterval interval = new FinalInterval(min, max);
		// TODO Auto-generated method stub
		return interval;
	}
}
