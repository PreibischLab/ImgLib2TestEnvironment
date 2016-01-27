package varun;

import java.io.File;
import java.util.ArrayList;
import java.util.Random;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class SortList {

	public static <T extends RealType<T>> PointSampleList<T> getList(RandomAccessibleInterval<T> img) {

		final RandomAccessible<T> infinite = Views.extendZero(img);

		final int n = img.numDimensions();
		long min[] = new long[n];
		long max[] = new long[n];

		for (int d = 0; d < n; ++d) {

			min[d] = img.min(d);
			max[d] = img.max(d);

		}

		FinalInterval interval = new FinalInterval(min, max);

		final IterableInterval<T> imgav = Views.interval(infinite, interval);

		final Cursor<T> first = imgav.cursor();

		// A point sample list with coordinates declared and initialized.
		PointSampleList<T> parent = new PointSampleList<T>(n);

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(n);

			cord.setPosition(first);

			parent.add(cord, first.get().copy());

		}

		return parent;

	}

	public static <T extends RealType<T>> void listtoImage(PointSampleList<T> list,
			RandomAccessibleInterval<T> imageout) {

		final Cursor<T> testCursor = list.localizingCursor().copyCursor();
		final RandomAccess<T> imageCursor = imageout.randomAccess();

		while (testCursor.hasNext()) {

			testCursor.fwd();

			imageCursor.setPosition(testCursor);
			imageCursor.get().set(testCursor.get());

		}

	}

	public static <T extends RealType<T>> void split(PointSampleList<T> list, int direction) {

		int n = list.numDimensions();

		if (list.dimension(direction) < 2)
			return;

		else {

			int meanIndex;

			final Cursor<T> listCursor = list.localizingCursor().copyCursor();

			PointSampleList<T> childA = new PointSampleList<T>(n);

			PointSampleList<T> childB = new PointSampleList<T>(n);

			// meanIndex = (int) (list.min(direction)+ ((list.max(direction) -
			// list.min(direction)) / 2)+ list.dimension(direction) % 2 ) ;

			meanIndex = (int) (list.min(direction) + ((list.dimension(direction) ) / 2)+ list.dimension(direction) % 2);

			// Here I split the list along the median in one direction

			while (listCursor.hasNext()) {

				listCursor.fwd();

				Point splitPoint = new Point(n);

				splitPoint.setPosition(meanIndex, direction);

				Point cord = new Point(n);

				cord.setPosition(listCursor);
				if (listCursor.getLongPosition(direction) < splitPoint.getLongPosition(direction)) {

					childA.add(cord, listCursor.get().copy());

				//	 System.out.println("childA: "+listCursor.get());

				} else

				{

					childB.add(cord, listCursor.get().copy());
				//	 System.out.println("childB: "+listCursor.get());
				}
				
				
			}

			split(childA, direction);

			split(childB, direction);

			mergeList(list, childA, childB);

		}

	}

	public static <T extends RealType<T>> void mergeList(PointSampleList<T> list, PointSampleList<T> listA,
			PointSampleList<T> listB) {

		Cursor<T> cursorA = listA.localizingCursor().copyCursor();

		Cursor<T> cursorB = listB.localizingCursor().copyCursor();

		Cursor<T> cursor = list.localizingCursor().copyCursor();

		cursor.fwd();
		cursorA.fwd();
		cursorB.fwd();

		
		
		
		while ( cursorA.hasNext() && cursorB.hasNext()) {

			if (cursorA.get().compareTo(cursorB.get()) < 0) {

				cursor.get().set(cursorA.get().copy());

				cursorA.fwd();
				cursor.fwd();

			} else {

				cursor.get().set(cursorB.get().copy());

				cursorB.fwd();
				cursor.fwd();
			}

		}

		while (cursorA.hasNext()) {

			cursor.get().set(cursorA.get().copy());

			cursorA.fwd();
			cursor.fwd();

		}

		while (cursorB.hasNext()) {

			cursor.get().set(cursorB.get().copy());

			cursorB.fwd();
			cursor.fwd();
		}

	}

	public static void main(String[] args) {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));

		Img<FloatType> imgout = new CellImgFactory<FloatType>().create(img, img.firstElement());

		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

		IterableInterval<FloatType> view = Views.interval(img, new long[] { 0, 0 }, new long[] { 1, 1 });

		final Cursor<FloatType> first = view.cursor();

		while (first.hasNext()) {
			first.fwd();
			Point cord = new Point(img.numDimensions());

			cord.setPosition(first);

			list.add(cord, first.get().copy());
			 System.out.println("Set of x co-ordinates Initial List : " +
			 cord.getDoublePosition(0));
			 System.out.println("Set of y co-ordinates Initial List : " +
			 cord.getDoublePosition(1));
			System.out.println("Values Initial list : " + first.get());

		}

		// list = getList(img);

	split(list, 0);
		// split(list, 1);

		Cursor<FloatType> testtwo = list.cursor();

		while (testtwo.hasNext()) {
			testtwo.fwd();
			Point newpoint = new Point(img.numDimensions());

			newpoint.setPosition(testtwo);

			 System.out.println("Set of x co-ordinates sorted List : " +
			 newpoint.getDoublePosition(0));
			 System.out.println("Set of y co-ordinates sorted List : " +
			 newpoint.getDoublePosition(1));
			System.out.println("Values sorted list : " + testtwo.get());

		}

	}
}
